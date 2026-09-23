#!/usr/bin/env python3
"""All-against-all TM-align pipeline: pTM filtering -> TM-align -> score matrix -> HDBSCAN ordering.

Runs end to end from a directory of ColabFold outputs to a cluster-ordered TM-score
matrix ready for plotting with plot_matrix.R.

    python tmalign_matrix.py --colabfold-dir ./colabfold_output

Every intermediate is written to --outdir so each stage remains inspectable, but no
stage requires manual editing of its input.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import hdbscan
import numpy as np
import pandas as pd

# One TM-align record. Chain names are captured as \S+ so the trailing
# "(to be superimposed onto Chain_2)" annotation is excluded, and \s* after each
# "TM-score=" absorbs the space that TM-align prints before the value.
RECORD_RE = re.compile(
    r"^Name of Chain_1:\s*(?P<chain_1>\S+).*?"
    r"^Name of Chain_2:\s*(?P<chain_2>\S+).*?"
    r"^TM-score=\s*(?P<tm_score_1>[\d.]+).*?LN=\s*(?P<length_1>\d+).*?"
    r"^TM-score=\s*(?P<tm_score_2>[\d.]+).*?LN=\s*(?P<length_2>\d+)",
    re.MULTILINE | re.DOTALL,
)


def collect_ptm(colabfold_dir: Path) -> pd.DataFrame:
    """Read the pTM score out of every ColabFold *_scores_*.json in a directory."""
    paths = sorted(colabfold_dir.glob("*.json"))
    if not paths:
        sys.exit(f"No .json files found in {colabfold_dir}")
    return pd.DataFrame(
        {"json_file": p.name, "ptm": json.loads(p.read_text())["ptm"]} for p in paths
    )


def select_models(ptm: pd.DataFrame, threshold: float, colabfold_dir: Path,
                  model_dir: Path) -> pd.DataFrame:
    """Copy the PDB of every model passing the pTM threshold into model_dir.

    ColabFold writes the score and structure of one model under matching names, e.g.
    job_scores_rank_001_....json / job_relaxed_rank_001_....pdb, so the PDB filename is
    recovered from the JSON filename by substitution.
    """
    selected = ptm.loc[ptm["ptm"] >= threshold].copy()
    selected["pdb_file"] = (selected["json_file"]
                            .str.replace("scores", "relaxed", regex=False)
                            .str.replace(".json", ".pdb", regex=False))

    missing = [f for f in selected["pdb_file"] if not (colabfold_dir / f).is_file()]
    if missing:
        sys.exit(f"{len(missing)} PDB file(s) named by a passing JSON are absent, "
                 f"first: {missing[0]}")

    model_dir.mkdir(parents=True, exist_ok=True)
    for f in selected["pdb_file"]:
        shutil.copy2(colabfold_dir / f, model_dir / f)
    return selected


def run_tmalign(model_dir: Path, chain_list: Path, results: Path, exe: str) -> None:
    """All-against-all alignment, equivalent to:

        ls <model_dir> > chain_list.txt
        TMalign -dir <model_dir>/ chain_list.txt > TMalign_results.txt
    """
    chain_list.write_text("".join(f"{p.name}\n" for p in sorted(model_dir.glob("*.pdb"))))
    with results.open("w") as out:
        subprocess.run([exe, "-dir", f"{model_dir}{os.sep}", str(chain_list)],
                       stdout=out, check=True)


def parse_tmalign(results: Path) -> pd.DataFrame:
    """Parse raw TM-align output into one row per aligned pair.

    Replaces the four-step text-munging chain of the original pipeline
    (Extract_TMscore_2Txt -> Replace_space -> Extract_TMScore_2CSV) with a single pass.
    """
    raw = results.read_text()
    records = [m.groupdict() for m in RECORD_RE.finditer(raw)]
    expected = raw.count("Name of Chain_1:")
    if len(records) != expected:
        sys.exit(f"Parsed {len(records)} of {expected} alignment records from {results}")

    columns = ["chain_1", "chain_2", "tm_score_1", "tm_score_2", "length_1", "length_2"]
    pairs = pd.DataFrame(records)[columns]
    pairs[["chain_1", "chain_2"]] = pairs[["chain_1", "chain_2"]].map(os.path.basename)
    return pairs.astype({"tm_score_1": float, "tm_score_2": float,
                         "length_1": int, "length_2": int})


def build_matrix(pairs: pd.DataFrame) -> pd.DataFrame:
    """Square, symmetric matrix of max(TM-score_1, TM-score_2) per pair.

    TM-align reports two scores per alignment because each is normalised by the length of
    a different chain; taking the maximum is symmetric in the two chains, so the (i, j)
    and (j, i) runs collapse to the same value. TM-align never aligns a chain against
    itself, so the diagonal is set here to its exact value of 1.0.
    """
    pairs = pairs.assign(tm_score=pairs[["tm_score_1", "tm_score_2"]].max(axis=1))
    names = sorted(set(pairs["chain_1"]) | set(pairs["chain_2"]))
    matrix = (pairs.pivot_table(index="chain_2", columns="chain_1", values="tm_score")
                   .reindex(index=names, columns=names))

    values = matrix.to_numpy(copy=True)
    np.fill_diagonal(values, 1.0)
    matrix = pd.DataFrame(values, index=names, columns=names)

    if matrix.isna().to_numpy().any():
        n = int(matrix.isna().to_numpy().sum())
        sys.exit(f"{n} pair(s) absent from the TM-align output; matrix is incomplete")
    return matrix


def cluster_and_order(matrix: pd.DataFrame, min_cluster_size: int,
                      drop_noise: bool) -> tuple[pd.Series, pd.DataFrame]:
    """Cluster models on their matrix rows and reorder both axes by cluster label.

    HDBSCAN is fitted on the full matrix, including the models it later labels as noise,
    then rows and columns are sorted by label. The sort is stable, so models within a
    cluster keep their alphabetical order and the diagonal stays intact.
    """
    labels = pd.Series(
        hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=None,
                        metric="euclidean").fit_predict(matrix.to_numpy()),
        index=matrix.index, name="hdbscan_label",
    )
    keep = labels[labels != -1] if drop_noise else labels
    order = keep.sort_values(kind="stable").index
    return labels, matrix.loc[order, order]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--colabfold-dir", type=Path, required=True,
                   help="directory of ColabFold *_scores_*.json and *_relaxed_*.pdb files")
    p.add_argument("--outdir", type=Path, default=Path("tmalign_output"))
    p.add_argument("--ptm-threshold", type=float, default=0.70,
                   help="minimum pTM to carry a model forward (default: %(default)s)")
    p.add_argument("--min-cluster-size", type=int, default=5,
                   help="HDBSCAN min_cluster_size (default: %(default)s)")
    p.add_argument("--keep-noise", action="store_true",
                   help="retain models HDBSCAN labels as noise (-1) in the plotted matrix")
    p.add_argument("--tmalign", default="TMalign", help="TM-align executable")
    p.add_argument("--reuse-alignments", action="store_true",
                   help="skip the TM-align run if its output is already present")
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    model_dir = args.outdir / "selected_models"
    chain_list = args.outdir / "chain_list.txt"
    results = args.outdir / "TMalign_results.txt"

    ptm = collect_ptm(args.colabfold_dir)
    ptm.to_csv(args.outdir / "ptm_scores.csv", index=False)

    selected = select_models(ptm, args.ptm_threshold, args.colabfold_dir, model_dir)
    selected.to_csv(args.outdir / "selected_models.csv", index=False)
    print(f"{len(selected)} of {len(ptm)} models pass pTM >= {args.ptm_threshold}")

    if not (args.reuse_alignments and results.is_file()):
        run_tmalign(model_dir, chain_list, results, args.tmalign)

    pairs = parse_tmalign(results)
    pairs.to_csv(args.outdir / "tmalign_pairwise.csv", index=False)

    matrix = build_matrix(pairs)
    matrix.to_csv(args.outdir / "tm_score_matrix.csv")

    labels, clustered = cluster_and_order(matrix, args.min_cluster_size,
                                          drop_noise=not args.keep_noise)
    labels.to_csv(args.outdir / "hdbscan_labels.csv")
    clustered.to_csv(args.outdir / "Clustered_TMalign_matrix.csv")

    counts = labels.value_counts()
    print(f"{len(matrix)} models, {(counts.index != -1).sum()} clusters, "
          f"{int(counts.get(-1, 0))} noise points")
    print(f"Plot with: Rscript plot_matrix.R {args.outdir / 'Clustered_TMalign_matrix.csv'}")


if __name__ == "__main__":
    main()
