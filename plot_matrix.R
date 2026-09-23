#!/usr/bin/env Rscript
# Heatmap of the cluster-ordered TM-score matrix written by tmalign_matrix.py.
#
#   Rscript plot_matrix.R [matrix.csv] [output.png]

library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else "tmalign_output/Clustered_TMalign_matrix.csv"
outfile <- if (length(args) >= 2) args[2] else "tm_score_heatmap.png"

# check.names = FALSE keeps the model filenames usable as column names.
tm_matrix <- as.matrix(read.csv(infile, header = TRUE, row.names = 1, check.names = FALSE))

# Reversed so that lower TM-scores are the lighter end of the palette.
palette <- rev(hcl.colors(100, palette = "Dark Mint"))

# Rows and columns are already ordered by HDBSCAN cluster, so pheatmap must not reorder.
pheatmap(tm_matrix,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE,
         color = palette, filename = outfile, width = 7, height = 6)

cat("Heatmap of", nrow(tm_matrix), "models written to", outfile, "\n")
