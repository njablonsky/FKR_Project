library(tidyverse)
library(ggplot2)
library(pheatmap)

# 1. Load data
dat <- read_csv("KRF_orthogroups.csv")

# Expect columns: OG, subfamily, Genome, Family
glimpse(dat)


# 2. Build family map 

fam_map <- dat %>%
  distinct(Genome, Family)

# Each genome must belong to exactly one family
dup_fams <- fam_map %>%
  count(Genome) %>%
  filter(n > 1)

if(nrow(dup_fams) > 0){
  stop("Some genomes have multiple family labels!")
}

all_genomes <- fam_map$Genome

# 3. Remove problematic OGs (paralog-rich)

copy_counts <- dat %>%
  count(subfamily, Genome)

bad_ogs <- copy_counts %>%
  group_by(subfamily) %>%
  summarise(
    max_copy = max(n),
    total = sum(n),
    .groups = "drop"
  ) %>%
  filter(max_copy > 3 | total > 100) %>%
  pull(subfamily)

dat_filt <- dat %>%
  filter(!subfamily %in% bad_ogs)

cat("OGs before:", length(unique(dat$subfamily)), "\n")
cat("OGs after :", length(unique(dat_filt$subfamily)), "\n")

# 4. Collapse to presence/absence
# One row = OG × Genome
collapsed <- dat_filt %>%
  distinct(subfamily, Genome)

# 5. Build complete presence/absence matrix

mat <- collapsed %>%
  mutate(present = 1) %>%
  complete(
    subfamily,
    Genome = all_genomes,
    fill = list(present = 0)
  ) %>%
  pivot_wider(
    names_from = Genome,
    values_from = present
  )

# 6. Sanity check: no missing genomes

missing1 <- setdiff(all_genomes, colnames(mat))
missing2 <- setdiff(colnames(mat), all_genomes)

if(length(missing1) > 0 | length(missing2) > 0){
  stop("Genome name mismatch between fam_map and matrix!")
}

missing1
missing2
# 7. Fisher test function

run_fisher <- function(mat, fam_map, target){
  
  target_genomes <- fam_map %>%
    filter(Family == target) %>%
    pull(Genome)
  
  other_genomes <- fam_map %>%
    filter(Family != target) %>%
    pull(Genome)
  
  nA <- length(target_genomes)
  nB <- length(other_genomes)
  
  res <- mat %>%
    rowwise() %>%
    mutate(
      
      # Counts
      a = sum(c_across(all_of(target_genomes))),
      b = sum(c_across(all_of(other_genomes))),
      
      na = nA - a,
      nb = nB - b,
      
      # Fisher test (run once)
      ft = list(
        fisher.test(matrix(c(a, na, b, nb), nrow = 2))
      ),
      
      p  = ft$p.value,
      OR = ft$estimate
      
    ) %>%
    ungroup() %>%
    select(-ft) %>%
    mutate(
      
      FDR = p.adjust(p, method = "BH"),
      Target = target,
      
      freq_target = a / nA,
      freq_other  = b / nB
    )
  
  return(res)
}

# 8. Run enrichment for all families

groups <- sort(unique(fam_map$Family))

results <- map_dfr(
  groups,
  ~run_fisher(mat, fam_map, .x)
)

# 9. Extract high-confidence hits

hits <- results %>%
  filter(
    FDR < 0.05,
    freq_target >= 0.8,
    freq_other <= 0.2
  ) %>%
  arrange(FDR)

# 10. Save output

write_csv(results, "OG_enrichment_all.csv")
write_csv(hits,    "OG_enrichment_hits.csv")

# 11. Plot results

top_ogs <- hits %>%
  arrange(FDR) %>%
  slice_head(n = 40) %>%
  pull(subfamily)

t_mat <- transpose(mat)
t_mat <- as.matrix(t_mat)

heat_dat <- mat %>%
  filter(subfamily %in% top_ogs) %>%
  column_to_rownames("subfamily")

pheatmap(
  heat_dat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  border_color = "grey85",
  color = c("white","firebrick4"),
  legend = FALSE,
  fontsize_row = 8,
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize_col = 6,
)
