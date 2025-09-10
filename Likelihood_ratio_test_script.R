# --- Load libraries ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "DEGreport"), ask = FALSE)

library(DESeq2)
library(DEGreport)
library(tidyverse)

# --- Import count matrix ---
data <- read.csv(
  "Counts.csv",
  header = TRUE,
  row.names = 1)

# Convert to numeric matrix
countData <- as.matrix(data)

# Define sample IDs
sample_names <- colnames(countData)

# --- Define metadata ---
Conditions <- c(
  "PMU","PMU","PMU","PMU","PMU","PMU","PMU","PMU",
  "PMU","PMU","PMU","PMU","PMU","PMU","PMU","PMU",  
  "UMU","UMU","UMU","UMU","UMU","UMU","UMU","UMU",
  "PMU","PMU","PMU","PMU","PMU","PMU","PMU",
  "UMU","UMU","UMU","UMU","UMU",
  "PMU","PMU","PMU","PMU","PMU","PMU","PMU","PMU",
  "UMU","UMU","UMU","UMU"
)

Time <- c(
  "Day1","Day2","Day4","Day28","Day1","Day28","Day4","Day2",
  "Day1","Day28","Day4","Day2","Day1","Day2","Day4","Day28",  
  "Day2","Day2","Day1","Day2","Day1","Day1","Day2","Day1",
  "Day4","Day4","Day4","Day4","Day4","Day4","Day4",
  "Day28","Day28","Day28","Day28","Day28",
  "Day4","Day4","Day4","Day4","Day4","Day4","Day4","Day4",
  "Day4","Day4","Day4","Day4"
)

Flowcell <- c(
  "7","7","7","7","7","7","7","7","7","7","7","7","7","7","7","7",
  "8","8","8","8","8","8","8","8","8","8","8","8","8","8","8",
  "9","9","9","9","9","9","9","9","9","9","9","9","9","9","9","9",
  "9"
)

# --- Build sample metadata ---
group.table <- data.frame(
  Condition = factor(Conditions, levels = c("PMU", "UMU")),
  Time      = factor(Time, levels = c("Day1", "Day2", "Day4", "Day28")),
  Flowcell  = factor(Flowcell),
  row.names = sample_names
)

# --- Create DESeq2 dataset ---
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = group.table,
  design    = ~ Flowcell + Condition * Time
)

# --- Fit model with likelihood ratio test (LRT) ---
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ Flowcell)

# Extract results
res_LRT <- results(dds_lrt)

# Number of genes tested
nrow(res_LRT)

# --- Filter significant genes ---
sig_res <- as.data.frame(res_LRT) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj) & padj < 0.001)

# List of significant gene IDs
sigLRT_genes <- sig_res$gene
length(sigLRT_genes)

# --- Sort by padj ---
resLRT_sorted <- sig_res[order(sig_res$padj), ]

# --- Normalised counts (VST) ---
vsd <- vst(dds_lrt, blind = FALSE)
vst_mat <- assay(vsd)

# Subset for significant genes only
vst_subset <- vst_mat[sigLRT_genes, ]

# --- Clustering with degPatterns ---
clusters <- degPatterns(
  vst_subset,
  metadata = group.table,
  time = "Time",
  col = "Condition",
  scale = TRUE,
  minc = 30)

# Extract clustering results
cluster_groups <- clusters$df

# Save each cluster to CSV
for (clust in unique(cluster_groups$cluster)) {
  group_df <- filter(cluster_groups, cluster == clust)
  write.csv(group_df, paste0("cluster", clust, "_newfull.csv"), row.names = FALSE)
}