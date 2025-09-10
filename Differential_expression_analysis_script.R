# ============================
# Differential Expression Analysis with DESeq2
# Includes contrasts, significant gene extraction,
# variance stabilization, and heatmap visualisation
# ============================

# --- Load required libraries ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

library(DESeq2)      # Core differential expression
library(pheatmap)    # Heatmap visualisation
library(tidyverse)   # Data manipulation

# --- Import count matrix ---
# Count matrix should have genes in rows and samples in columns
# "Counts.csv" should be prepared accordingly
data <- read.csv(
  "Counts.csv",
  header = TRUE,
  row.names = 1)

# Convert to numeric matrix
countData <- as.matrix(data)

# Extract sample IDs (column names = samples)
sample_names <- colnames(countData)

# --- Define metadata (sample annotation) ---
# "Flowcell" = optional technical covariate if needed
# "Conditions_Time" = biological condition + timepoint

Flowcell <- c("6","6","6","6","6","6","6","6","6","6",
              "7","7","7","7","7","7","7","7","7","7",
              "7","7","7","7","7","7","7","7","7","8","8","8","8","8","8","8",
              "8","8","8","8","8","8","8","8","8","8","8","9","9","9","9",
              "9","9","9","9","9","9","9","9","9","9","9","9",
              "9","9","9","9","9","9","9","9","9","9","9","9","9","9")

Conditions_Time <- c("Naive","Naive","Naive","Naive","Naive","Naive","Naive","Naive","Naive","Naive","Naive",
                     "Naive","Naive","PMUDay1","PMUDay1","PMUDay2","PMUDay4","PMUDay28","PMUDay28","PMUDay4",
                     "PMUDay2","PMUDay1","PMUDay28","PMUDay4","PMUDay2","PMUDay1","PMUDay2","PMUDay4",
                     "PMUDay28","UMUDay2","UMUDay2","UMUDay1","UMUDay2","UMUDay1","UMUDay1","UMUDay2",
                     "UMUDay1","PMUDay4","Naive","PMUDay4","PMUDay4","PMUDay4","Naive","PMUDay4","PMUDay4",
                     "Naive","PMUDay4","Naive","PMUDay56","PMUDay56","PMUDay56","PMUDay56",
                     "PMUDay56","UMUDay28","UMUDay28","UMUDay28","Naive","UMUDay28",
                     "UMUDay28","PMUDay4","PMUDay4","Naive","PMUDay4","PMUDay4","PMUDay4","PMUDay4","PMUDay4",
                     "PMUDay4","Naive","UMUDay4","UMUDay4","Naive","UMUDay4","Naive",
                     "UMUDay4","Naive","Naive")

# Build colData (metadata required by DESeq2)
colData <- data.frame(
  Conditions_Time = Conditions_Time,
  Flowcell = Flowcell
)
rownames(colData) <- colnames(countData)


# --- Create DESeq2 dataset ---
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ Flowcell + Conditions_Time   # account for batch (Flowcell) + condition/time
)

# --- Run DESeq2 pipeline ---
des <- DESeq(dds)

# --- Define contrasts (pairwise comparisons) ---
# PMU Day1 vs UMU Day1
res_day1_pmu_umu <- results(des, contrast = c("Conditions_Time", "PMUDay1", "UMUDay1"))

# PMU Day2 vs uPMU Day2
res_day2_pmu_umu <- results(des, contrast = c("Conditions_Time", "PMUDay2", "UMUDay2"))

# PMU Day4 vs uPMU Day4
res_day4_pmu_umu <- results(des, contrast = c("Conditions_Time", "PMUDay4", "UMUDay4"))

# PMU Day28 vs uPMU Day28
res_day28_pmu_umu <- results(des, contrast = c("Conditions_Time", "PMUDay28", "UMUDay28"))

# Naive vs PMU Day1
res_day1_pmu_naive <- results(des, contrast = c("Conditions_Time","PMUDay1","Naive"))

# Naive vs PMU Day2
res_day2_pmu_naive <- results(des, contrast = c("Conditions_Time","PMUDay2","Naive"))

# Naive vs PMU Day4
res_day4_pmu_naive <- results(des, contrast = c("Conditions_Time","PMUDay4","Naive"))

# Naive vs PMU Day28
res_day28_pmu_naive <- results(des, contrast = c("Conditions_Time","PMUDay28","Naive"))

# Naive vs uPMU Day1
res_day1_umu_naive <- results(des, contrast = c( "Conditions_Time","UMUDay1","Naive"))

# Naive vs uPMU Day2
res_day2_umu_naive <- results(des, contrast = c( "Conditions_Time","UMUDay2","Naive"))

# Naive vs uPMU Day4=default
res_day4_umu_naive <- results(des, contrast=c("Conditions_Time","UMUDay4","Naive"))

# Naive vs uPMU Day28
res_day28_umu_naive <- results(des, contrast = c("Conditions_Time","UMUDay28","Naive"))

# Naive vs PMU Day56
res_day56_pmu_naive <- results(des, contrast = c("Conditions_Time","PMUDay56","Naive"))

# Within condition , across time
# PMU1 vs PMU2
res_day1_day2_pmu <- results(des, contrast = c("Conditions_Time", "PMUDay2", "PMUDay1"))

# PMU1 vs PMU4
res_day1_day4_pmu <- results(des, contrast = c("Conditions_Time", "PMUDay4", "PMUDay1"))

# PMU1 vs PMU28
res_day1_day28_pmu <- results(des, contrast = c("Conditions_Time", "PMUDay28", "PMUDay1"))

# PMU2 vs PMU4
res_day2_day4_pmu <- results(des, contrast = c("Conditions_Time", "PMUDay4", "PMUDay2"))

# PMU2 vs PMU28
res_day2_day28_pmu <- results(des, contrast = c("Conditions_Time", "PMUDay28", "PMUDay2"))

# PMU4 vs PMU28
res_day4_day28_pmu <- results(des, contrast = c("Conditions_Time", "PMUDay28", "PMUDay4"))

# PMU1 vs PMU56
res_day1_day56_pmu <- results(des, contrast = c("Conditions_Time", "PMUDay56", "PMUDay1"))

# PMU2 vs PMU56
res_day2_day56_pmu <- results(des, contrast = c("Conditions_Time", "PMUDay56", "PMUDay2"))

# PMU4 vs PMU 56
res_day4_day56_pmu <- results(des, contrast = c("Conditions_Time", "PMUDay56", "PMUDay4"))

# PMU28 vs PMU56
res_day28_day56_pmu <- results(des, contrast = c("Conditions_Time", "PMUDay56", "PMUDay28"))

# UMU1 vs UMU2
res_day1_day2_umu <- results(des, contrast = c("Conditions_Time", "UMUDay2", "UMUDay1"))

# UMU2 vs UMU4
res_day2_day4_umu <- results(des, contrast = c("Conditions_Time", "UMUDay4", "UMUDay2"))

# UMU4 vs UMU28
res_day4_day28_umu <- results(des, contrast = c("Conditions_Time", "UMUDay28", "UMUDay4"))

# UMU2 vs UMU28
res_day2_day28_umu <- results(des, contrast = c("Conditions_Time", "UMUDay28", "UMUDay2"))

# UMU1 vs UMU4
res_day1_day4_umu <- results(des, contrast = c("Conditions_Time", "UMUDay4", "UMUDay1"))

# UMU1 vs UMU28
res_day1_day28_umu <- results(des, contrast = c("Conditions_Time", "UMUDay28", "UMUDay1"))

# --- Extract significant genes (padj < 0.01) ---
sig_genes_day1 <- rownames(subset(res_day1_pmu_umu, padj < 0.01 & !is.na(padj)))
sig_genes_day2 <- rownames(subset(res_day2_pmu_umu, padj < 0.01 & !is.na(padj)))
sig_genes_day4 <- rownames(subset(res_day4_pmu_umu, padj < 0.01 & !is.na(padj)))
sig_genes_day28 <- rownames(subset(res_day28_pmu_umu, padj < 0.01 & !is.na(padj)))
sig_genes_day1pmun <- rownames(subset(res_day1_pmu_naive, padj < 0.01 & !is.na(padj)))
sig_genes_day2pmun <- rownames(subset(res_day2_pmu_naive, padj < 0.01 & !is.na(padj)))
sig_genes_day4pmun <- rownames(subset(res_day4_pmu_naive, padj < 0.01 & !is.na(padj)))
sig_genes_day28pmun <- rownames(subset(res_day28_pmu_naive, padj < 0.01 & !is.na(padj)))
sig_genes_day1umun <- rownames(subset(res_day1_umu_naive, padj < 0.01 & !is.na(padj)))
sig_genes_day2umun <- rownames(subset(res_day2_umu_naive, padj < 0.01 & !is.na(padj)))
sig_genes_day4umun <- rownames(subset(res_day4_umu_naive, padj < 0.01 & !is.na(padj)))
sig_genes_day28umun <- rownames(subset(res_day28_umu_naive, padj < 0.01 & !is.na(padj)))
sig_genes_day56pmun <- rownames(subset(res_day56_pmu_naive, padj < 0.01 & !is.na(padj)))
##
sig_genes_day1pmu2 <- rownames(subset(res_day1_day2_pmu, padj < 0.01 & !is.na(padj)))
sig_genes_day1pmu4 <- rownames(subset(res_day1_day4_pmu, padj < 0.01 & !is.na(padj)))
sig_genes_day1pmu28 <- rownames(subset(res_day1_day28_pmu, padj < 0.01 & !is.na(padj)))
sig_genes_day2pmu4 <- rownames(subset(res_day2_day4_pmu, padj < 0.01 & !is.na(padj)))
sig_genes_day2pmu28 <- rownames(subset(res_day2_day28_pmu, padj < 0.01 & !is.na(padj)))
sig_genes_day4pmu28 <- rownames(subset(res_day4_day28_pmu, padj < 0.01 & !is.na(padj)))
sig_genes_day1pmu56 <- rownames(subset(res_day1_day56_pmu, padj < 0.01 & !is.na(padj)))
sig_genes_day2pmu56 <- rownames(subset(res_day2_day56_pmu, padj < 0.01 & !is.na(padj)))
sig_genes_day4pmu56 <- rownames(subset(res_day4_day56_pmu, padj < 0.01 & !is.na(padj)))
sig_genes_day28pmu56 <- rownames(subset(res_day28_day56_pmu, padj < 0.01 & !is.na(padj)))
##
sig_genes_day1umu28 <- rownames(subset(res_day1_day28_umu, padj < 0.01 & !is.na(padj)))
sig_genes_day1umu2 <- rownames(subset(res_day1_day2_umu, padj < 0.01 & !is.na(padj)))
sig_genes_day1umu4 <- rownames(subset(res_day1_day4_umu, padj < 0.01 & !is.na(padj)))
sig_genes_day2umu4 <- rownames(subset(res_day2_day4_umu, padj < 0.01 & !is.na(padj)))
sig_genes_day2umu28 <- rownames(subset(res_day2_day28_umu, padj < 0.01 & !is.na(padj)))
sig_genes_day4umu28 <- rownames(subset(res_day4_day28_umu, padj < 0.01 & !is.na(padj)))

# Combine into one list of DEGs
gene_list <- c(sig_genes_day1, sig_genes_day2,sig_genes_day4,sig_genes_day28,sig_genes_day1pmun,
               sig_genes_day2pmun,sig_genes_day4pmun,sig_genes_day28pmun,sig_genes_day56pmun,
               sig_genes_day1umun,sig_genes_day2umun, sig_genes_day4umun,sig_genes_day28umun,
               sig_genes_day1pmu2,sig_genes_day1pmu4,sig_genes_day1pmu28,sig_genes_day2pmu4,
               sig_genes_day2pmu28,sig_genes_day4pmu28,sig_genes_day1pmu56,sig_genes_day2pmu56,
               sig_genes_day4pmu56,sig_genes_day28pmu56,sig_genes_day1umu28,sig_genes_day1umu2,
               sig_genes_day1umu4,sig_genes_day2umu4,sig_genes_day2umu28,sig_genes_day4umu28)

# Total number of unique DEGs
length(gene_list)
# [1] 14765

# --- Variance Stabilizing Transformation (VST) ---
vsd <- vst(des, blind = FALSE)
norm_mat <- assay(vsd)   # normalized expression matrix

# Subset for DEGs
heatmap_mat <- norm_mat[gene_list, ]

# --- Heatmap visualization ---
# Define custom annotation colors
my_colors <- list(
  Condition = c(
    "Naive" = "grey",  
    "PMUDay1" = "#87CEEB",  
    "PMUDay2" = "pink",  
    "PMUDay4" = "purple",  
    "PMUDay28" = "lightgreen",  
    "PMUDay56" = "#6B8E23",
    "UMUDay1" = "#F0E442",  
    "UMUDay2" = "tomato", 
    "UMUDay4" = "deeppink2",  
    "UMUDay28" = "brown"   
    
  ))

# Create annotation for samples (which become rows after transpose)
annotation_row <- colData[, "Conditions_Time", drop = FALSE]
colnames(annotation_row) <- "Condition"
rownames(annotation_row) <- colnames(heatmap_mat)  # match to sample names

# Convert to factor with proper levels
annotation_row$Condition <- factor(annotation_row$Condition, levels = names(my_colors$Condition))

# Heatmap of DEGs
pheat<-pheatmap(t(heatmap_mat),  # 👈 Transpose to flip axes
                scale = "column",   # Scale across genes now (columns)
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                annotation_row = annotation_row,# Annotation goes with rows now
                annotation_colors = my_colors,
                show_rownames = FALSE,             # hide sample names 
                show_colnames = F,            # Hide gene names if too many
                fontsize = 10,
                color = colorRampPalette(c(" blue ", " white ", " red "))(20))

# Manual cluster analysis using cutree()
col_clusters <- cutree(pheat$tree_col, k = 2) 
genes_cluster1 <- names(col_clusters[col_clusters == 1])
cluster1_matrix <- heatmap_mat[genes_cluster1, ]

# Sub-heatmap for cluster 1
subpheat<-pheatmap(t(cluster1_matrix),
                   scale = "column",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   annotation_row = annotation_row,
                   annotation_colors = my_colors,
                   show_rownames = FALSE,
                   show_colnames = F,
                   color = colorRampPalette(c("blue", "white", "red"))(20))

# Export gene list for cluster 1
write.csv(data.frame(Gene = genes_cluster1),
          file = "cluster1_genes.csv",
          row.names = FALSE)