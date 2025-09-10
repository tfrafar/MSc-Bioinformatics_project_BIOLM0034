
# Install BiocManager if missing
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Load required libraries
library(vroom)
library(tidyverse)
library(ggplot2)
library(DESeq2)

# ================================
# 1. Load and prepare count data
# ================================

# Ensure rownames = gene IDs and colnames = sample IDs
data <- read.csv("Counts.csv", header = TRUE, row.names = 1)

# Convert to numeric matrix for DESeq2
countData <- as.matrix(data)

# Quick preview of the data
head(data)

# ================================
# 2. Define metadata
# ================================
# Experimental design variables:
# - Flowcell: sequencing batch information
# - Batches: alternative batch info (can replace Flowcell if needed)
# - Conditions_Time: experimental conditions (treatment + timepoints)

Batches <- c("SS8","SS8","SS8","SS8","SS8","SS8","SS9","SS9","SS9","SS9","SS14","SS14","SS14",
             "SS14","SS14","SS14","SS14","SS14","SS14","SS14","SS14","SS14","SS14",
             "SS14","SS14","SS14","SS14","SS14","SS14","SS15","SS15","SS15","SS15","SS15","SS15","SS15",
             "SS15","SS16","SS16","SS16","SS16","SS16","SS16",
             "SS16","SS16","SS16","SS16","SS18","SS18","SS18","SS18","SS18","SS18",
             "SS18","SS18","SS18","SS18","SS18","SS18","SS19","SS19","SS19",
             "SS19","SS19","SS19","SS19","SS19","SS19","SS19","SS20","SS20",
             "SS20","SS20","SS20","SS20","SS21","SS21")

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

# Build colData (sample metadata for DESeq2)
# NOTE: Replace Flowcell with Batches if you want batch correction by Batches instead
colData <- data.frame(Conditions_Time, Flowcell)
rownames(colData) <- colnames(countData)

# ================================
# 3. Create DESeq2 dataset & run DESeq
# ================================
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ Flowcell + Conditions_Time   # or ~ Batches + Conditions_Time
)

des <- DESeq(dds)

# ================================
# 4. Variance stabilizing transformation (VST)
# ================================
vsdata <- vst(des, blind = FALSE)

# Extract PCA data (first two PCs)
pcaData <- plotPCA(vsdata, intgroup = "Conditions_Time", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# ================================
# 5. Define condition names and manual colors
# ================================
condition_names <- c(
  "Naive",
  "PMUDay1", "PMUDay2", "PMUDay4","PMUDay7","PMUDay28", "PMUDay56",
  "UMUDay1", "UMUDay2", "UMUDay4", "UMUDay28"
)

manual_colors <- c(
  "grey",       # Naive
  "#87CEEB",    # PMUDay1
  "pink",       # PMUDay2
  "purple",     # PMUDay4
  "cyan",       # PMUDay7
  "lightgreen", # PMUDay28
  "#6B8E23",    # PMUDay56
  "#F0E442",    # UMUDay1
  "tomato",     # UMUDay2
  "deeppink2",  # UMUDay4
  "brown"       # UMUDay28
)
names(manual_colors) <- condition_names

# Force factor order for legend
pcaData$Conditions_Time <- factor(as.character(pcaData$Conditions_Time), levels = condition_names)

# ================================
# 6. Shapes for Flowcell (or Batches)
# ================================
pcaData$Flowcell <- factor(as.character(pcaData$Flowcell))
flowcell_levels <- levels(pcaData$Flowcell)

manual_shapes <- c(21, 22, 23, 24)[seq_along(flowcell_levels)]
names(manual_shapes) <- flowcell_levels

# To use Batches instead, swap with:
# pcaData$Batches <- factor(as.character(pcaData$Batches))
# batch_levels <- levels(pcaData$Batches)
# manual_shapes <- c(21, 22, 23, 24, 25, 26, 27, 28, 29)[seq_along(batch_levels)]
# names(manual_shapes) <- batch_levels

# ================================
# 7. PCA plots
# ================================

# PCA plot with both Condition and Flowcell
ggplot(pcaData, aes(x = PC1, y = PC2, fill = Conditions_Time, shape = Flowcell)) +
  geom_point(size = 4, color = "black") +
  scale_fill_manual(values = manual_colors, na.value = "grey50") +
  scale_shape_manual(values = manual_shapes) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of RNA-seq Samples by Condition and Flowcell",
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)"),
    fill = "Condition",
    shape = "Flowcell" # change to "Batch" if using Batches
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4, shape = 21)),
    shape = guide_legend(override.aes = list(size = 4))
  ) +
  theme(
    panel.background = element_rect(fill = "#F5F5F5", color = NA),
    panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )

# Simple PCA with just conditions
ggplot(pcaData, aes(x = PC1, y = PC2, color = Conditions_Time)) +
  geom_point(size = 3) +
  scale_color_manual(values = manual_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of RNA-seq Samples by Condition (color only)",
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)"),
    color = "Condition"
  ) +
  theme(
    panel.background = element_rect(fill = "#F5F5F5", color = NA),
    panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right")
  