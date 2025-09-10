# MSc-Bioinformatics_project_BIOLM0034

This repository contains R scripts for performing RNA-seq data analysis, including normalisation, exploratory visualisation, differential expression, likelihood ratio tests, clustering, and GO enrichment.

1. PCA Analysis (PCA script.R)

Input: Normalised variance-stabilized counts (vst).
Purpose: Visualise sample clustering based on conditions, time points, and batch/flowcell effects.
Output: PCA plots with colours for conditions and shapes for batches.

2. Differential Expression Analysis (Differential_expression_analysis_script.R)

Input: Raw count matrix (Counts.csv) and metadata.
Purpose: Run pairwise contrasts between conditions/timepoints using DESeq2.
Output:
-Differential expression results for each contrast.
-Lists of significant genes (padj < 0.01).

3. Likelihood Ratio Test (LRT) (Likelihood_ratio_test_script.R)

Input: DESeq2 dataset.
Purpose: Identify genes with significant changes across conditions/timepoints (global test).
Output:
-Significant gene list.
-Clustering of expression patterns across samples using DEGreport::degPatterns.
-Cluster gene lists exported to .csv.

4. GO Enrichment Barplot (Barplot_script.R)

Input: Pre-computed GO enrichment results (GO.csv).
Purpose: Plot the top 20 enriched GO terms ranked by –log10(p-value).
Output: Publication-quality barplots.
