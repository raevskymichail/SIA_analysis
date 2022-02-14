#!/usr/bin/env Rscript

library("MutEffect")
source("MutEffect/R/preprocessing.R", chdir = TRUE)
source("MutEffect/R/visualization.R", chdir = TRUE)
library("optparse")

if (!require("pacman")) install.packages("pacman")
pacman::p_load("readr", "grid", "ggthemes", "scales", "FactoMineR", "factoextra", "ggplot2", "ggpubr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("preprocessCore")) BiocManager::install("preprocessCore")

# ## Usage
# ```bash
# cd ./compare_BT_vs_Norms
# Rscript ./scripts/plot_pca_gene_exp_by_source.R --exp-matrix="data/exp_data/processed/Oncobox_BT_exp_matrix.txt" \
#                              --output-dir="analyses/Oncobox_BT_vs_Norms/plots/pca_gene_exp_by_source" \
#                              --output-filename="pca_plot_by_source" \
#                              --samples-list="data/exp_data/Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
# Rscript ./scripts/plot_pca_gene_exp_by_source.R --exp-matrix="data/exp_data/processed/GTEx_BT_w_Oncobox_BT_and_norms.txt" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/pca_gene_exp_by_source" \
#                              --output-filename="pca_plot_by_source" \
#                              --samples-list="data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
# # FIX - Plot only 7 ANTE samples that pass QC
#
# Rscript ./scripts/plot_pca_gene_exp_by_source.R --exp-matrix="data/exp_data/processed/GTEx_BT_w_Oncobox_BT_and_norms.txt" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/pca_gene_exp_by_source" \
#                              --output-filename="pca_plot_by_source_passed_QC" \
#                              --samples-list="data/exp_data/all_samples_passed_QC_w_6_ANTE.txt"
# ```

option_list <- list(
  make_option(c("-e", "--exp-matrix"),
              dest = "exp_matrix",
              type = "character",
              metavar = "character",
              default = "data/exp_data/processed/exp_matrix.txt",
              help = "Path to processed expression matrix [default= %default]"),
  make_option(c("s-", "--samples-list"),
              dest = "samples_list",
              type = "character",
              metavar = "character",
              default = NULL,
              help = "Path to a list of TCGA sample codes to use in analysis"),
  make_option(c("-o", "--output-dir"),
              dest = "output_dir",
              type = "character",
              metavar = "character",
              help = "Path to output directory"),
  make_option(c("-f", "--output-filename"),
              dest = "output_filename",
              type = "character",
              metavar = "character",
              help = "Common name for output files"),
  make_option(c("-m", "--sample-annotations"),
              dest = "sample_annotations",
              type = "character",
              metavar = "character",
              default = "data/exp_data/exp_samples_annotation.csv",
              help = "Path to sample annotations database [default= %default]")
)

# opt <- list(exp_matrix = "data/exp_data/processed/Oncobox_BT_exp_matrix.txt",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv",
#             samples_list = "data/exp_data/all_samples_passed_QC.txt",
#             output_filename = "pca_plot_by_source")
# opt <- list(exp_matrix = "data/exp_data/processed/GTEx_BT_w_Oncobox_BT_and_norms.txt",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv",
#             samples_list = "data/exp_data/all_samples_passed_QC.txt",
#             output_dir = "analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/pca_gene_exp_by_source_passed_QC",
#             output_filename = "pca_plot_by_source")

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$output_dir)) {
  dir.create(file.path(getwd(), opt$output_dir), recursive = TRUE)
}

# Read and normaliza the expression matrix and select needed samples
exp_matrix <- load_exp_matrix(exp_matrix_path = opt$exp_matrix,
                              drop_dupl_genes = TRUE,
                              sample_list_path = opt$samples_list,
                              quantile_normalization = TRUE,
                              log10_scalling = TRUE)

# Obtain sample annotations
sample_ann <- get_exp_samples_annotations(sample_ann_table_path = opt$sample_annotations,
                                          sample_list_path = opt$samples_list)

filter_non_common_samples(exp_matrix, sample_ann)

exp_matrix <- t(exp_matrix[, names(sample_ann)])

# Run PCA - Gene Expression
pca_plot <- plot_pca(data = exp_matrix, metadata = sample_ann)

export_analysis_plot(filename = opt$output_filename,
                     plot = pca_plot,
                     path = opt$output_dir,
                     scale = 1,
                     width = 210 / 1.80,
                     height = 297 / 2.65,
                     units = "mm",
                     load_fonts = FALSE)
