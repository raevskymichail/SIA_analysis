#!/usr/bin/env Rscript

library("MutEffect")
source("MutEffect/R/preprocessing.R", chdir = TRUE)
source("MutEffect/R/visualization.R", chdir = TRUE)
library("optparse")

if (!require("pacman")) install.packages("pacman")
pacman::p_load("readr", "readxl", "grid", "ggthemes", "scales", "FactoMineR", "factoextra", "ggplot2", "ggpubr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("preprocessCore")) BiocManager::install("preprocessCore")

# ## Usage
# ```bash
# cd ./compare_BT_vs_Norms
# Rscript ./scripts/plot_pca_pal_by_source.R --pathways-activations-db="data/pathways_activation_scores/Oncobox_BT_vs_Norm.xlsx" \
#                              --output-dir="analyses/Oncobox_BT_vs_Norms/plots/pca_pal_by_source" \
#                              --output-filename="pca_plot_by_source" \
#                              --samples-list="data/exp_data/Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
# Rscript ./scripts/plot_pca_pal_by_source.R --pathways-activations-db="data/pathways_activation_scores/GTEx_BT_w_Oncobox_BT_and_norms_PAL.xlsx" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/pca_pal_by_source" \
#                              --output-filename="pca_plot_by_source" \
#                              --samples-list="data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
## Plot only significant pathways
# Rscript ./scripts/plot_pca_pal_by_source.R --pathways-activations-db="data/pathways_activation_scores/GTEx_BT_w_Oncobox_BT_and_norms_PAL_significant.xlsx" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/pca_pal_by_source/only_sign_pws" \
#                              --output-filename="pca_plot_by_source" \
#                              --samples-list="data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
## Plot only significant pathways with 10+ genes
# Rscript ./scripts/plot_pca_pal_by_source.R --pathways-activations-db="data/pathways_activation_scores/GTEx_BT_w_Oncobox_BT_and_norms_PAL_10_plus_genes.xlsx" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/pca_pal_by_source/only_sign_pws_10_plus_genes" \
#                              --output-filename="pca_plot_by_source" \
#                              --samples-list="data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
## Plot only significant pathways with 10+ genes + Plot only 7 ANTE samples that pass QC
# Rscript ./scripts/plot_pca_pal_by_source.R --pathways-activations-db="data/pathways_activation_scores/GTEx_BT_w_Oncobox_BT_and_norms_PAL_10_plus_genes.xlsx" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/pca_pal_by_source/only_sign_pws_10_plus_genes" \
#                              --output-filename="pca_plot_by_source_passed_QC" \
#                              --samples-list="data/exp_data/all_samples_passed_QC_w_6_ANTE.txt"
# ```

option_list <- list(
  make_option(c("-d", "--pathways-activations-db"),
              dest = "pathways_activations_db",
              type = "character",
              metavar = "character",
              help = "Path to sample annotations (drug scores) database [default= %default]"),
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
  make_option(c("-r", "--return-ranks"),
              dest = "return_ranks",
              action = "store_true",
              default = FALSE,
              help = "Flag indicating whether Drug Score / Pathway Activation values should be presented as ranks [default= %default]"),
  make_option(c("-m", "--sample-annotations"),
              dest = "sample_annotations",
              type = "character",
              metavar = "character",
              default = "data/exp_data/exp_samples_annotation.csv",
              help = "Path to sample annotations database [default= %default]")
)

# opt <- list(pathways_activations_db = "data/pathways_activation_scores/GTEx_BT_w_Oncobox_BT_and_norms_PAL.xlsx",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv",
#             output_dir = "analyses/Oncobox_BT_vs_Norms/plots/pca_pal_by_source",
#             output_filename = "pca_plot_by_source",
#             samples_list = "data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt",
#             return_ranks = FALSE)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$output_dir)) {
  dir.create(file.path(getwd(), opt$output_dir), recursive = TRUE)
}


# Obtain sample annotations (Pathways Activations)
pal_matrix <- get_pathways_activation_sample_annotations(pathways_activation_db_path = opt$pathways_activations_db,
                                                         return_ranks = opt$return_ranks,
                                                         return_dataframe = TRUE)

# Obtain sample annotations
sample_ann <- get_exp_samples_annotations(sample_ann_table_path = opt$sample_annotations,
                                          sample_list_path = opt$samples_list)

# Drop non-common samples
common_samples <- intersect(names(sample_ann), colnames(pal_matrix))
sample_ann <- sample_ann[common_samples]
pal_matrix <- pal_matrix[, common_samples]

# Run PCA - PAL Scores
pca_plot <- plot_pca(data = t(pal_matrix), metadata = sample_ann)

export_analysis_plot(filename = opt$output_filename,
                     plot = pca_plot,
                     path = opt$output_dir,
                     scale = 1,
                     width = 210 / 1.80,
                     height = 297 / 2.65,
                     units = "mm",
                     load_fonts = FALSE)
