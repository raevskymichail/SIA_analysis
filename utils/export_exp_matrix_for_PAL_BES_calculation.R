#!/usr/bin/env Rscript

library("MutEffect")
source("MutEffect/R/preprocessing.R", chdir = TRUE)
library("dplyr")
library("readr")
library("preprocessCore")
library("psych")
library("optparse")

# ## Usage
# ```bash
# cd ./compare_BT_vs_Norms
# Rscript ./utils/export_exp_matrix_for_PAL_BES_calculation.R --exp-matrix="data/exp_data/processed/Oncobox_BT_exp_matrix.txt" \
#                                                             --samples-list="data/exp_data/Oncobox_BT_and_Norms_sample_names.txt" \
#                                                             --output-matrix="data/exp_data/processed/Oncobox_BT_exp_matrix_for_PAL_BES.tsv"
#
# Rscript ./utils/export_exp_matrix_for_PAL_BES_calculation.R --exp-matrix="data/exp_data/processed/GTEx_BT_and_Oncobox_BT_norm_exp_matrix.txt" \
#                                                             --samples-list="data/exp_data/GTEx_BT_and_Oncobox_BT_norm_sample_names.txt" \
#                                                             --output-matrix="data/exp_data/processed/GTEx_BT_and_Oncobox_BT_norm_exp_matrix_for_PAL_BES.tsv"
#
# Rscript ./utils/export_exp_matrix_for_PAL_BES_calculation.R --exp-matrix="data/exp_data/processed/GTEx_BT_w_Oncobox_BT_and_norms.txt" \
#                                                             --samples-list="data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt" \
#                                                             --output-matrix="data/exp_data/processed/GTEx_BT_w_Oncobox_BT_and_norms_exp_matrix_for_PAL_BES.tsv"
#
# Rscript ./utils/export_exp_matrix_for_PAL_BES_calculation.R --exp-matrix="data/exp_data/processed/GTEx_BT_w_Oncobox_BT_and_norms.txt" \
#                                                             --samples-list="data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt" \
#                                                             --output-matrix="data/exp_data/processed/Oncobox_BT_and_GTEx_norms_exp_matrix_for_PAL_BES.tsv"
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
  make_option(c("-o", "--output-matrix"),
              dest = "output_matrix",
              type = "character",
              metavar = "character",
              help = "Path to output an expression matrix"),
  make_option(c("-m", "--sample-annotations"),
              dest = "sample_annotations",
              type = "character",
              metavar = "character",
              default = "data/exp_data/exp_samples_annotation.csv",
              help = "Path to sample annotations database [default= %default]")
)

# opt <- list(exp_matrix = "data/exp_data/processed/Oncobox_BT_exp_matrix.txt",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv")
# opt <- list(exp_matrix = "data/exp_data/processed/GTEx_BT_and_Oncobox_BT_norm_exp_matrix.txt",
#             samples_list = "data/exp_data/GTEx_BT_and_Oncobox_BT_norm_sample_names.txt",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv")

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Read and normaliza the expression matrix and select needed samples
exp_matrix <- load_exp_matrix(exp_matrix_path = opt$exp_matrix,
                              drop_dupl_genes = TRUE,
                              sample_list_path = opt$samples_list,
                              quantile_normalization = TRUE,
                              log10_scalling = FALSE)

# Obtain sample annotations
sample_ann <- get_exp_samples_annotations(sample_ann_table_path = opt$sample_annotations,
                                          sample_list_path = opt$samples_list)

filter_non_common_samples(exp_matrix, sample_ann)


# Calculate pseudo samples (geom. mean)
# control_samples <- names(sample_ann)[sample_ann == "ANTE normal brain"]
control_samples <- names(sample_ann)[sample_ann != "SIA"]
case_samples <- names(sample_ann)[sample_ann == "SIA"]
# case_samples <- names(sample_ann)[sample_ann == "GTEx spinal cord"]
# case_samples <- names(sample_ann)[sample_ann != "ANTE normal brain"]

exp_matrix_control <- exp_matrix[colnames(exp_matrix) %in% control_samples]
exp_matrix_case <- exp_matrix[colnames(exp_matrix) %in% case_samples]

pseudo_control_sample <- geometric.mean(t(exp_matrix_control), na.rm = TRUE)
pseudo_control_sample <- as.data.frame(pseudo_control_sample)

pseudo_case_sample <- geometric.mean(t(exp_matrix_case), na.rm = TRUE)
pseudo_case_sample <- as.data.frame(pseudo_case_sample)

colnames(pseudo_control_sample) <- c("Norm_geomean")
colnames(pseudo_case_sample) <- c("Case_geomean")
exp_matrix <- cbind(exp_matrix, pseudo_control_sample)
exp_matrix <- cbind(exp_matrix, pseudo_case_sample)

# Format columns
colnames(exp_matrix)[colnames(exp_matrix) %in% control_samples] <- paste0("Norm_", control_samples)
colnames(exp_matrix)[colnames(exp_matrix) %in% case_samples] <- paste0("Case_", case_samples)

exp_matrix <- exp_matrix[complete.cases(exp_matrix),] # drop NAs
exp_matrix$SYMBOL <- rownames(exp_matrix)
exp_matrix <- exp_matrix[, c("SYMBOL", paste0("Norm_", control_samples), paste0("Case_", case_samples), "Norm_geomean", "Case_geomean")]


# Export prepared expression matrix
write.table(exp_matrix, opt$output_matrix, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
