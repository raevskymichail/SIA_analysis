#!/usr/bin/env /usr/bin/Rscript

library("MutEffect")
source("MutEffect/R/preprocessing.R", chdir = TRUE)
source("MutEffect/R/visualization.R", chdir = TRUE)
library("DESeq2")
library("EnhancedVolcano")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("gridExtra")
library("grid")
library("readr")
library("readxl")
library("optparse")

# ## Usage
# ```bash
# cd ./compare_mut_gene_exp
# Rscript ./scripts/plot_dif_exp_analysis.R --exp-matrix="data/exp_data/processed/Oncobox_BT_exp_matrix.txt" \
#                              --output-dir="analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano" \
#                              --output-filename="degs_volcano" \
#                              --samples-list="data/exp_data/Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
# Rscript ./scripts/plot_dif_exp_analysis.R --exp-matrix="data/exp_data/processed/GTEx_BT_and_Oncobox_BT_norm_exp_matrix.txt" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_Norms/plots/DEGs_Volcano" \
#                              --output-filename="degs_volcano" \
#                              --samples-list="data/exp_data/GTEx_BT_and_Oncobox_BT_norm_sample_names.txt"
#
# Rscript ./scripts/plot_dif_exp_analysis.R --exp-matrix="data/exp_data/processed/GTEx_BT_w_Oncobox_BT_and_norms.txt" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_norms/plots/DEGs_Volcano" \
#                              --output-filename="degs_volcano" \
#                              --samples-list="data/exp_data/Oncobox_BT_and_GTEx_norms_sample_names_for_PCA.txt"
# ```

option_list <- list(
  make_option(c("-e", "--exp-matrix"),
              dest = "exp_matrix",
              type = "character",
              metavar = "character",
              default = "data/exp_matrix.txt",
              help = "Path to processed expression matrix [default= %default]"),
  make_option(c("s-", "--samples-list"),
              dest = "samples_list",
              type = "character",
              metavar = "character",
              default = NULL,
              help = "Path to a list of TCGA sample codes to use in analysis"),
  make_option(c("s-", "--condition-name"),
              dest = "condition_name",
              type = "character",
              metavar = "character",
              default = NULL,
              help = "a string, specifying expected comparison for `DESeq2`. For ex. `condition_Mutated_vs_Control` where `Mutated` and `Control` are values, presented in `--sample-annotations`"),
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
#             output_dir = "analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano",
#             output_filename = "degs_volcano",
#             condition_name = ""condition_SIA_vs_ANTE_normal_brain"",
#             samples_list = "data/exp_data/Oncobox_BT_and_Norms_sample_names_for_PCA.txt")
#
# opt <- list(exp_matrix = "data/exp_data/processed/GTEx_BT_w_Oncobox_BT_and_norms.txt",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv",
#             output_dir = "analyses/Oncobox_BT_vs_GTEx_norms/plots/DEGs_Volcano",
#             output_filename = "degs_volcano",
#             condition_name = NULL,
#             samples_list = "data/exp_data/Oncobox_BT_and_GTEx_norms_sample_names_for_PCA.txt")

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$output_dir)) {
  dir.create(file.path(getwd(), opt$output_dir), recursive = TRUE)
}

# Read and normaliza the expression matrix and select needed samples
exp_matrix <- load_exp_matrix(exp_matrix_path = opt$exp_matrix,
                              drop_dupl_genes = TRUE,
                              sample_list_path = opt$samples_list,
                              quantile_normalization = FALSE,
                              log10_scalling = FALSE)

# Obtain sample annotations
sample_ann <- get_exp_samples_annotations(sample_ann_table_path = opt$sample_annotations,
                                          sample_list_path = opt$samples_list)
sample_ann <- gsub(" ", "_", sample_ann)

filter_non_common_samples(exp_matrix, sample_ann)
exp_matrix <- exp_matrix[, names(sample_ann)]

# Perform DEGs analysis
res_dif_exp_genes <- dif_exp_analysis(exp_matrix, sample_ann, condition_name = opt$condition_name, shring_results = FALSE)
# res_dif_exp_genes <- res_dif_exp_genes[complete.cases(res_dif_exp_genes),]
# res_dif_exp_genes <- res_dif_exp_genes[!res_dif_exp_genes$pvalue %in% c(0),]

# Plot DEG Volcano Plot
volcano_plot <- EnhancedVolcano(res_dif_exp_genes,
                                title = "Differential Expression Analysis",
                                subtitle = NULL,
                                caption = paste0("Total = ", nrow(res_dif_exp_genes), " genes"),
                                lab = NA,
# lab = rownames(res_dif_exp_genes),
                                pCutoff = 0.05,
                                FCcutoff = 2,
                                xlab = bquote(italic(.(~Log[2] ~ "Fold Change"))),
                                ylab = bquote(italic(.(~-Log[10] ~ "P"[adj]))),
                                pointSize = 2,
                                labSize = 3,
                                colAlpha = 0.75,
                                x = "log2FoldChange",
                                y = "padj",
                                legendPosition = "bottom") +
                                      theme(text = element_text(size = 14)) +
                                      theme(plot.title = element_text(hjust = 0.5))
# xlim = c(-5, 5),
# ylim = c(0, 6),
# labCol = "dark grey",
# labFace = "bold",
# col = col_pallete,

export_analysis_plot(
  filename = opt$output_filename,
# filename = "degs_volcano_wo_inf_pval",
  plot = volcano_plot,
  path = opt$output_dir,
  scale = 1,
  width = 200,
  height = 145,
  units = "mm")

saveRDS(res_dif_exp_genes, file = file.path(opt$output_dir, paste0(opt$output_filename, "_results.RDS")))