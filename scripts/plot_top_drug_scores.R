#!/usr/bin/env Rscript

library("MutEffect")
source("MutEffect/R/preprocessing.R", chdir = TRUE)
source("MutEffect/R/visualization.R", chdir = TRUE)
library("optparse")

if (!require("pacman")) install.packages("pacman")
pacman::p_load("readr", "readxl", "reshape2", "grid", "ggthemes", "scales", "FactoMineR", "factoextra", "ggplot2", "ggpubr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("preprocessCore")) BiocManager::install("preprocessCore")

# ## Usage
# ```bash
# cd ./compare_BT_vs_Norms
# Rscript ./scripts/plot_top_drug_scores.R --drug-scores-db="data/drug_scores/drugs_scores_on_DEGs/Oncobox_BT_DEGs_Drug_Scores.xlsx" \
#                              --output-dir="analyses/Oncobox_BT_vs_Norms/plots/top_drug_scores" \
#                              --output-filename="top_20_drugs_scores" \
#                              --top-n-drugs="20"
#
# Rscript ./scripts/plot_top_drug_scores.R --drug-scores-db="data/drug_scores/drugs_scores_on_DEGs/Oncobox_BT_and_GTEx_norms_DEGs_Drug_Scores.xlsx" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_norms/plots/top_drug_scores" \
#                              --output-filename="top_20_drugs_scores" \
#                              --top-n-drugs="20"
# ```

option_list <- list(
  make_option(c("-d", "--drug-scores-db"),
              dest = "drug_scores_db",
              type = "character",
              metavar = "character",
              help = "Path to sample annotations (drug scores) database (.xlsx) [default= %default]"),
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
  make_option(c("-r", "--top-n-drugs"),
              dest = "top_n_drugs",
              default = 20,
              help = "Parameter indicating quantaty of drugs to show up in ranking chart [default= %default]"),
  make_option(c("-m", "--sample-annotations"),
              dest = "sample_annotations",
              type = "character",
              metavar = "character",
              default = "data/exp_data/exp_samples_annotation.csv",
              help = "Path to sample annotations database [default= %default]")
)

# opt <- list(drug_scores_db = "data/drug_scores/drugs_scores_on_DEGs/Oncobox_BT_DEGs_Drug_Scores.xlsx",
#             output_dir = "analyses/Oncobox_BT_vs_Norms/plots/top_drug_scores",
#             output_filename = "top_20_drugs_scores",
#             top_n_drugs = 20,
#             samples_list = "data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt")
#
# opt <- list(drug_scores_db = "data/drug_scores/drugs_scores_on_DEGs/Oncobox_BT_and_GTEx_norms_DEGs_Drug_Scores.xlsx",
#             output_dir = "analyses/Oncobox_BT_vs_GTEx_norms/plots/top_drug_scores",
#             output_filename = "top_20_drugs_scores",
#             top_n_drugs = 20,
#             samples_list = "data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt")

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$output_dir)) {
  dir.create(file.path(getwd(), opt$output_dir), recursive = TRUE)
}


# Obtain sample annotations (Drug Scores)
drug_scores <- read_excel(opt$drug_scores_db, sheet = "DS_BALANCED")
drug_scores <- as.data.frame(drug_scores)
drug_scores <- drug_scores[order(drug_scores$Case_geomean, decreasing = TRUE),]
drug_scores$Case_geomean <- as.numeric(round(drug_scores$Case_geomean, 2))
drug_scores <- drug_scores[1:opt$top_n_drugs,]

# # Obtain sample annotations
# sample_ann <- get_exp_samples_annotations(sample_ann_table_path = opt$sample_annotations,
#                                           sample_list_path = opt$samples_list)

# # Drop non-common samples
# common_samples <- intersect(names(sample_ann), colnames(pal_matrix))
# sample_ann <- sample_ann[common_samples]
# pal_matrix <- pal_matrix[, common_samples]

# Plot top drugs scores
labels <- formatC(rev(round(drug_scores$Case_geomean, 1)),
                  format = "f", digits = 1)

top_drugs_plot <- ggbarplot(drug_scores,
  x = "Drug",
  y = "Case_geomean",
  xlab = "Drug",
  ylab = "Drug Score",
  # fill = "#619cff",
  # color = "#619cff", # "#04bb3b" - green, "#fb746c" - red, "#619cff" - blue
  fill = "#fb746c",
  color = "#fb746c", # "#04bb3b" - green, "#fb746c" - red, "#619cff" - blue
  sort.val = "asc",
  sort.by.groups = FALSE,
  rotate = TRUE,
  label = labels,
  lab.col = "white",
  lab.size = 4,
  lab.pos = "out",
  lab.vjust = 0.5,
  lab.hjust = 1.1,
  lab.nb.digits = NULL,
  ggtheme = theme_minimal(),
) + theme_Publication() + theme(axis.title.y = element_blank())

export_analysis_plot(filename = opt$output_filename,
                     plot = drug_occurence_plot,
                     path = opt$output_dir,
                     scale = 1,
                     width = 210 / 1.7,
                     height = 297 / 2.25,
                     units = "mm",
                     load_fonts = FALSE)

##########################
# For plot with whiskers #
##########################

# rownames(drug_scores) <- drug_scores$Drug
# drug_scores_case_samples <- drug_scores[,grepl("[Cc]ase_", colnames(drug_scores))]
# drug_scores_case_geomean <- drug_scores_case_samples$Case_geomean
# drug_scores_case_samples$Case_geomean <- NULL
# drug_scores_case_samples["BT-20_S5_R1_001"] <- NULL
# colnames(drug_scores_case_samples) <- gsub("Case_", "", colnames(drug_scores_case_samples))
# colnames(drug_scores_case_samples) <- gsub("_S[0-9]+.*", "", colnames(drug_scores_case_samples))
# colnames(drug_scores_case_samples) <- gsub("_", "-", colnames(drug_scores_case_samples))
# drug_scores_case_samples <- drug_scores_case_samples[, !duplicated(colnames(drug_scores_case_samples))]

# # Convert to long format
# drug_scores_case_samples$Drug <- rownames(drug_scores_case_samples)
# drug_scores <- melt(drug_scores_case_samples, id.vars = "Drug")
# colnames(drug_scores) <- c("Drug", "SampleID", "Drug Score")

# # Plot top drugs scores
# labels <- formatC(rev(round(drug_scores_case_geomean, 1)),
#                   format = "f", digits = 1) # For a plot with whiskers

# top_drugs_plot <- ggbarplot(drug_scores,
#   x = "Drug",
#   y = "Drug Score",
#   # fill = "#619cff",
#   # color = "#619cff", # "#04bb3b" - green, "#fb746c" - red, "#619cff" - blue
#   # fill = "#fb746c",
#   # color = "#fb746c", # "#04bb3b" - green, "#fb746c" - red, "#619cff" - blue
#   sort.val = "asc",
#   sort.by.groups = FALSE,
#   rotate = TRUE,
#   label = labels,
#   lab.col = "white",
#   lab.size = 4,
#   lab.pos = "out",
#   lab.vjust = 0.5,
#   lab.hjust = 1.1,
#   lab.nb.digits = NULL,
#   add = c("mean_sd"),
#   ggtheme = theme_minimal(),
# ) + theme_Publication() + theme(axis.title.y = element_blank())