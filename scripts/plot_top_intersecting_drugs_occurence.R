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
# Rscript ./scripts/plot_top_intersecting_drugs_occurence.R --drug-scores-db="data/drug_scores/drugs_scores_on_DEGs/Oncobox_BT_DEGs_Drug_Scores.xlsx" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/top_intersecting_drugs_occurence" \
#                              --output-filename="top_20_intersecting_drugs_occurence" \
#                              --top-n-drugs="20"
#
# Rscript ./scripts/plot_top_intersecting_drugs_occurence.R --drug-scores-db="data/drug_scores/drugs_scores_on_DEGs/Oncobox_BT_and_GTEx_norms_DEGs_Drug_Scores.xlsx" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/top_intersecting_drugs_occurence" \
#                              --output-filename="top_20_intersecting_drugs_occurence" \
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
#             output_filename = "top_20_intersecting_drugs_occurence",
#             output_dir = "analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/top_intersecting_drugs_occurence",
#             top_n_drugs = 20,
#             samples_list = "data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt")

# opt <- list(drug_scores_db = "data/drug_scores/drugs_scores_on_DEGs/Oncobox_BT_and_GTEx_norms_DEGs_Drug_Scores.xlsx",
#             output_filename = "top_20_intersecting_drugs_occurence",
#             output_dir = "analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/top_intersecting_drugs_occurence",
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

# Define list of intersecting drugs
intersecting_drugs <- c("Thalidomide",
"Pomalidomide",
"Regorafenib",
"Lenvatinib",
"Nintedanib (BIBF 1120)",
"Sorafenib",
"Dovitinib",
"Sunitinib",
"Copanlisib",
"Tivozanib",
"Pazopanib",
"Perifosine",
"Foretinib",
"Dasatinib",
"Midostaurin",
"Imatinib",
"Erdafitinib",
"Denosumab",
"Alpelisib")

drug_scores <- drug_scores[drug_scores$Drug %in% intersecting_drugs,]
rownames(drug_scores) <- drug_scores$Drug
drug_scores_case_samples <- drug_scores[,grepl("[Cc]ase_", colnames(drug_scores))]
drug_scores_case_samples$Case_geomean <- NULL
drug_scores_case_samples["BT-20_S5_R1_001"] <- NULL
colnames(drug_scores_case_samples) <- gsub("Case_", "", colnames(drug_scores_case_samples))
colnames(drug_scores_case_samples) <- gsub("_S[0-9]+.*", "", colnames(drug_scores_case_samples))
colnames(drug_scores_case_samples) <- gsub("_", "-", colnames(drug_scores_case_samples))
drug_scores_case_samples <- drug_scores_case_samples[, !duplicated(colnames(drug_scores_case_samples))]

intersecting_drugs_occurence <- data.frame(Drug = rownames(drug_scores_case_samples),
                                           Occurence = rowSums(drug_scores_case_samples > 0.1))

# Plot intersecting drug occurence
labels <- formatC(sort(intersecting_drugs_occurence$Occurence),
                  format = "f", digits = 0)

drug_occurence_plot <- ggbarplot(intersecting_drugs_occurence,
  x = "Drug",
  y = "Occurence",
  xlab = "Drug",
  ylab = "Occurence",
  # fill = "#619cff",
  # color = "#619cff", # "#04bb3b" - green, "#fb746c" - red, "#619cff" - blue
  fill = "grey",
  color = "grey", # "#04bb3b" - green, "#fb746c" - red, "#619cff" - blue
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


# # Plot Number of Intersecting Drugs per Sample_ID

# labels <- formatC(sort(intersecting_drugs_occurence$Occurence),
#                   format = "f", digits = 0)

# drug_occurence_plot <- ggbarplot(intersecting_drugs_occurence,
#   x = "Sample_ID",
#   y = "Occurence",
#   xlab = "Sample_ID",
#   ylab = "Number of Drugs",
#   # fill = "#619cff",
#   # color = "#619cff", # "#04bb3b" - green, "#fb746c" - red, "#619cff" - blue
#   fill = "grey",
#   color = "grey", # "#04bb3b" - green, "#fb746c" - red, "#619cff" - blue
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
#   ggtheme = theme_minimal(),
# ) + theme_Publication() + theme(axis.title.y = element_blank())

# export_analysis_plot(filename = "n_intersecting_drugs_per_sample_id",
#                      plot = drug_occurence_plot,
#                      path = opt$output_dir,
#                      scale = 1,
#                      width = 210 / 1.7,
#                      height = 297 / 2.25,
#                      units = "mm",
#                      load_fonts = FALSE)