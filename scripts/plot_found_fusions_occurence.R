#!/usr/bin/env Rscript

library("MutEffect")
source("MutEffect/R/preprocessing.R", chdir = TRUE)
source("MutEffect/R/visualization.R", chdir = TRUE)
library("optparse")

if (!require("pacman")) install.packages("pacman")
pacman::p_load("readr", "readxl", "dplyr", "grid", "ggthemes", "scales", "FactoMineR", "factoextra", "ggplot2", "ggpubr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("preprocessCore")) BiocManager::install("preprocessCore")

# ## Usage
# ```bash
# cd ./compare_BT_vs_Norms
# Rscript ./scripts/plot_found_fusions_occurence.R --fusions-annotations="data/fusions_annotation/fus_ann_ref_splice_validated.xlsx" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/found_fusions_occurence" \
#                              --output-filename="found_fusions_occurence"
# ```

option_list <- list(
  make_option(c("-d", "--fusions-annotations"),
              dest = "fusions_annotations",
              type = "character",
              metavar = "character",
              default = "data/fusions_annotation/fus_ann_ref_splice_validated.xlsx",
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
  make_option(c("-m", "--sample-annotations"),
              dest = "sample_annotations",
              type = "character",
              metavar = "character",
              default = "data/exp_data/exp_samples_annotation.csv",
              help = "Path to sample annotations database [default= %default]")
)

# opt <- list(fusions_annotations = "data/fusions_annotation/fus_ann_ref_splice_validated.xlsx",
#             output_dir = "analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/found_fusions_occurence",
#             output_filename = "found_fusions_occurence",
#             samples_list = "data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt")

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$output_dir)) {
  dir.create(file.path(getwd(), opt$output_dir), recursive = TRUE)
}

# Obtain sample annotations (Found Fusions)
fusions_annotation <- read_excel(opt$fusions_annotation)
fusions_annotation <- as.data.frame(fusions_annotation)
fusions_annotation <- fusions_annotation[fusions_annotation$VALIDATION == "TRUE",]
fusions_annotation$SAMPLE_ID <- gsub("_S[0-9]+.*", "", fusions_annotation$SAMPLE_ID)
fusions_annotation$SAMPLE_ID <- gsub("_", "-", fusions_annotation$SAMPLE_ID)

fusions_occurence <- fusions_annotation %>% group_by(FUSION) %>% summarize(Occurence = length(FUSION))
# fusions_occurence <- fusions_annotation %>% group_by(SAMPLE_ID) %>% summarize(Occurence = length(SAMPLE_ID))

# Plot intersecting drug occurence
labels <- formatC(sort(fusions_occurence$Occurence),
                  format = "f", digits = 0)

fusions_occurence_plot <- ggbarplot(fusions_occurence,
  x = "FUSION",
  y = "Occurence",
  xlab = "Fusion",
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
                     plot = fusions_occurence_plot,
                     path = opt$output_dir,
                     scale = 1,
                     width = 210 / 1.7,
                     height = 297 / 2.25,
                     units = "mm",
                     load_fonts = FALSE)

# # Plot Number of Fusions per Sample_ID

# labels <- formatC(sort(fusions_occurence$Occurence),
#                   format = "f", digits = 0)

# fusions_occurence_plot <- ggbarplot(fusions_occurence,
#   x = "SAMPLE_ID",
#   y = "Occurence",
#   xlab = "Sample_ID",
#   ylab = "Fusions",
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

# export_analysis_plot(filename = "n_fusions_per_sample_id",
#                      plot = fusions_occurence_plot,
#                      path = opt$output_dir,
#                      scale = 1,
#                      width = 210 / 1.7,
#                      height = 297 / 2.25,
#                      units = "mm",
#                      load_fonts = FALSE)
