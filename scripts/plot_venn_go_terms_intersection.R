#!/usr/bin/env /usr/bin/Rscript

library("MutEffect")
source("MutEffect/R/preprocessing.R", chdir = TRUE)
source("MutEffect/R/visualization.R", chdir = TRUE)
library("DESeq2")
library("ggplot2")
library("ggpubr")
library("VennDiagram")
library("dplyr")
library("gridExtra")
library("grid")
library("readr")
library("optparse")

# ## Usage
# ```bash
# cd ./compare_mut_gene_exp
# Rscript ./scripts/plot_venn_dif_genes_intersection.R --go-enrich-results-set-a-path="analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --go-enrich-results-set-b-path="analyses/GTEx_BT_vs_Oncobox_BT_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --set-a-label="SIA" \
#                              --set-b-label="GTEx spinal cord" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/venn_dif_genes_intersection" \
#                              --output-filename="venn_dif_genes_intersection" \
#                              --samples-list="data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
# Rscript ./scripts/plot_venn_dif_genes_intersection.R --go-enrich-results-set-a-path="analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --go-enrich-results-set-b-path="analyses/GTEx_BT_vs_Oncobox_BT_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --set-a-label="SIA" \
#                              --set-b-label="GTEx spinal cord" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/venn_dif_genes_intersection/FC_greater_5" \
#                              --output-filename="venn_dif_genes_intersection" \
#                              --abs-log-FC-treshold="5" \
#                              --samples-list="data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
# ```

option_list <- list(
  make_option(c("--go-enrich-results-set-a-path"),
              dest = "go_enrich_results_set_a_path",
              type = "character",
              help = "Path to .RDS object corresponding to the GO enrichment results data frame from the first (set A) experiment"),
  make_option(c("--go-enrich-results-set-b-path"),
              dest = "go_enrich_results_set_b_path",
              type = "character",
              help = "Path to .RDS object corresponding to the GO enrichment results data frame from the first (set B) experiment"),
  make_option(c("-s", "--samples-list"),
              dest = "samples_list",
              type = "character",
              metavar = "character",
              default = NULL,
              help = "Path to a list of TCGA sample codes to use in analysis"),
  make_option(c("-a", "--set-a-label"),
              dest = "set_a_label",
              type = "character",
              metavar = "character",
              default = NULL,
              help = "String specifying label which will be used for DEGs for the first (set A) experiment"),
  make_option(c("-b", "--set-b-label"),
              dest = "set_b_label",
              type = "character",
              metavar = "character",
              default = NULL,
              help = "String specifying label which will be used for DEGs for the second (set b) experiment"),
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

# opt <- list(go_enrich_results_set_a_path = "analyses/Oncobox_BT_vs_Norms/plots/DEGs_GO_enrichment/go_enrichment.RDS",
#             go_enrich_results_set_b_path = "analyses/GTEx_BT_vs_Oncobox_BT_Norms/plots/DEGs_GO_enrichment/go_enrichment.RDS",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv",
#             set_a_label = "SIA",
#             set_b_label = "GTEx spinal cord",
#             output_dir = "analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/venn_go_enrich_intersection",
#             output_filename = "venn_go_enrich_intersection",
#             samples_list = "data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt")

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$output_dir)) {
  dir.create(file.path(getwd(), opt$output_dir), recursive = TRUE)
}

# Load results of GO enrichment analysis
go_enrich_results_set_a <- readRDS(opt$go_enrich_results_set_a_path)
go_enrich_results_set_a <- go_enrich_results_set_a$ID[go_enrich_results_set_a$p.adjust < 0.05]

go_enrich_results_set_b <- readRDS(opt$go_enrich_results_set_b_path)
go_enrich_results_set_b <- go_enrich_results_set_b$ID[go_enrich_results_set_b$p.adjust < 0.05]

# Plot Venn diagram

custom_pallete <- c("#619CFF", "#FB746C") # Blue, Red

venn.diagram(
         x = list(go_enrich_results_set_a, go_enrich_results_set_b),
         category.names = c(opt$set_a_label, opt$set_b_label),
         filename = file.path(opt$output_dir, paste0(opt$output_filename, ".tiff")),
## Output features
 imagetype = "tiff",
 height = 297 / 3,
 width = 210 / 2,
 units = "mm",
 resolution = 300,
 compression = "lzw",

## Circles
         lwd = 2,
         lty = "blank",
         fill = custom_pallete,

         cex = 1,
         fontfamily = "sans",

## Set names
         cat.cex = 1,
         main.cex = 1,
         sub.cex = 0.5,
# Test
# hyper.test = TRUE,
# total.population = 30000
)


## Save obtained intersections

go_enrich_results_common_GO_terms <- intersect(go_enrich_results_set_a, go_enrich_results_set_b)

saveRDS(go_enrich_results_common_GO_terms, file = file.path(opt$output_dir, paste0(opt$output_filename, "common_GO_terms.RDS")))