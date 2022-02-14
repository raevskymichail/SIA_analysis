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
# Rscript ./scripts/plot_venn_dif_genes_intersection.R --dif-exp-results-set-a-path="analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --dif-exp-results-set-b-path="analyses/GTEx_BT_vs_Oncobox_BT_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --set-a-label="SIA" \
#                              --set-b-label="GTEx spinal cord" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/venn_dif_genes_intersection" \
#                              --output-filename="venn_dif_genes_intersection" \
#                              --samples-list="data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
# Rscript ./scripts/plot_venn_dif_genes_intersection.R --dif-exp-results-set-a-path="analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --dif-exp-results-set-b-path="analyses/GTEx_BT_vs_Oncobox_BT_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --set-a-label="SIA" \
#                              --set-b-label="GTEx spinal cord" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_and_Norms/plots/venn_dif_genes_intersection/FC_greater_5" \
#                              --output-filename="venn_dif_genes_intersection" \
#                              --samples-list="data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
# Rscript ./scripts/plot_venn_dif_genes_intersection.R --dif-exp-results-set-a-path="analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --dif-exp-results-set-b-path="analyses/Oncobox_BT_vs_GTEx_norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --set-a-label="SIA vs ANTE normal brain" \
#                              --set-b-label="SIA vs GTEx spinal cord" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection" \
#                              --output-filename="venn_dif_genes_intersection" \
#                              --samples-list="data/exp_data/Oncobox_BT_and_GTEx_norms_sample_names_for_PCA.txt"

# Rscript ./scripts/plot_venn_dif_genes_intersection.R --dif-exp-results-set-a-path="analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --dif-exp-results-set-b-path="analyses/Oncobox_BT_vs_GTEx_norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --set-a-label="SIA vs ANTE normal brain" \
#                              --set-b-label="SIA vs GTEx spinal cord" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection/FC_greater_5" \
#                              --output-filename="venn_dif_genes_intersection" \
#                              --abs-log-FC-treshold="5" \
#                              --samples-list="data/exp_data/Oncobox_BT_and_GTEx_norms_sample_names_for_PCA.txt"

# Rscript ./scripts/plot_venn_dif_genes_intersection.R --dif-exp-results-set-a-path="analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --dif-exp-results-set-b-path="analyses/Oncobox_BT_vs_GTEx_norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --set-a-label="SIA vs ANTE normal brain" \
#                              --set-b-label="SIA vs GTEx spinal cord" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection/FC_greater_2" \
#                              --output-filename="venn_dif_genes_intersection" \
#                              --abs-log-FC-treshold="2" \
#                              --samples-list="data/exp_data/Oncobox_BT_and_GTEx_norms_sample_names_for_PCA.txt"
# ```

option_list <- list(
  make_option(c("--dif-exp-results-set-a-path"),
              dest = "dif_exp_results_set_a_path",
              type = "character",
              help = "Path to .RDS object corresponding to the DESeq2 results data frame from the first (set A) experiment"),
  make_option(c("--dif-exp-results-set-b-path"),
              dest = "dif_exp_results_set_b_path",
              type = "character",
              help = "Path to .RDS object corresponding to the DESeq2 results data frame from the first (set B) experiment"),
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
  make_option(c("-t", "--abs-log-FC-treshold"),
              dest = "abs_log_FC_treshold",
              type = "numeric",
              metavar = "numeric",
              default = 0,
              help = "An absolute value of log2FoldChange treshold for subsetting of differentially expressed genes"),
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

# opt <- list(dif_exp_results_set_a_path = "analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS",
#             dif_exp_results_set_b_path = "analyses/Oncobox_BT_vs_GTEx_norms/plots/DEGs_Volcano/degs_volcano_results.RDS",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv",
#             set_a_label = "SIA vs ANTE normal brain",
#             set_b_label = "SIA vs GTEx spinal cord",
#             abs_log_FC_treshold = 2,
#             output_dir = "analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection/FC_greater_2",
#             output_filename = "venn_dif_genes_intersection",
#             samples_list = "data/exp_data/GTEx_Oncobox_BT_and_Norms_sample_names_for_PCA.txt")

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$output_dir)) {
  dir.create(file.path(getwd(), opt$output_dir), recursive = TRUE)
}

# Load results of Differential Genes Analysis
dif_exp_results_set_a <- readRDS(opt$dif_exp_results_set_a_path)
dif_exp_results_set_a <- as.data.frame(dif_exp_results_set_a)
dif_exp_results_set_a <- dif_exp_results_set_a[complete.cases(dif_exp_results_set_a),]

dif_exp_results_set_b <- readRDS(opt$dif_exp_results_set_b_path)
dif_exp_results_set_b <- as.data.frame(dif_exp_results_set_b)
dif_exp_results_set_b <- dif_exp_results_set_b[complete.cases(dif_exp_results_set_b),]

# Subsets positive and negative sets of DEGs
dif_exp_results_set_a <- dif_exp_results_set_a[dif_exp_results_set_a$padj < 0.05 & abs(dif_exp_results_set_a$log2FoldChange) > opt$abs_log_FC_treshold,]
dif_genes_set_a_significant <- rownames(dif_exp_results_set_a)
dif_genes_positive_set_a_significant <- rownames(dif_exp_results_set_a[dif_exp_results_set_a$log2FoldChange > 0,])
dif_genes_negative_set_a_significant <- rownames(dif_exp_results_set_a[dif_exp_results_set_a$log2FoldChange < 0,])

dif_exp_results_set_b <- dif_exp_results_set_b[dif_exp_results_set_b$padj < 0.05 & abs(dif_exp_results_set_b$log2FoldChange) > opt$abs_log_FC_treshold,]
dif_genes_set_b_significant <- rownames(dif_exp_results_set_b)
dif_genes_positive_set_b_significant <- rownames(dif_exp_results_set_b[dif_exp_results_set_b$log2FoldChange > 0,])
dif_genes_negative_set_b_significant <- rownames(dif_exp_results_set_b[dif_exp_results_set_b$log2FoldChange < 0,])

# Plot Venn diagram

custom_pallete <- c("#619CFF", "#FB746C") # Blue, Red
positive_FC_label <- bquote(italic(.(~Log[2] ~ "Fold Change") > .(opt$abs_log_FC_treshold)))
negative_FC_label <- bquote(italic(.(~Log[2] ~ "Fold Change") < .(- opt$abs_log_FC_treshold)))

# Venn - all DEGs

venn.diagram(
         x = list(dif_genes_set_a_significant, dif_genes_set_b_significant),
         category.names = c(opt$set_a_label, opt$set_b_label),
         filename = file.path(opt$output_dir, paste0(opt$output_filename, "_all_DEGs.tiff")),
         main = "all DEGs",
#  sub = "",

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


# Venn - positive DEGs

venn.diagram(
         x = list(dif_genes_positive_set_a_significant, dif_genes_positive_set_b_significant),
         category.names = c(opt$set_a_label, opt$set_b_label),
         filename = file.path(opt$output_dir, paste0(opt$output_filename, "_positive_DEGs.tiff")),
         main = positive_FC_label,
#  sub = "",

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

# Venn - negative DEGs

venn.diagram(
         x = list(dif_genes_negative_set_a_significant, dif_genes_negative_set_b_significant),
         category.names = c(opt$set_a_label, opt$set_b_label),
         filename = file.path(opt$output_dir, paste0(opt$output_filename, "_negative_DEGs.tiff")),
         main = negative_FC_label,
#  sub = "",

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

dif_exp_results_common_DEGs <- merge(dif_exp_results_set_a, dif_exp_results_set_b,
                                             by = 0, all = FALSE) # merge by rownames - (SYMBOL genes)
dif_exp_results_common_positive_DEGs <- dif_exp_results_common_DEGs[dif_exp_results_common_DEGs$log2FoldChange.x > 0 & dif_exp_results_common_DEGs$log2FoldChange.y > 0,]
dif_exp_results_common_negative_DEGs <- dif_exp_results_common_DEGs[dif_exp_results_common_DEGs$log2FoldChange.x < 0 & dif_exp_results_common_DEGs$log2FoldChange.y < 0,]

# dif_exp_results_common_DEGs <- dif_exp_results_common_DEGs[complete.cases(dif_exp_results_common_DEGs),]
# dif_exp_results_common_positive_DEGs <- dif_exp_results_common_positive_DEGs[complete.cases(dif_exp_results_common_positive_DEGs),]
# dif_exp_results_common_negative_DEGs <- dif_exp_results_common_negative_DEGs[complete.cases(dif_exp_results_common_negative_DEGs),]

rownames(dif_exp_results_common_DEGs) <- dif_exp_results_common_DEGs$Row.names
rownames(dif_exp_results_common_positive_DEGs) <- dif_exp_results_common_positive_DEGs$Row.names
rownames(dif_exp_results_common_negative_DEGs) <- dif_exp_results_common_negative_DEGs$Row.names

saveRDS(dif_exp_results_common_DEGs, file = file.path(opt$output_dir, paste0(opt$output_filename, "_all_common_DEGs.RDS")))
saveRDS(dif_exp_results_common_positive_DEGs, file = file.path(opt$output_dir, paste0(opt$output_filename, "_positive_common_DEGs.RDS")))
saveRDS(dif_exp_results_common_negative_DEGs, file = file.path(opt$output_dir, paste0(opt$output_filename, "_negative_common_DEGs.RDS")))

write.csv(dif_exp_results_common_DEGs, file = file.path(opt$output_dir, paste0(opt$output_filename, "_all_common_DEGs.csv")), quote = FALSE, row.names = FALSE)
write.csv(dif_exp_results_common_positive_DEGs, file = file.path(opt$output_dir, paste0(opt$output_filename, "_positive_common_DEGs.csv")), quote = FALSE, row.names = FALSE)
write.csv(dif_exp_results_common_negative_DEGs, file = file.path(opt$output_dir, paste0(opt$output_filename, "_negative_common_DEGs.csv")), quote = FALSE, row.names = FALSE)