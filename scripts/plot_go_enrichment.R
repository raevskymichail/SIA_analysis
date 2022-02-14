#!/usr/bin/env /usr/bin/Rscript

library("MutEffect")
source("MutEffect/R/preprocessing.R", chdir = TRUE)
source("MutEffect/R/visualization.R", chdir = TRUE)
library("DESeq2")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("gridExtra")
library("grid")
library("readr")
library("optparse")

# ## Usage
# ```bash
# cd ./compare_mut_gene_exp
# Rscript ./scripts/plot_go_enrichment.R --dif-exp-results-path="analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --output-dir="analyses/Oncobox_BT_vs_Norms/plots/DEGs_GO_enrichment" \
#                              --output-filename="go_enrichment" \
#                              --samples-list="data/exp_data/Oncobox_BT_and_Norms_sample_names_for_PCA.txt"
#
# Rscript ./scripts/plot_go_enrichment.R --dif-exp-results-path="analyses/GTEx_BT_vs_Oncobox_BT_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --output-dir="analyses/GTEx_BT_vs_Oncobox_BT_Norms/plots/DEGs_GO_enrichment" \
#                              --output-filename="go_enrichment" \
#                              --samples-list="data/exp_data/GTEx_BT_and_Oncobox_BT_norm_sample_names.txt"
#
# Rscript ./scripts/plot_go_enrichment.R --dif-exp-results-path="analyses/Oncobox_BT_vs_GTEx_norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_norms/plots/DEGs_GO_enrichment" \
#                              --output-filename="go_enrichment" \
#                              --samples-list="data/exp_data/Oncobox_BT_and_GTEx_norms_sample_names_for_PCA.txt"
#
# Rscript ./scripts/plot_go_enrichment.R --dif-exp-results-path="analyses/Oncobox_BT_vs_GTEx_norms/plots/DEGs_Volcano/degs_volcano_results.RDS" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_norms/plots/DEGs_GO_enrichment/FC_greater_5" \
#                              --output-filename="go_enrichment" \
#                              --abs-log-FC-treshold="5" \
#                              --samples-list="data/exp_data/Oncobox_BT_and_GTEx_norms_sample_names_for_PCA.txt"
#
## GO enrichment of dif genes intersection
#
# Rscript ./scripts/plot_go_enrichment.R --dif-exp-results-path="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection/venn_dif_genes_intersection_negative_common_DEGs.RDS" \
#                              --dif-exp-results-as-intersection \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/go_dif_genes_intersetion" \
#                              --output-filename="go_enrichment_negative_degs_intersection"
#
# Rscript ./scripts/plot_go_enrichment.R --dif-exp-results-path="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection/venn_dif_genes_intersection_positive_common_DEGs.RDS" \
#                              --dif-exp-results-as-intersection \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/go_dif_genes_intersetion" \
#                              --output-filename="go_enrichment_positive_degs_intersection"
#
# Rscript ./scripts/plot_go_enrichment.R --dif-exp-results-path="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection/venn_dif_genes_intersection_negative_common_DEGs.RDS" \
#                              --dif-exp-results-as-intersection \
#                              --abs-log-FC-treshold="5" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/go_dif_genes_intersetion/FC_greater_5" \
#                              --output-filename="go_enrichment_negative_degs_intersection"
#
# Rscript ./scripts/plot_go_enrichment.R --dif-exp-results-path="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection/venn_dif_genes_intersection_positive_common_DEGs.RDS" \
#                              --dif-exp-results-as-intersection \
#                              --abs-log-FC-treshold="5" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/go_dif_genes_intersetion/FC_greater_5" \
#                              --output-filename="go_enrichment_positive_degs_intersection"
#
# Rscript ./scripts/plot_go_enrichment.R --dif-exp-results-path="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection/venn_dif_genes_intersection_negative_common_DEGs.RDS" \
#                              --dif-exp-results-as-intersection \
#                              --abs-log-FC-treshold="2" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/go_dif_genes_intersetion/FC_greater_2" \
#                              --output-filename="go_enrichment_negative_degs_intersection"

# Rscript ./scripts/plot_go_enrichment.R --dif-exp-results-path="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection/venn_dif_genes_intersection_positive_common_DEGs.RDS" \
#                              --dif-exp-results-as-intersection \
#                              --abs-log-FC-treshold="2" \
#                              --output-dir="analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/go_dif_genes_intersetion/FC_greater_2" \
#                              --output-filename="go_enrichment_positive_degs_intersection"
# ```

option_list <- list(
  make_option(c("--dif-exp-results-path"),
              dest = "dif_exp_results_path",
              type = "character",
              help = "Path to .RDS object corresponding to the DESeq2 results data frame"),
  make_option(c("-s", "--samples-list"),
              dest = "samples_list",
              type = "character",
              metavar = "character",
              default = NULL,
              help = "Path to a list of TCGA sample codes to use in analysis"),
  make_option(c("-t", "--abs-log-FC-treshold"),
              dest = "abs_log_FC_treshold",
              type = "numeric",
              metavar = "numeric",
              default = 0,
              help = "An absolute value of log2FoldChange treshold for subsetting of differentially expressed genes"),
  make_option(c("-g", "--hgnc-set-path"),
              dest = "hgnc_set_path",
              type = "character",
              metavar = "character",
              default = "data/exp_data/hgnc_complete_set.txt",
              help = "Path to a list of HGNC Gene Symbols to use in analysis [default= %default]"),
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
  make_option(c("-r", "--dif-exp-results-as-intersection"),
              dest = "dif_exp_results_as_intersection",
              action = "store_true",
              default = FALSE,
              help = "Flag indicating whether Differential Expression results table should be processed as DEGs intersection table [default= %default]"),
make_option(c("-m", "--sample-annotations"),
              dest = "sample_annotations",
              type = "character",
              metavar = "character",
              default = "data/exp_data/exp_samples_annotation.csv",
              help = "Path to sample annotations database [default= %default]")
)

# opt <- list(dif_exp_results_path = "analyses/Oncobox_BT_vs_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv",
#             hgnc_set_path = "data/exp_data/hgnc_complete_set.txt",
#             output_dir = "analyses/Oncobox_BT_vs_Norms/plots/DEGs_GO_enrichment",
#             output_filename = "go_enrichment",
#             abs_log_FC_treshold = 0,
#             samples_list = "data/exp_data/Oncobox_BT_and_Norms_sample_names_for_PCA.txt")

# opt <- list(dif_exp_results_path = "analyses/GTEx_BT_vs_Oncobox_BT_Norms/plots/DEGs_Volcano/degs_volcano_results.RDS",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv",
#             hgnc_set_path = "data/exp_data/hgnc_complete_set.txt",
#             output_dir = "analyses/GTEx_BT_vs_Oncobox_BT_Norms/plots/DEGs_GO_enrichment",
#             output_filename = "go_enrichment",
#             abs_log_FC_treshold = 0,
#             samples_list = "data/exp_data/GTEx_BT_and_Oncobox_BT_norm_sample_names.txt")

# opt <- list(dif_exp_results_path = "analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/venn_dif_genes_intersection/venn_dif_genes_intersection_negative_common_DEGs.RDS",
#             sample_annotations = "data/exp_data/exp_samples_annotation.csv",
#             hgnc_set_path = "data/exp_data/hgnc_complete_set.txt",
#             output_dir = "analyses/Oncobox_BT_vs_GTEx_vs_ANTE_norms/plots/go_dif_genes_intersetion",
#             output_filename = "go_enrichment_negative_degs_intersection",
#             abs_log_FC_treshold = 0,
#             dif_exp_results_as_intersection = TRUE,
#             samples_list = NULL)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$output_dir)) {
  dir.create(file.path(getwd(), opt$output_dir), recursive = TRUE)
}

# Load results of Differential Genes Analysis
dif_exp_results <- readRDS(opt$dif_exp_results_path)
dif_exp_results <- as.data.frame(dif_exp_results)

if (opt$dif_exp_results_as_intersection == TRUE) {
  dif_genes_significant <- rownames(dif_exp_results[dif_exp_results$padj.x < 0.05 & abs(dif_exp_results$log2FoldChange.x) > opt$abs_log_FC_treshold,])
} else {
  dif_genes_significant <- rownames(dif_exp_results[dif_exp_results$padj < 0.05 & abs(dif_exp_results$log2FoldChange) > opt$abs_log_FC_treshold,])
}

# Calculate GO Enrichment
hgnc <- read_tsv(opt$hgnc_set_path, col_names = TRUE)
entrez <- as.character(hgnc[!is.na(match(hgnc$symbol, dif_genes_significant)), "entrez_id"]$entrez_id)
go_terms <- enrichGO(entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", readable = TRUE)

# Plot Go Enrichment results
go_enrich_dotplot <- dotplot(go_terms, showCategory = 35)

if (opt$dif_exp_results_as_intersection == TRUE) {
  foldChange <- dif_exp_results[dif_genes_significant,]$log2FoldChange.x
} else {
  foldChange <- dif_exp_results[dif_genes_significant,]$log2FoldChange
}
names(foldChange) <- entrez
go_enrich_heatplot <- heatplot(go_terms, foldChange = foldChange, showCategory = 35)

# cowplot::plot_grid(go_enrich_dotplot, go_enrich_heatplot, ncol = 1, labels = c("A", "B"))

export_analysis_plot(filename = paste0(opt$output_filename, "_dotplot"),
                     plot = go_enrich_dotplot,
                     path = opt$output_dir,
                     scale = 1,
                     width = 210 / 0.98,
                     height = 297 / 1.5,
                     units = "mm",
                     load_fonts = FALSE, cairo_pdf_device = TRUE)

export_analysis_plot(filename = paste0(opt$output_filename, "_heatplot"),
                     plot = go_enrich_heatplot,
                     path = opt$output_dir,
                     scale = 1,
                     width = 210 / 0.98,
                     height = 297 / 1.8,
                     units = "mm",
                     load_fonts = FALSE, cairo_pdf_device = TRUE)

saveRDS(go_terms, file = file.path(opt$output_dir, paste0(opt$output_filename, ".RDS")))
