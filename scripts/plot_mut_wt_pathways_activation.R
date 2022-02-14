#!/usr/bin/env /usr/bin/Rscript

library("MutEffect")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("optparse")

# ## Usage
# ```bash
# cd ./compare_mut_gene_exp
# /usr/bin/Rscript ./scripts/main/plot_mut_wt_pathways_activation.R --gene="EGFR" \
#                              --pathways-activations-db="data/pathways_activation_scores/EGF_PWs.xlsx" \
#                              --output-dir="analyses/EGFR/any/lung_cancer" \
#                              --output-filename="EGFR_mut_wt_pathways_activation" \
#                              --samples-list="data/TCGA_lung_cancer_sample_names.txt"
# ```

option_list <- list(
  make_option(c("-m", "--mutations-sample-annotations"),
              dest = "mutations_sample_annotations",
              type = "character",
              metavar = "character",
              default = "data/TCGA_VCFs_found_AA_changes.xlsx",
              help = "Path to sample annotations (mutations - AA changes) database [default= %default]"),
  make_option(c("-d", "--pathways-activations-db"),
              dest = "pathways_activations_db",
              type = "character",
              metavar = "character",
              help = "Path to sample annotations (drug scores) database [default= %default]"),
  make_option(c("-g", "--gene"),
              type = "character",
              metavar = "character",
              help = "HGNC name of a mutated gene (ex. `BRAF`) or a distinct AA change (ex. `BRAF p.V600E`)"),
  make_option(c("--aa-change"),
              dest = "aa_change",
              type = "character",
              metavar = "character",
              help = "HGNC name of a distinct AA change (ex. `BRAF p.V600E`)"),
  make_option(c("-r", "--return-ranks"),
              dest = "return_ranks",
              action = "store_true",
              default = FALSE,
              help = "Flag indicating whether Drug Score / Pathway Activation values should be presented as ranks [default= %default]"),
  make_option(c("-s", "--samples-list"),
              dest = "samples_list",
              type = "character",
              metavar = "character",
              help = "Path to a list of TCGA sample codes to use in analysis"),
  make_option(c("--set-control-wo-any-genes-mut"),
              dest = "set_control_wo_any_genes_mut",
              action = "store_true",
              default = FALSE,
              help = "Flag indicating whether to set controls as samples without any mutations in `check_gene_mutated` AND without any mutations in other genes (BRAF, *RAS, etc.)"),
  make_option(c("--set-control-as-norms"),
              dest = "set_control_as_norms",
              action = "store_true",
              default = FALSE,
              help = "Flag indicating whether to sets controls as adjacent normal samples"),
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
  make_option(c("--label-height"),
              dest = "label_height",
              default = 75,
              help = "Height value for a stats comparison label [default= %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$output_dir)) {
  dir.create(file.path(getwd(), opt$output_dir), recursive = TRUE)
}

# Obtain sample annotations (mutations - AA changes)
sample_mut_ann <- get_mut_sample_annotations(sample_ann_table_path = opt$mutations_sample_annotations,
                                             check_gene_mutated = opt$gene,
                                             check_aa_change = opt$aa_change,
                                             sample_list_path = opt$samples_list,
                                             set_control_wo_any_genes_mut = opt$set_control_wo_any_genes_mut,
                                             set_control_as_norms = opt$set_control_as_norms)

# Obtain sample annotations (Pathways Activations)
sample_pw_ann_list <- get_pathways_activation_sample_annotations(pathways_activation_db_path = opt$pathways_activations_db,
                                                                 return_ranks = opt$return_ranks)
sample_pw_ann <- unlist(sample_pw_ann_list)
sample_pw_ann <- data.frame(sample_names = gsub(".*\\.", "", names(sample_pw_ann)),
                            pathways_names = gsub("\\..*", "", names(sample_pw_ann)),
                            pw_activation_score = as.numeric(sample_pw_ann))

# Drop non-common samples
common_samples <- intersect(names(sample_mut_ann), sample_pw_ann$sample_names)
sample_mut_ann <- sample_mut_ann[common_samples]
sample_pw_ann <- sample_pw_ann[sample_pw_ann$sample_names %in% common_samples,]

# Plot the data

gene_data <- data.frame(sample_names = sample_pw_ann$sample_names,
                        is_mutated = sample_mut_ann[sample_pw_ann$sample_names],
                        pathways_names = gsub(" - ", "\n", sample_pw_ann$pathways_names),
                        pw_activation_score = sample_pw_ann$pw_activation_score)

# Violine plot

if (is.null(opt$aa_change)) {
  title_label <- paste("Is", opt$gene, "mutated")
} else {
  title_label <- paste("Is", opt$gene, opt$aa_change, "mutated")
}

sample_size <- gene_data %>% dplyr::count(is_mutated)

violine_plot <- gene_data %>%
                  left_join(sample_size) %>%
                  mutate(x_w_sample_size = paste0(is_mutated, "\n", "(n = ", n, ")")) %>%
                  ggplot(aes(x = pathways_names, y = pw_activation_score, fill = is_mutated)) +
                  geom_violin(position = position_dodge(width = 0.9, preserve = "total")) +
                  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9), alpha = 0.4) +
                  geom_boxplot(width = 0.15, alpha = 0.75, color = "grey",
                               position = position_dodge(width = 0.9, preserve = "total")) +
                  stat_compare_means(method = "wilcox.test", label.x = 1.35, label.y = opt$label_height) +
                  labs(title = "Pathways Activation",
                      x = "Pathway Name",
                      y = "Activation Score") +
                      theme(text = element_text(size = 14)) +
                      theme_Publication() +
                      scale_fill_Publication(name = title_label) +
                      theme(legend.position = "top", legend.title = element_text(size = 11),
                            legend.key.size = unit(0.4, "cm"))

export_analysis_plot(
  filename = opt$output_filename,
  plot = violine_plot,
  path = opt$output_dir,
  scale = 1,
  width = 6,
  height = 6,
  units = "in")