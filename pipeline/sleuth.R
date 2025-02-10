#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(sleuth)
  library(dplyr)
  library(ggplot2)
  library(tidyverse)
})

# Define command line options
option_list <- list(
  make_option(c("-r", "--root"), type="character", default=NULL, help="Root path for analysis"),
  make_option(c("-m", "--metadata"), type="character", default="metadata.tsv", help="Metadata filename"),
  make_option(c("-b", "--bootstrap"), type="character", default="06_bootstrap", help="Bootstrap directory name"),
  make_option(c("-p", "--pseudocount"), type="numeric", default=0.5, help="Pseudocount value"),
  make_option(c("-a", "--alpha"), type="numeric", default=0.05, help="Alpha threshold"),
  make_option(c("-f", "--fc"), type="numeric", default=10, help="Fold change threshold"),
  make_option(c("--counts"), type="character", default="normalized_estimated_counts.tsv", help="Counts filename"),
  make_option(c("--tpm"), type="character", default="normalized_tpm.tsv", help="TPM filename"),
  make_option(c("--wald"), type="character", default="Wald_results.tsv", help="Wald results filename"),
  make_option(c("--lrt"), type="character", default="LRT_results.tsv", help="LRT results filename"),
  make_option(c("--plotting"), type="character", default="plotting.tsv", help="Plotting data filename"),
  make_option(c("--wald-plot"), type="character", default="Wald_volcano.png", help="Wald plot filename"),
  make_option(c("--lrt-plot"), type="character", default="LRT_volcano.png", help="LRT plot filename"),
  make_option(c("--parameter"), type="character", default="strain", help="Design formula parameter")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$root)) {
  stop("Root path must be specified. Use --help for more information.")
}

# Set paths
root_path <- opt$root
sample_ids <- dir(file.path(root_path, opt$bootstrap))
data_path <- file.path(root_path, opt$bootstrap, sample_ids)

# Set variables
pseudocount <- opt$pseudocount
alpha <- opt$alpha
FC_threshold <- opt$fc

# Set filenames
counts_filename <- opt$counts
tpm_filename <- opt$tpm
wald_filename <- opt$wald
lrt_filename <- opt$lrt
plotting_filename <- opt$plotting
wald_plot_filename <- opt$wald_plot
lrt_plot_filename <- opt$lrt_plot

# Validate and set plot filenames with defaults if not specified
if (is.null(opt$`wald-plot`) || opt$`wald-plot` == "") {
  wald_plot_filename <- "Wald_volcano.png"
} else {
  wald_plot_filename <- opt$`wald-plot`
}

if (is.null(opt$`lrt-plot`) || opt$`lrt-plot` == "") {
  lrt_plot_filename <- "LRT_volcano.png"
} else {
  lrt_plot_filename <- opt$`lrt-plot`
}

# Set parameter
parameter <- opt$parameter

# Define design formula
design_formula <- reformulate(parameter)

text_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.text = element_text(family = "Helvetica", size = 15, colour = "black", face = "bold"),
  legend.title = element_text(family = "Helvetica", size = 15, colour = "black", face = "bold"),
  legend.key = element_blank(),
  strip.text.x = element_text(size = 15, face = "bold"),
  plot.title = element_text(family = "Helvetica", size = 15, colour = "black", face = "bold", hjust = 0.5),
  axis.title = element_text(family = "Helvetica", size = 20, colour = "black", face = "bold"),
  axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), size = 18, face = "bold"),
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 18, face = "bold"),
  axis.text.y = element_text(family = "Helvetica", size = 15, face = "bold", color = "black"),
  axis.text.x = element_text(family = "Helvetica", size = 15, face = "bold", colour = "black", angle = 0),
  panel.background = element_rect(colour = "black"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
  axis.text = element_text(face = "bold")
)

create_volcano_plot <- function(plotting, alpha, FC_threshold, groups,
                              fc_column = "mean_log2fc", pval_column = "wt_qval",
                              sig_column = "wt_sig", plot_title = "Denim - Wald Test",
                              text_size = 4) {
  required_cols <- c(fc_column, pval_column, sig_column)
  missing_cols <- !required_cols %in% colnames(plotting)
  if (any(missing_cols)) {
    stop(sprintf("Missing required columns: %s", paste(required_cols[missing_cols], collapse = ", ")))
  }
  
  max_abs_x <- max(abs(plotting[[fc_column]]))
  min_abs_y <- min(abs(plotting[[pval_column]]))
  max_y <- max(-log10(plotting[[pval_column]]))
  
  x_margin <- max_abs_x * 0.1
  y_margin <- max_y * 0.1
  
  p <- plotting %>% 
    ggplot(aes(x = .data[[fc_column]], y = -log10(.data[[pval_column]]), col = .data[[sig_column]])) +
    geom_point(size = 3.5) +
    xlim(-max_abs_x, max_abs_x) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 1) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "dark grey", linewidth = 1.0) +
    annotate("text", x = -max_abs_x, y = -log10(alpha) + y_margin / 2.5,
             label = bquote("q-value" <= .(sprintf("%.3f", alpha))),
             hjust = 0, fontface = "bold", size = text_size, color = "dark grey") +
    geom_vline(xintercept = log2(FC_threshold), linetype = "dashed", color = "dark grey", linewidth = 1) +
    geom_vline(xintercept = -log2(FC_threshold), linetype = "dashed", color = "dark grey", linewidth = 1) +
    annotate("text", x = log2(FC_threshold) / 1.3, y = y_margin * 1.5,
             label = bquote("FC" >= .(sprintf("%.0f", FC_threshold))),
             hjust = -0.0, fontface = "bold", size = text_size, angle = 90, color = "dark grey") +
    annotate("text", x = -log2(FC_threshold) / 1.3, y = y_margin * 1.5,
             label = bquote("FC" <= .(sprintf("%.0f", -FC_threshold))),
             hjust = -0.0, fontface = "bold", size = text_size, angle = 90, color = "dark grey") +
    geom_segment(aes(x = log2(FC_threshold) + x_margin,
                    xend = log2(FC_threshold) + (3 * x_margin),
                    y = max_y, yend = max_y), 
                arrow = arrow(length = unit(0.1, "inches")), color = "dark grey") + 
    annotate("text", x = log2(FC_threshold) + x_margin, y = max_y - y_margin * 0.75,
             label = paste0("enriched in\n", groups[1]),
             hjust = 0.0, fontface = "bold", size = text_size, color = "dark grey") +
    geom_segment(aes(x = -log2(FC_threshold) - x_margin,
                    xend = -log2(FC_threshold) - (3 * x_margin),
                    y = max_y, yend = max_y), 
                arrow = arrow(length = unit(0.1, "inches")), color = "dark grey") + 
    annotate("text", x = -log2(FC_threshold) - x_margin, y = max_y - y_margin * 0.75,
             label = paste("enriched in\n", groups[2]),
             hjust = 1.0, fontface = "bold", size = text_size, color = "dark grey") +
    scale_color_manual(values = c("FALSE" = "dark grey", "TRUE" = "blue")) +
    labs(title = plot_title,
         x = expression("log"[2]*"(FC)"),
         y = expression("-log"[10]*"(q-value)"),
         color = "Significant") +
    theme_bw() +
    theme(aspect.ratio = 1, legend.position = "none") +
    text_theme
  
  return(p)
}

# Main analysis pipeline
message("Loading metadata...")
metadata <- read_tsv(file.path(root_path, opt$metadata),
                     show_col_types = FALSE) %>%
  mutate(path = file.path(root_path, opt$bootstrap, sample))

message("Preparing sleuth object...")
sleuth_object <- sleuth_prep(metadata,
                            extra_bootstrap_summary = TRUE,
                            read_bootstrap_tpm = TRUE)

message("Fitting models...")
sleuth_object <- sleuth_fit(sleuth_object, design_formula, fit_name = "full")
sleuth_object <- sleuth_fit(sleuth_object, ~1, fit_name = "reduced")
sleuth_object <- sleuth_lrt(sleuth_object, "reduced", "full")

model_info <- capture.output(models(sleuth_object))
coef_start <- grep("coefficients:", model_info)[1]
key_coef <- trimws(gsub("\\t", "", model_info[coef_start + 2]))

sleuth_object <- sleuth_wt(sleuth_object, key_coef, which_model = "full")

message("Extracting results...")
wald_results <- sleuth_results(sleuth_object, key_coef, "wt", show_all = TRUE)
lrt_results <- sleuth_results(sleuth_object, "reduced:full", "lrt", show_all = TRUE)

message("Saving count and TPM results...")
sleuth_to_matrix(sleuth_object, "obs_norm", "est_counts") %>%
  as_tibble(rownames = "target_id") %>%
  arrange(str_extract(target_id, "^[^0-9]+"),
          as.numeric(str_extract(target_id, "\\d+"))) %>%
  write_tsv(counts_filename)

sleuth_to_matrix(sleuth_object, "obs_norm", "tpm") %>%
  as_tibble(rownames = "target_id") %>%
  arrange(str_extract(target_id, "^[^0-9]+"),
          as.numeric(str_extract(target_id, "\\d+"))) %>%
  write_tsv(tpm_filename)

message("Saving statistical results...")
wald_results %>% write_tsv(wald_filename)
lrt_results %>% write_tsv(lrt_filename)

groups <- metadata %>%
  distinct(!!sym(parameter)) %>%
  pull()

message("Preparing plotting data...")
plotting <- sleuth_to_matrix(sleuth_object, "obs_norm", "est_counts") %>%
  as_tibble(rownames = "target_id") %>%
  mutate(across(-target_id, ~if_else(. == 0, . + pseudocount, .))) %>%
  pivot_longer(cols = -target_id,
               names_to = "condition",
               values_to = "value") %>%
  left_join(
    metadata %>% 
      select(sample, !!sym(parameter)) %>%
      rename(condition = sample, group = !!sym(parameter)),
    by = "condition"
  ) %>%
  group_by(target_id, group) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            .groups = 'drop') %>%
  pivot_wider(names_from = group,
              values_from = c(mean, sd),
              names_sep = "_") %>%
  mutate(
    mean_log2fc = log2(!!sym(paste0("mean_", groups[1])) / !!sym(paste0("mean_", groups[2]))),
    sd_log2fc = log2(!!sym(paste0("sd_", groups[1])) / !!sym(paste0("sd_", groups[2])))
  ) %>%
  left_join(wald_results %>%
              rename_with(~paste0("wt_", .), -target_id),
            by = "target_id") %>%
  left_join(lrt_results %>%
              rename_with(~paste0("lrt_", .), -target_id),
            by = "target_id") %>%
  arrange(str_extract(target_id, "^[^0-9]+"),
          as.numeric(str_extract(target_id, "\\d+")))

plotting %>% write_tsv(plotting_filename)

plotting <- plotting %>% 
  mutate(wt_sig = ifelse(wt_qval <= alpha & abs(mean_log2fc) >= log2(FC_threshold), TRUE, FALSE)) %>%
  mutate(lrt_sig = ifelse(lrt_qval <= alpha & abs(mean_log2fc) >= log2(FC_threshold), TRUE, FALSE)) %>%
  select(target_id, mean_log2fc, wt_qval, lrt_qval, wt_sig, lrt_sig)

options(repr.plot.width = 6, repr.plot.height = 6)

message("Creating and saving plots...")
p <- create_volcano_plot(plotting = plotting,
                        alpha = alpha,
                        FC_threshold = FC_threshold,
                        groups = groups,
                        fc_column = "mean_log2fc",
                        pval_column = "wt_qval",
                        sig_column = "wt_sig",
                        plot_title = "Denim - Wald Test")

ggsave(filename = wald_plot_filename,
       plot = p,
       dpi = 1200,
       width = 6,
       height = 6)

p <- create_volcano_plot(plotting = plotting,
                        alpha = alpha,
                        FC_threshold = FC_threshold,
                        groups = groups,
                        fc_column = "mean_log2fc",
                        pval_column = "lrt_qval",
                        sig_column = "lrt_sig",
                        plot_title = "Denim - Likelihood Ratio Test")

ggsave(filename = lrt_plot_filename,
       plot = p,
       dpi = 1200,
       width = 6,
       height = 6)
