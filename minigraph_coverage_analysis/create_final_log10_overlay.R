#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)

# Read spread bin coverage with proper node distribution
cat("Reading spread bin coverage data...\n")
spread_data <- read.csv("bin_coverage_spread.tsv", sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE)

# Clean chromosome names
spread_data$chromosome <- gsub("CHM13#", "", spread_data$chromosome)

# Get all chromosomes
all_chroms <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# Filter to actual chromosomes
data <- spread_data %>%
  filter(chromosome %in% all_chroms) %>%
  mutate(
    # Calculate ACTUAL COVERAGE for each component (not fractions!)
    total_coverage = total_bp / 100000,
    ref_coverage = ref_bp / 100000,
    nonref_coverage = nonref_bp / 100000,
    nested_coverage = nested_bp / 100000,

    # Log10 transform all coverage values (add small value to handle zeros)
    log10_total = log10(total_coverage + 0.01),
    log10_ref = log10(ref_coverage + 0.01),
    log10_nonref = log10(nonref_coverage + 0.01),
    log10_nested = log10(nested_coverage + 0.01),

    # Position in Mb
    bin_mb = bin_start / 1e6
  )

# Order chromosomes
data$chromosome <- factor(data$chromosome, levels = all_chroms)

# Calculate summary statistics for subtitle
total_gb <- sum(data$total_bp) / 1e9
nonref_gb <- sum(data$nonref_bp) / 1e9
nonref_pct <- 100 * sum(data$nonref_bp) / sum(data$total_bp)
max_coverage <- max(data$total_coverage)
n_bins_high <- sum(data$total_coverage > 10)

# Print statistics
cat("\nSummary for visualization:\n")
cat(sprintf("  Total sequence: %.2f Gb\n", total_gb))
cat(sprintf("  Non-reference: %.2f Gb (%.1f%%)\n", nonref_gb, nonref_pct))
cat(sprintf("  Maximum coverage: %.1fx\n", max_coverage))
cat(sprintf("  Bins with >10x coverage: %d\n", n_bins_high))

# Create consistent-scale visualization
cat("\nCreating final log10 overlay visualization...\n")

theme_final <- theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 5),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9),
    plot.caption = element_text(size = 7, hjust = 0, color = "gray40"),
    strip.text.y = element_text(size = 6, angle = 0),
    strip.background = element_rect(fill = "gray95", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2),
    panel.spacing = unit(0.05, "lines"),
    legend.position = "top",
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, -5, 0)
  )

# OVERLAID VERSION - All on same log10 coverage scale
overlay_data <- data %>%
  select(chromosome, bin_mb,
         `Total sequence` = log10_total,
         `Non-reference sequence` = log10_nonref,
         `Nested sequence` = log10_nested) %>%
  pivot_longer(cols = -c(chromosome, bin_mb),
               names_to = "metric",
               values_to = "log10_coverage")

overlay_data$metric <- factor(overlay_data$metric,
                              levels = c("Total sequence",
                                       "Non-reference sequence",
                                       "Nested sequence"))

p_overlay <- ggplot(overlay_data, aes(x = bin_mb, y = log10_coverage,
                                      color = metric)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  facet_grid(chromosome ~ ., scales = "fixed", switch = "y") +

  scale_color_manual(values = c("Total sequence" = "darkblue",
                               "Non-reference sequence" = "darkorange",
                               "Nested sequence" = "darkviolet")) +

  scale_x_continuous(breaks = seq(0, 250, 50), expand = c(0.01, 0)) +

  scale_y_continuous(limits = c(-2, 2.5),
                    breaks = c(-2, -1, 0, 1, 2),
                    labels = c("0.01x", "0.1x", "1x", "10x", "100x"),
                    name = "Coverage (log10 scale)") +

  labs(title = "HPRCv2: Minigraph Sequence vs CHM13",
       subtitle = sprintf("Coverage of minigraph nodes in 100kb bins across CHM13 reference | Total: %.2f Gb | Non-reference: %.2f Gb (%.1f%%) | Max: %.0fx",
                         total_gb, nonref_gb, nonref_pct, max_coverage),
       caption = paste("Each line shows sequence depth in 100kb bins. Most bins have ~1x coverage (reference only).",
                      sprintf("High coverage regions indicate structural variation. %d bins exceed 10x coverage.", n_bins_high),
                      "\nBlue = all minigraph sequences | Orange = non-CHM13 sequences | Purple = nested sequences (>1 hop from reference)",
                      sep = "\n"),
       x = "Position (Mb)") +
  theme_final

# Save with higher quality
pdf("HPRCv2_minigraph_vs_CHM13_log10.pdf", width = 14, height = 16)
print(p_overlay)
dev.off()

# Also create a simpler version without nested for clarity
cat("\nCreating simplified version (without nested)...\n")

overlay_simple <- data %>%
  select(chromosome, bin_mb,
         `Total sequence` = log10_total,
         `Non-reference sequence` = log10_nonref) %>%
  pivot_longer(cols = -c(chromosome, bin_mb),
               names_to = "metric",
               values_to = "log10_coverage")

p_simple <- ggplot(overlay_simple, aes(x = bin_mb, y = log10_coverage,
                                       color = metric)) +
  geom_line(linewidth = 0.5, alpha = 0.9) +
  facet_grid(chromosome ~ ., scales = "fixed", switch = "y") +

  scale_color_manual(values = c("Total sequence" = "darkblue",
                               "Non-reference sequence" = "darkorange")) +

  scale_x_continuous(breaks = seq(0, 250, 50), expand = c(0.01, 0)) +

  scale_y_continuous(limits = c(-2, 2.5),
                    breaks = c(-2, -1, 0, 1, 2),
                    labels = c("0.01x", "0.1x", "1x", "10x", "100x"),
                    name = "Coverage (log10 scale)") +

  labs(title = "HPRCv2: Minigraph Sequence vs CHM13",
       subtitle = sprintf("Coverage of minigraph nodes in 100kb bins across CHM13 reference | Total: %.2f Gb | Non-reference: %.2f Gb (%.1f%%)",
                         total_gb, nonref_gb, nonref_pct),
       caption = "Blue = total minigraph sequence depth | Orange = non-CHM13 sequence depth\nMost bins show ~1x coverage (CHM13 reference only). Peaks indicate structural variation hotspots.",
       x = "Position (Mb)") +
  theme_final

pdf("HPRCv2_minigraph_vs_CHM13_simple.pdf", width = 14, height = 16)
print(p_simple)
dev.off()

cat("\nDone! Created final visualizations:\n")
cat("  - HPRCv2_minigraph_vs_CHM13_log10.pdf (with nested sequences)\n")
cat("  - HPRCv2_minigraph_vs_CHM13_simple.pdf (simplified version)\n")
cat("\nTitle and descriptions updated as requested.\n")