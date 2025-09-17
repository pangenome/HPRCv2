options(scipen = 10000)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Read the data
num_haplo <- 466
num_sample <- 234
data <- read_tsv("/home/guarracino/Desktop/hprc25272.CHM13.w100k-xm3-l3000.tsv.gz")
# Extract chromosome information from the chrom column
data <- data %>%
  mutate(
    chromosome = gsub("CHM13#0#(chr[^\\s]+).*", "\\1", chrom),
    # Extract just the number/letter from the chromosome
    chrom_num = gsub("chr([0-9XYM]+)", "\\1", chromosome),
    # Create numeric position for sorting (X, Y, M at the end)
    chrom_order = case_when(
      chrom_num == "X" ~ 23,
      chrom_num == "Y" ~ 24,
      chrom_num == "M" ~ 25,
      TRUE ~ as.numeric(chrom_num)
    )
  )

# Reorder the factor levels to ensure chromosomes are in the correct order
data$chromosome <- factor(data$chromosome, 
                          levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))

# Create a position variable for x-axis that represents the bin midpoint
data <- data %>%
  mutate(position = (start + end) / 2 / 1e6)  # Convert to Mbp

# First reshape the data to long format for the first combined plot (alignments)
data_alignments <- data %>%
  pivot_longer(
    cols = c(num_alignments, num_alignments_merged),
    names_to = "metric",
    values_to = "value"
  )

# Create a more informative legend for alignments
data_alignments$metric <- factor(data_alignments$metric,
                                 levels = c("num_alignments", "num_alignments_merged"),
                                 labels = c("Alignments", "Alignments Merged"))

#data_alignments <- data_alignments %>% filter(metric == "Alignments Merged")

# Create combined alignments plot with more x ticks and bigger text
p_combined_alignments <- ggplot(data_alignments, aes(x = position, y = value, color = metric, group = metric))

p_combined_alignments <- p_combined_alignments +
  geom_line() +
  facet_wrap(~ chromosome, scales = "free", ncol = 5) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +  # More x-axis ticks
  scale_color_manual(values = c("blue", "red")) +
  labs(
    title = "Alignments across CHM13 (100kbp windows)",
    x = "Position (Mbp)",
    y = "Count (log10)",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 13),  # Larger facet text
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom",
    legend.title = element_text(size = 15, face = "bold"),  # Larger legend title
    legend.text = element_text(size = 13),  # Larger legend values
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14)
  )

# Reshape data for the second combined plot (haplotypes and samples)
data_haplo_samples <- data %>%
  pivot_longer(
    cols = c(num_haplotypes, num_samples),
    names_to = "metric",
    values_to = "value"
  )

# Create a more informative legend for haplotypes and samples
data_haplo_samples$metric <- factor(data_haplo_samples$metric,
                                    levels = c("num_haplotypes", "num_samples"),
                                    labels = c("Haplotypes", "Samples"))

# Create combined haplotypes and samples plot with more x ticks and bigger text
p_combined_haplo_samples <- ggplot(data_haplo_samples, aes(x = position, y = value, color = metric, group = metric))

p_combined_haplo_samples <- p_combined_haplo_samples +
  geom_line() +
  facet_wrap(~ chromosome, scales = "free", ncol = 5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(0, NA)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6),
                     expand = expansion(mult = c(0.09, 0.15))) +  # Add padding to x-axis
  scale_color_manual(values = c("darkgreen", "purple")) +
  labs(
    title = "Samples/Haplotypes across CHM13 (100kbp windows)",
    x = "Position (Mbp)",
    y = "Count",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 13),  # Larger facet text
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom",
    legend.title = element_text(size = 15, face = "bold"),  # Larger legend title
    legend.text = element_text(size = 13),  # Larger legend values
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14)
  ) +
  # Add horizontal dashed reference lines for all facets
  geom_hline(yintercept = num_haplo, linetype = "dashed", color = "darkgreen", alpha = 0.7) +
  geom_hline(yintercept = num_sample, linetype = "dashed", color = "purple", alpha = 0.7) +
  # Add text annotations for all facets
  annotate("text", x = 0, y = num_haplo, label = paste(num_haplo), 
           color = "darkgreen", hjust = +0.9, vjust = +1.5, size = 3) +
  annotate("text", x = 0, y = num_sample, label = paste(num_sample), 
           color = "purple", hjust = +0.9, vjust = +1.5, size = 3)

# Get the levels from the original data to ensure consistent factor levels
chrom_levels <- levels(data_haplo_samples$chromosome)

# Create a dataframe for the chrX horizontal line
chrX_line_data <- data.frame(
  yintercept = 348,
  chromosome = "chrX"
)
# Apply the same factor levels as the original data
chrX_line_data$chromosome <- factor(chrX_line_data$chromosome, levels = chrom_levels)

# Create a dataframe for the chrX annotation text
chrX_text_data <- data.frame(
  x = 0,
  y = 348,
  label = "348",
  chromosome = "chrX"
)
# Apply the same factor levels as the original data
chrX_text_data$chromosome <- factor(chrX_text_data$chromosome, levels = chrom_levels)

# Add the chrX-specific horizontal line and text annotation
p_combined_haplo_samples <- p_combined_haplo_samples +
  geom_hline(data = chrX_line_data, 
             aes(yintercept = yintercept), 
             linetype = "dashed", color = "black", alpha = 0.7) +
  geom_text(data = chrX_text_data, 
            aes(x = x, y = y, label = label), 
            inherit.aes = FALSE,  # Don't inherit aesthetics from the main plot
            color = "black", hjust = +0.9, vjust = +1.5, size = 3)

p_num_chromosomes <- ggplot(data, aes(x = position, y = num_chromosomes))

p_num_chromosomes <- p_num_chromosomes +
  geom_line(color = "darkorange", linewidth = 0.8) +
  facet_wrap(~ chromosome, scales = "free", ncol = 5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(0, NA)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6),
                     expand = expansion(mult = c(0.09, 0.15))) +  # Add padding to x-axis
  labs(
    title = "Number of Unique Chromosomes per Region across CHM13 (100kbp windows)",
    x = "Position (Mbp)",
    y = "Number of Chromosomes"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 13),  # Larger facet text
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14)
  ) #+
  # Add horizontal reference lines for expected number of autosomes (22) and all chromosomes (24 without chrM, or 25 with chrM)
  #geom_hline(yintercept = 22, linetype = "dashed", color = "gray50", alpha = 0.7) +
  #geom_hline(yintercept = 24, linetype = "dashed", color = "gray30", alpha = 0.7) #+
  # Add text annotations for reference lines
  #annotate("text", x = 0, y = 22, label = "22 (autosomes)", 
  #         color = "gray50", hjust = +0.9, vjust = +1.5, size = 3) +
  #annotate("text", x = 0, y = 24, label = "24 (all chr)", 
  #         color = "gray30", hjust = +0.9, vjust = +1.5, size = 3)

# Print the plots
print(p_combined_alignments)
print(p_combined_haplo_samples)
print(p_num_chromosomes)

width=16
height=9
dpi=300

ggsave(
  filename = "p_combined_alignments.pdf",
  plot = p_combined_alignments,
  width = width,
  height = height,
  dpi = dpi,
  units = "in",
  bg = "white"
)
ggsave(
  filename = "p_combined_haplo_samples.pdf",
  plot = p_combined_haplo_samples,
  width = width,
  height = height,
  dpi = dpi,
  units = "in",
  bg = "white"
)
ggsave(
  filename = "p_num_chromosomes.pdf",
  plot = p_num_chromosomes,
  width = width,
  height = height,
  dpi = dpi,
  units = "in",
  bg = "white"
)
