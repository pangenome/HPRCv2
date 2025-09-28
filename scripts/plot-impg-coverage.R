options(scipen = 10000)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Optional: Set BED file path here (set to NULL if no BED file)
bed_file_path <- '/home/guarracino/Dropbox/git/HPRCv2/data/chm13-annotations.bed'  # Change this to your BED file path

# Function to read and process BED file
read_bed_regions <- function(bed_path) {
  if (!is.null(bed_path) && file.exists(bed_path)) {
    bed_data <- read_tsv(bed_path, 
                         col_names = c("chromosome", "start", "end", "name", "score", "strand"),
                         col_types = cols(.default = "c", start = "d", end = "d"))
    
    # Convert positions to Mbp to match the plot scale
    bed_data <- bed_data %>%
      mutate(
        start_mbp = start / 1e6,
        end_mbp = end / 1e6,
        # Ensure chromosome names match the format in the main data
        chromosome = ifelse(grepl("^chr", chromosome), chromosome, paste0("chr", chromosome)),
        # Create a factor for the name column for consistent coloring
        name = as.factor(name)
      )
    
    return(bed_data)
  } else {
    return(NULL)
  }
}

# Read the BED file
bed_regions <- read_bed_regions(bed_file_path) %>%
  filter(name != 'PHR-sex') %>%
  mutate(name = if_else(name == 'PHR-acro', 'PHR', name)) %>%
  mutate(name = if_else(name %in% c('PAR1', 'PAR2'), 'PAR', name)) %>%
  mutate(name = if_else(name %in% c('XTR1', 'XTR2'), 'XTR', name)) %>%
  mutate(name = if_else(name == 'Centromere', 'CEN', name))

# Create color palette for BED regions if they exist
if (!is.null(bed_regions)) {
  unique_labels <- unique(bed_regions$name)
  n_labels <- length(unique_labels)
  
  # Create a color palette - you can customize these colors
  # Using a colorblind-friendly palette
  if (n_labels <= 8) {
    bed_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                    "#0072B2", "#D55E00", "#CC79A7")[1:n_labels]
  } else {
    # For more than 8 labels, use a larger palette
    bed_colors <- scales::hue_pal()(n_labels)
  }
  
  names(bed_colors) <- unique_labels
  
  # Print the color mapping for reference
  cat("\nBED region color mapping:\n")
  for (i in 1:length(bed_colors)) {
    cat(names(bed_colors)[i], ":", bed_colors[i], "\n")
  }
}

# Read the data
window_size <- '100kb'
num_haplo <- 466
num_sample <- 234
data <- read_tsv(paste0("/home/guarracino/Desktop/hprc25272.CHM13.w", window_size, "-xm5-id095-l5k.tsv.gz"))

# Parse the chroms-num_haplotypes column to extract chromosome information
parse_chroms_column <- function(chroms_str) {
  if (is.na(chroms_str) || chroms_str == "") {
    return(list(num_chromosomes = 0, chromosomes = NA))
  }
  
  # Split by comma
  pairs <- strsplit(chroms_str, ",")[[1]]
  
  # Extract chromosome names (everything before the last hyphen)
  chrom_names <- sapply(pairs, function(x) {
    parts <- strsplit(x, "-")[[1]]
    # Join all but the last part (in case chromosome name contains hyphen)
    paste(parts[-length(parts)], collapse = "-")
  })
  
  # Remove any whitespace
  chrom_names <- trimws(chrom_names)
  
  return(list(
    num_chromosomes = length(unique(chrom_names)),
    chromosomes = paste(unique(chrom_names), collapse = ",")
  ))
}

# Apply parsing to the entire dataset
parsed_data <- lapply(data$`chroms-num_haplotypes`, parse_chroms_column)
data$num_chromosomes <- sapply(parsed_data, function(x) x$num_chromosomes)
data$chromosomes <- sapply(parsed_data, function(x) x$chromosomes)

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

# If BED regions exist, ensure they have the same factor levels
if (!is.null(bed_regions)) {
  bed_regions$chromosome <- factor(bed_regions$chromosome, 
                                   levels = levels(data$chromosome))
}

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
                                 labels = c("Alignments not merged", "Alignments"))

# Hide not merged alignments
data_alignments <- data_alignments %>% filter(metric == "Alignments")

# Create combined alignments plot with BED regions
p_combined_alignments <- ggplot(data_alignments, aes(x = position, y = value, color = metric, group = metric))

# Add BED regions if available (use log scale appropriate limits)
if (!is.null(bed_regions)) {
  p_combined_alignments <- p_combined_alignments +
    geom_rect(data = bed_regions,
              aes(xmin = start_mbp, xmax = end_mbp, ymin = 0.1, ymax = 1e8, fill = name),
              inherit.aes = FALSE,
              #fill = "grey50",
              alpha = 0.2) +
    scale_fill_manual(values = bed_colors, name = "Region")
}

p_combined_alignments <- p_combined_alignments +
  geom_line() +
  facet_wrap(~ chromosome, scales = "free", ncol = 5) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_color_manual(values = c("black", "red")) +
  guides(color = "none") +  # To hide the Metric legend
  labs(
    title = paste0("Alignments across CHM13 (", window_size, " windows)"),
    x = "Position (Mbp)",
    y = "Count (log10)",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 13),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 13),
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

# Create combined haplotypes and samples plot with BED regions
p_combined_haplo_samples <- ggplot(data_haplo_samples, aes(x = position, y = value, color = metric, group = metric))

# Add BED regions if available
if (!is.null(bed_regions)) {
  p_combined_haplo_samples <- p_combined_haplo_samples +
    geom_rect(data = bed_regions,
              aes(xmin = start_mbp, xmax = end_mbp, ymin = 0, ymax = 500, fill = name),
              inherit.aes = FALSE,
              #fill = "grey50",
              alpha = 0.2) +
    scale_fill_manual(values = bed_colors, name = "Region")
}

p_combined_haplo_samples <- p_combined_haplo_samples +
  geom_line() +
  facet_wrap(~ chromosome, scales = "free", ncol = 5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(0, NA)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6),
                     expand = expansion(mult = c(0.09, 0.15))) +
  scale_color_manual(values = c("darkgreen", "purple")) +
  labs(
    title = paste0("Samples/Haplotypes across CHM13 (", window_size, " windows)"),
    x = "Position (Mbp)",
    y = "Count",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 13),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 13),
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
            inherit.aes = FALSE,
            color = "black", hjust = +0.9, vjust = +1.5, size = 3)

# Create p_num_chromosomes plot with BED regions
p_num_chromosomes <- ggplot(data, aes(x = position, y = num_chromosomes))

# Add BED regions if available
if (!is.null(bed_regions)) {
  p_num_chromosomes <- p_num_chromosomes +
    geom_rect(data = bed_regions,
              aes(xmin = start_mbp, xmax = end_mbp, ymin = 0, ymax = 25, fill = name),
              inherit.aes = FALSE,
              #fill = "grey50",
              alpha = 0.2) +
    scale_fill_manual(values = bed_colors, name = "Region")
}

p_num_chromosomes <- p_num_chromosomes +
  geom_line(color = "darkorange", linewidth = 0.8) +
  facet_wrap(~ chromosome, scales = "free_x", ncol = 5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(0, NA)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6),
                     expand = expansion(mult = c(0.09, 0.15))) +
  labs(
    title = paste0("Number of unique chromosomes per region across CHM13 (", window_size, "windows)"),
    x = "Position (Mbp)",
    y = "Number of chromosomes"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 13),
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 15, face = "bold")
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7)

# Create p_num_chromosomes_wide plot with BED regions
p_num_chromosomes_wide <- ggplot(data, aes(x = position, y = num_chromosomes))

# Add BED regions if available (as horizontal bands for wide format)
if (!is.null(bed_regions)) {
  p_num_chromosomes_wide <- p_num_chromosomes_wide +
    geom_rect(data = bed_regions,
              aes(xmin = start_mbp, xmax = end_mbp, ymin = 0, ymax = 25, fill = name),
              inherit.aes = FALSE,
              #fill = "grey50",
              alpha = 0.2) +
    scale_fill_manual(values = bed_colors, name = "Region")
}

p_num_chromosomes_wide <- p_num_chromosomes_wide +
  geom_line(color = "darkorange", linewidth = 0.5, alpha = 0.9) +
  facet_grid(chromosome ~ ., scales = "free_x", switch = "y") +
  
  scale_x_continuous(
    breaks = seq(0, 300, 50),
    expand = c(0.01, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 25),
    breaks = seq(0, 25, 5)
  ) +
  
  # Add horizontal reference lines
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.5, linewidth = 0.3) +
  
  labs(
    title = paste0("Number of unique chromosomes per region across CHM13 (", window_size, " windows)"),
    x = "Position (Mbp)",
    y = "Number of chromosomes"
  ) +
  
  theme_minimal() +
  theme(
    strip.text.y = element_text(size = 8, angle = 0),
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.placement = "outside",
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 8, hjust = 0, color = "gray40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2),
    panel.spacing = unit(0.05, "lines"),
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold")
  )

# Print summary of BED regions if loaded
if (!is.null(bed_regions)) {
  cat("\nBED regions loaded successfully!\n")
  cat("Number of regions:", nrow(bed_regions), "\n")
  cat("Chromosomes covered:", paste(unique(bed_regions$chromosome), collapse = ", "), "\n\n")
} else {
  cat("\nNo BED file loaded. To add region annotations, set bed_file_path to your BED file.\n\n")
}

# Print the plots
print(p_combined_alignments)
print(p_combined_haplo_samples)
print(p_num_chromosomes)
print(p_num_chromosomes_wide)

width=16
height=9
dpi=300

# Save the plots
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
ggsave(
  filename = "p_num_chromosomes_wide.pdf",
  plot = p_num_chromosomes_wide,
  width = width,
  height = height,
  dpi = dpi,
  units = "in",
  bg = "white"
)

# Function to create a single chromosome plot with specified range
plot_single_chromosome <- function(data, 
                                   target_chr = "chr1", 
                                   start_mbp = NULL, 
                                   end_mbp = NULL,
                                   bed_regions = NULL) {
  
  # Filter for the specified chromosome
  data_chr <- data %>%
    filter(chromosome == target_chr)
  
  # If no data for this chromosome, return NULL with warning
  if (nrow(data_chr) == 0) {
    warning(paste("No data found for chromosome:", target_chr))
    return(NULL)
  }
  
  # Get the actual data range for this chromosome
  data_min <- min(data_chr$position)
  data_max <- max(data_chr$position)
  
  # Adjust start_mbp and end_mbp to actual data range
  if (!is.null(start_mbp)) {
    # If start_mbp is less than minimum, use minimum
    if (start_mbp < data_min) {
      message(paste0("Adjusting start from ", start_mbp, " to data minimum ", round(data_min, 2), " Mbp"))
      start_mbp <- data_min
    }
  } else {
    start_mbp <- data_min
  }
  
  if (!is.null(end_mbp)) {
    # If end_mbp exceeds maximum, use maximum
    if (end_mbp > data_max) {
      message(paste0("Adjusting end from ", end_mbp, " to data maximum ", round(data_max, 2), " Mbp"))
      end_mbp <- data_max
    }
  } else {
    end_mbp <- data_max
  }
  
  # Apply position range filter
  data_chr <- data_chr %>%
    filter(position >= start_mbp & position <= end_mbp)
  
  # Filter BED regions if they exist
  if (!is.null(bed_regions)) {
    bed_regions <- bed_regions %>%
      filter(chromosome == target_chr) %>%
      filter(end_mbp >= start_mbp & start_mbp <= end_mbp) %>%
      # Clip BED regions to the specified range
      mutate(
        start_mbp = pmax(start_mbp, !!start_mbp),
        end_mbp = pmin(end_mbp, !!end_mbp)
      )
  }
  
  # Create the plot
  p <- ggplot(data_chr, aes(x = position, y = num_chromosomes))
  
  # Add BED regions if available with colors
  if (!is.null(bed_regions) && nrow(bed_regions) > 0) {
    p <- p +
      geom_rect(data = bed_regions,
                aes(xmin = start_mbp, xmax = end_mbp, ymin = 0, ymax = 25, fill = name),
                inherit.aes = FALSE,
                alpha = 0.3)
    
    # Add color scale if bed_colors provided
    if (exists("bed_colors")) {
      p <- p + scale_fill_manual(values = bed_colors, name = "Region")
    }
  }
  
  p <- p +
    geom_line(color = "darkorange", linewidth = 0.8, alpha = 0.9) +
    
    # Set x-axis limits based on specified range
    scale_x_continuous(
      limits = c(start_mbp, end_mbp),
      breaks = scales::pretty_breaks(n = 10),
      expand = c(0.01, 0)
    ) +
    
    scale_y_continuous(
      limits = c(0, 25),
      breaks = seq(0, 25, 5)
    ) +
    
    # Add horizontal reference line at 1
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.5, linewidth = 0.3) +
    
    labs(
      title = paste0("Number of Unique Chromosomes per Region - ", target_chr, 
                     " (", round(start_mbp, 1), "-", round(end_mbp, 1), " Mbp)"),
      x = "Position (Mbp)",
      y = "Number of Chromosomes"
    ) +
    
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2),
      legend.position = "right"
    )
  
  return(p)
}

# General function to create chromosome matching heatmap
plot_chromosome_matching_heatmap <- function(data, 
                                             target_chr = "chr1", 
                                             start_mbp = NULL, 
                                             end_mbp = NULL,
                                             bed_regions = NULL,
                                             highlight_regions = TRUE,
                                             output_filename = NULL) {
  
  # Filter for the target chromosome
  chr_data <- data %>%
    filter(chromosome == target_chr) %>%
    mutate(position = (start + end) / 2 / 1e6)
  
  # Check if data exists for this chromosome
  if (nrow(chr_data) == 0) {
    warning(paste("No data found for chromosome:", target_chr))
    return(NULL)
  }
  
  # Get actual data range
  data_min <- min(chr_data$position)
  data_max <- max(chr_data$position)
  
  # Set position range with defaults
  if (is.null(start_mbp)) {
    start_mbp <- data_min
  } else {
    # Ensure start_mbp is within data range
    if (start_mbp < data_min) {
      message(paste0("Adjusting start from ", start_mbp, " to data minimum ", round(data_min, 2), " Mbp"))
      start_mbp <- data_min
    }
    if (start_mbp > data_max) {
      warning(paste0("start_mbp (", start_mbp, ") is beyond data range. Using data minimum."))
      start_mbp <- data_min
    }
  }
  
  if (is.null(end_mbp)) {
    end_mbp <- data_max
  } else {
    # Ensure end_mbp is within data range
    if (end_mbp > data_max) {
      message(paste0("Adjusting end from ", end_mbp, " to data maximum ", round(data_max, 2), " Mbp"))
      end_mbp <- data_max
    }
    if (end_mbp < data_min) {
      warning(paste0("end_mbp (", end_mbp, ") is before data range. Using data maximum."))
      end_mbp <- data_max
    }
  }
  
  # Validate that start < end
  if (start_mbp >= end_mbp) {
    warning("start_mbp must be less than end_mbp. Swapping values.")
    temp <- start_mbp
    start_mbp <- end_mbp
    end_mbp <- temp
  }
  
  # Filter by position range
  chr_data_filtered <- chr_data %>%
    filter(position >= start_mbp & position <= end_mbp)
  
  if (nrow(chr_data_filtered) == 0) {
    warning(paste0("No data in specified range [", round(start_mbp, 2), ", ", round(end_mbp, 2), "] Mbp"))
    return(NULL)
  }
  
  # Parse the comma-separated chromosomes list into a binary matrix
  parse_chromosomes_binary <- function(data) {
    # Get all possible chromosomes
    all_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
    
    # Initialize binary matrix
    binary_matrix <- matrix(0, 
                            nrow = nrow(data), 
                            ncol = length(all_chroms),
                            dimnames = list(NULL, all_chroms))
    
    # Fill in the binary matrix
    for(i in 1:nrow(data)) {
      if(!is.na(data$chromosomes[i])) {
        # Split the comma-separated list
        chroms_list <- strsplit(data$chromosomes[i], ",")[[1]]
        # Mark present chromosomes as 1
        for(chr in chroms_list) {
          chr_clean <- trimws(chr)
          if(chr_clean %in% all_chroms) {
            binary_matrix[i, chr_clean] <- 1
          }
        }
      }
    }
    
    return(binary_matrix)
  }
  
  # Create binary matrix
  chr_binary <- parse_chromosomes_binary(chr_data_filtered)
  
  # Combine with position data for plotting
  chr_heatmap_data <- cbind(chr_data_filtered %>% select(position), chr_binary) %>%
    as.data.frame() %>%
    pivot_longer(cols = -position, names_to = "matching_chromosome", values_to = "present")
  
  # Calculate appropriate x-axis breaks based on the range
  x_range <- end_mbp - start_mbp
  if (x_range <= 5) {
    x_breaks <- seq(floor(start_mbp*2)/2, ceiling(end_mbp*2)/2, by = 0.5)
  } else if (x_range <= 10) {
    x_breaks <- seq(floor(start_mbp), ceiling(end_mbp), by = 1)
  } else if (x_range <= 50) {
    x_breaks <- seq(floor(start_mbp/5)*5, ceiling(end_mbp/5)*5, by = 5)
  } else if (x_range <= 100) {
    x_breaks <- seq(floor(start_mbp/10)*10, ceiling(end_mbp/10)*10, by = 10)
  } else {
    x_breaks <- seq(floor(start_mbp/50)*50, ceiling(end_mbp/50)*50, by = 50)
  }
  
  # Filter to ensure breaks are within the actual range
  x_breaks <- x_breaks[x_breaks >= start_mbp & x_breaks <= end_mbp]
  
  # Create the binary heatmap
  p_heatmap <- ggplot(chr_heatmap_data, 
                      aes(x = position, y = matching_chromosome, fill = factor(present))) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "white", "1" = "darkblue"),
                      labels = c("No match", "Match"),
                      name = "Status") +
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(start_mbp, end_mbp),
                       breaks = x_breaks) +
    scale_y_discrete(limits = rev(paste0("chr", c(1:22, "X", "Y", "M")))) +
    labs(
      title = paste0("Chromosome Matching Pattern for ", target_chr, 
                     " (", round(start_mbp, 1), "-", round(end_mbp, 1), " Mbp)"),
      subtitle = paste0(window_size, " windows (", nrow(chr_data_filtered), " windows in range)"),
      x = paste0("Position on ", target_chr, " (Mbp)"),
      y = "Matching Chromosome"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "right"
    )
  
  # Add BED region annotations if provided
  if (!is.null(bed_regions) && highlight_regions) {
    # Store the input parameters with different names to avoid naming conflicts
    input_start <- start_mbp
    input_end <- end_mbp
    
    chr_bed <- bed_regions %>% 
      filter(chromosome == target_chr) %>%
      # Filter for regions that overlap with the specified range
      filter(end_mbp >= input_start & start_mbp <= input_end) %>%
      # Clip the regions to the specified range
      mutate(
        start_mbp_clipped = pmax(start_mbp, input_start),
        end_mbp_clipped = pmin(end_mbp, input_end),
        mid_pos = (start_mbp_clipped + end_mbp_clipped) / 2
      )
    
    if (nrow(chr_bed) > 0) {
      # Add vertical lines for region boundaries
      p_heatmap <- p_heatmap +
        geom_vline(data = chr_bed, aes(xintercept = start_mbp_clipped), 
                   color = "red", linetype = "dashed", alpha = 0.5, size = 0.3) +
        geom_vline(data = chr_bed, aes(xintercept = end_mbp_clipped), 
                   color = "red", linetype = "dashed", alpha = 0.5, size = 0.3)
      
      # Add region labels at the top
      # Calculate y position for labels (above the heatmap)
      y_label_pos <- length(paste0("chr", c(1:22, "X", "Y", "M"))) + 0.5
      
      p_heatmap <- p_heatmap +
        geom_text(data = chr_bed,
                  aes(x = mid_pos, y = y_label_pos, label = name),
                  angle = 0, vjust = 0, size = 3, color = "red",
                  inherit.aes = FALSE)
    }
  }
  
  # Print summary statistics
  cat("\n=== Chromosome Matching Summary for", target_chr, "===\n")
  cat("Position range:", round(start_mbp, 2), "-", round(end_mbp, 2), "Mbp\n")
  cat("Number of windows:", nrow(chr_data_filtered), "\n")
  
  # Calculate matching frequency
  match_freq <- colSums(chr_binary)
  match_freq <- match_freq[match_freq > 0]
  match_freq <- sort(match_freq, decreasing = TRUE)
  
  cat("\nTop matching chromosomes:\n")
  if (length(match_freq) > 0) {
    top_matches <- head(match_freq, 10)
    for(i in 1:length(top_matches)) {
      cat(sprintf("  %s: %d windows (%.1f%%)\n", 
                  names(top_matches)[i], 
                  top_matches[i], 
                  100 * top_matches[i] / nrow(chr_data_filtered)))
    }
  }
  
  # If BED regions exist, summarize matches within them
  if (!is.null(bed_regions) && highlight_regions) {
    input_start <- start_mbp
    input_end <- end_mbp
    
    chr_bed <- bed_regions %>% 
      filter(chromosome == target_chr) %>%
      filter(end_mbp >= input_start & start_mbp <= input_end)
    
    if (nrow(chr_bed) > 0) {
      cat("\n--- Matches within annotated regions ---\n")
      for(i in 1:nrow(chr_bed)) {
        # Use the actual region boundaries (not clipped) for data filtering
        region_start <- max(chr_bed$start_mbp[i], input_start)
        region_end <- min(chr_bed$end_mbp[i], input_end)
        
        region_data <- chr_data_filtered %>%
          filter(position >= region_start & position <= region_end)
        
        if(nrow(region_data) > 0) {
          all_chroms_in_region <- unlist(strsplit(region_data$chromosomes, ","))
          all_chroms_in_region <- trimws(all_chroms_in_region)
          all_chroms_in_region <- all_chroms_in_region[!is.na(all_chroms_in_region)]
          
          if(length(all_chroms_in_region) > 0) {
            chrom_freq_region <- table(all_chroms_in_region)
            chrom_freq_region <- sort(chrom_freq_region, decreasing = TRUE)
            
            cat("\n", chr_bed$name[i], " (", round(region_start, 2), 
                "-", round(region_end, 2), " Mbp):\n", sep="")
            
            top_in_region <- head(chrom_freq_region, 5)
            for(j in 1:length(top_in_region)) {
              cat(sprintf("    %s: %d occurrences\n", 
                          names(top_in_region)[j], 
                          top_in_region[j]))
            }
          }
        }
      }
    }
  }

  return(p_heatmap)
}

plot_single_chromosome(data, 
                       target_chr = "chrY", 
                       start_mbp = 0, 
                       end_mbp = 999,
                       bed_regions = bed_regions)

plot_chromosome_matching_heatmap(data, 
                                 target_chr = "chrM",
                                 #start_mbp = 11,
                                 #end_mbp = 12,
                                 bed_regions = bed_regions)


# Loop through all chromosomes and save heatmaps
for (chr in paste0("chr", c(1:22, "X", "Y", "M"))) {
  p <- plot_chromosome_matching_heatmap(data, 
                                        target_chr = chr,
                                        bed_regions = bed_regions,
                                        save_plot = TRUE)
  if (!is.null(p)) print(p)
}


plot_chromosome_matching_heatmap2 <- function(data, 
                                             target_chr = "chr1", 
                                             start_mbp = NULL, 
                                             end_mbp = NULL,
                                             bed_regions = NULL,
                                             highlight_regions = TRUE,
                                             show_numbers = "auto",  # "auto", "always", "never"
                                             number_threshold = 30,  # max windows to show numbers
                                             output_filename = NULL) {
  
  # Filter for the target chromosome
  chr_data <- data %>%
    filter(chromosome == target_chr) %>%
    mutate(position = (start + end) / 2 / 1e6)
  
  # Check if data exists for this chromosome
  if (nrow(chr_data) == 0) {
    warning(paste("No data found for chromosome:", target_chr))
    return(NULL)
  }
  
  # Get actual data range
  data_min <- min(chr_data$position)
  data_max <- max(chr_data$position)
  
  # Set position range with defaults
  if (is.null(start_mbp)) {
    start_mbp <- data_min
  } else {
    # Ensure start_mbp is within data range
    if (start_mbp < data_min) {
      message(paste0("Adjusting start from ", start_mbp, " to data minimum ", round(data_min, 2), " Mbp"))
      start_mbp <- data_min
    }
    if (start_mbp > data_max) {
      warning(paste0("start_mbp (", start_mbp, ") is beyond data range. Using data minimum."))
      start_mbp <- data_min
    }
  }
  
  if (is.null(end_mbp)) {
    end_mbp <- data_max
  } else {
    # Ensure end_mbp is within data range
    if (end_mbp > data_max) {
      message(paste0("Adjusting end from ", end_mbp, " to data maximum ", round(data_max, 2), " Mbp"))
      end_mbp <- data_max
    }
    if (end_mbp < data_min) {
      warning(paste0("end_mbp (", end_mbp, ") is before data range. Using data maximum."))
      end_mbp <- data_max
    }
  }
  
  # Validate that start < end
  if (start_mbp >= end_mbp) {
    warning("start_mbp must be less than end_mbp. Swapping values.")
    temp <- start_mbp
    start_mbp <- end_mbp
    end_mbp <- temp
  }
  
  # Filter by position range
  chr_data_filtered <- chr_data %>%
    filter(position >= start_mbp & position <= end_mbp)
  
  if (nrow(chr_data_filtered) == 0) {
    warning(paste0("No data in specified range [", round(start_mbp, 2), ", ", round(end_mbp, 2), "] Mbp"))
    return(NULL)
  }
  
  # Parse the chroms-num_haplotypes column to get actual counts
  parse_chromosomes_counts <- function(data) {
    # Get all possible chromosomes
    all_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
    
    # Initialize count matrix
    count_matrix <- matrix(0, 
                           nrow = nrow(data), 
                           ncol = length(all_chroms),
                           dimnames = list(NULL, all_chroms))
    
    # Fill in the count matrix
    for(i in 1:nrow(data)) {
      if(!is.na(data$`chroms-num_haplotypes`[i]) && data$`chroms-num_haplotypes`[i] != "") {
        # Split the comma-separated list
        pairs <- strsplit(data$`chroms-num_haplotypes`[i], ",")[[1]]
        
        # Parse each pair to get chromosome and count
        for(pair in pairs) {
          # Split by the last hyphen to separate chromosome from count
          parts <- strsplit(trimws(pair), "-")[[1]]
          if(length(parts) >= 2) {
            # The chromosome is everything except the last part
            chr_name <- paste(parts[-length(parts)], collapse = "-")
            # The count is the last part
            count <- as.numeric(parts[length(parts)])
            
            if(chr_name %in% all_chroms && !is.na(count)) {
              count_matrix[i, chr_name] <- count
            }
          }
        }
      }
    }
    
    return(count_matrix)
  }
  
  # Create count matrix
  chr_counts <- parse_chromosomes_counts(chr_data_filtered)
  
  # Combine with position data for plotting
  chr_heatmap_data <- cbind(chr_data_filtered %>% select(position), chr_counts) %>%
    as.data.frame() %>%
    pivot_longer(cols = -position, names_to = "matching_chromosome", values_to = "count")
  
  # Determine whether to show numbers
  show_text <- FALSE
  if (show_numbers == "always") {
    show_text <- TRUE
  } else if (show_numbers == "auto") {
    # Show numbers if we have 30 or fewer windows
    show_text <- nrow(chr_data_filtered) <= number_threshold
  }
  
  # Calculate appropriate x-axis breaks based on the range
  x_range <- end_mbp - start_mbp
  if (x_range <= 5) {
    x_breaks <- seq(floor(start_mbp*2)/2, ceiling(end_mbp*2)/2, by = 0.5)
  } else if (x_range <= 10) {
    x_breaks <- seq(floor(start_mbp), ceiling(end_mbp), by = 1)
  } else if (x_range <= 50) {
    x_breaks <- seq(floor(start_mbp/5)*5, ceiling(end_mbp/5)*5, by = 5)
  } else if (x_range <= 100) {
    x_breaks <- seq(floor(start_mbp/10)*10, ceiling(end_mbp/10)*10, by = 10)
  } else {
    x_breaks <- seq(floor(start_mbp/50)*50, ceiling(end_mbp/50)*50, by = 50)
  }
  
  # Filter to ensure breaks are within the actual range
  x_breaks <- x_breaks[x_breaks >= start_mbp & x_breaks <= end_mbp]
  
  # Get the maximum count for color scale
  max_count <- max(chr_heatmap_data$count, na.rm = TRUE)
  
  # Create custom color scale
  # White for 0, then a gradient from a visible blue to dark blue
  if (max_count > 0) {
    # Create breaks for the color scale
    if (max_count == 1) {
      color_breaks <- c(0, 1)
      color_values <- c(0, 1)
      colors <- c("white", "#2166ac")
    } else {
      # Create a color scale that starts at 0 (white) and has visible colors for values > 0
      color_breaks <- c(0, 1, max_count)
      color_values <- c(0, 0.4, 1)  # 0 maps to white, 1 maps to 40% of the gradient, max to 100%
      colors <- c("white", "#6baed6", "#08306b")
    }
  } else {
    colors <- c("white", "white")
    color_breaks <- c(0, 1)
    color_values <- c(0, 1)
  }
  
  # Create the heatmap
  p_heatmap <- ggplot(chr_heatmap_data, 
                      aes(x = position, y = matching_chromosome, fill = count))
  
  # Add the tile layer
  p_heatmap <- p_heatmap + geom_tile()
  
  # Add text labels if appropriate
  if (show_text) {
    # Create a subset with only non-zero values for text
    text_data <- chr_heatmap_data %>%
      filter(count > 0)
    
    p_heatmap <- p_heatmap +
      geom_text(data = text_data,
                aes(x = position, y = matching_chromosome, label = count),
                color = "white", size = 2.5, fontface = "bold")
  }
  
  # Apply color scale
  if (max_count > 0) {
    p_heatmap <- p_heatmap +
      scale_fill_gradientn(
        colors = colors,
        values = color_values,
        breaks = pretty(c(0, max_count), n = 5),
        limits = c(0, max_count),
        name = "Haplotype\nCount",
        na.value = "white"
      )
  } else {
    p_heatmap <- p_heatmap +
      scale_fill_gradient(low = "white", high = "white", 
                          name = "Haplotype\nCount",
                          limits = c(0, 1))
  }
  
  p_heatmap <- p_heatmap +
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(start_mbp, end_mbp),
                       breaks = x_breaks) +
    scale_y_discrete(limits = rev(paste0("chr", c(1:22, "X", "Y", "M")))) +
    labs(
      title = paste0("Chromosome Matching Pattern for ", target_chr, 
                     " (", round(start_mbp, 1), "-", round(end_mbp, 1), " Mbp)"),
      subtitle = paste0(window_size, " windows (", nrow(chr_data_filtered), 
                        " windows, max ", max_count, " haplotypes",
                        ifelse(show_text, ", showing counts", ""), ")"),
      x = paste0("Position on ", target_chr, " (Mbp)"),
      y = "Matching Chromosome"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "right",
      legend.key.height = unit(1, "cm"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  # Add BED region annotations if provided
  if (!is.null(bed_regions) && highlight_regions) {
    # Store the input parameters with different names to avoid naming conflicts
    input_start <- start_mbp
    input_end <- end_mbp
    
    chr_bed <- bed_regions %>% 
      filter(chromosome == target_chr) %>%
      # Filter for regions that overlap with the specified range
      filter(end_mbp >= input_start & start_mbp <= input_end) %>%
      # Clip the regions to the specified range
      mutate(
        start_mbp_clipped = pmax(start_mbp, input_start),
        end_mbp_clipped = pmin(end_mbp, input_end),
        mid_pos = (start_mbp_clipped + end_mbp_clipped) / 2
      )
    
    if (nrow(chr_bed) > 0) {
      # Add vertical lines for region boundaries
      p_heatmap <- p_heatmap +
        geom_vline(data = chr_bed, aes(xintercept = start_mbp_clipped), 
                   color = "red", linetype = "dashed", alpha = 0.7, size = 0.5) +
        geom_vline(data = chr_bed, aes(xintercept = end_mbp_clipped), 
                   color = "red", linetype = "dashed", alpha = 0.7, size = 0.5)
      
      # Add region labels at the top
      # Calculate y position for labels (above the heatmap)
      y_label_pos <- length(paste0("chr", c(1:22, "X", "Y", "M"))) + 0.5
      
      p_heatmap <- p_heatmap +
        geom_text(data = chr_bed,
                  aes(x = mid_pos, y = y_label_pos, label = name),
                  angle = 0, vjust = 0, size = 3, color = "red",
                  inherit.aes = FALSE)
    }
  }
  
  # Print summary statistics
  cat("\n=== Chromosome Matching Summary for", target_chr, "===\n")
  cat("Position range:", round(start_mbp, 2), "-", round(end_mbp, 2), "Mbp\n")
  cat("Number of windows:", nrow(chr_data_filtered), "\n")
  cat("Maximum haplotype count:", max_count, "\n")
  if (show_text) {
    cat("Showing count numbers in cells\n")
  }
  
  # Calculate total haplotypes per chromosome
  haplo_totals <- colSums(chr_counts)
  haplo_totals <- haplo_totals[haplo_totals > 0]
  haplo_totals <- sort(haplo_totals, decreasing = TRUE)
  
  cat("\nTop matching chromosomes (by total haplotype count):\n")
  if (length(haplo_totals) > 0) {
    top_matches <- head(haplo_totals, 10)
    for(i in 1:length(top_matches)) {
      # Also calculate how many windows have matches
      windows_with_match <- sum(chr_counts[, names(top_matches)[i]] > 0)
      avg_haplo <- round(top_matches[i] / windows_with_match, 1)
      
      cat(sprintf("  %s: %d total haplotypes across %d windows (avg %.1f per window)\n", 
                  names(top_matches)[i], 
                  top_matches[i],
                  windows_with_match,
                  avg_haplo))
    }
  }
  
  # If BED regions exist, summarize matches within them
  if (!is.null(bed_regions) && highlight_regions) {
    input_start <- start_mbp
    input_end <- end_mbp
    
    chr_bed <- bed_regions %>% 
      filter(chromosome == target_chr) %>%
      filter(end_mbp >= input_start & start_mbp <= input_end)
    
    if (nrow(chr_bed) > 0) {
      cat("\n--- Matches within annotated regions ---\n")
      for(i in 1:nrow(chr_bed)) {
        # Use the actual region boundaries (not clipped) for data filtering
        region_start <- max(chr_bed$start_mbp[i], input_start)
        region_end <- min(chr_bed$end_mbp[i], input_end)
        
        region_indices <- which(chr_data_filtered$position >= region_start & 
                                  chr_data_filtered$position <= region_end)
        
        if(length(region_indices) > 0) {
          # Get counts for this region
          region_counts <- chr_counts[region_indices, , drop = FALSE]
          region_totals <- colSums(region_counts)
          region_totals <- region_totals[region_totals > 0]
          region_totals <- sort(region_totals, decreasing = TRUE)
          
          cat("\n", chr_bed$name[i], " (", round(region_start, 2), 
              "-", round(region_end, 2), " Mbp, ", 
              length(region_indices), " windows):\n", sep="")
          
          if(length(region_totals) > 0) {
            top_in_region <- head(region_totals, 5)
            for(j in 1:length(top_in_region)) {
              windows_with_match <- sum(region_counts[, names(top_in_region)[j]] > 0)
              cat(sprintf("    %s: %d haplotypes in %d windows\n", 
                          names(top_in_region)[j], 
                          top_in_region[j],
                          windows_with_match))
            }
          } else {
            cat("    No matches in this region\n")
          }
        }
      }
    }
  }
  
  return(p_heatmap)
}

# Full chromosome with gradient (no numbers shown due to many windows)
plot_chromosome_matching_heatmap2(data, 
                                 target_chr = "chrM",
                                 bed_regions = NULL)

# Zoomed region with automatic number display (if <= 30 windows)
plot_chromosome_matching_heatmap2(data, 
                                 target_chr = "chrY",
                                 start_mbp = 10,
                                 end_mbp = 13,
                                 bed_regions = bed_regions)

# Force showing numbers regardless of window count
plot_chromosome_matching_heatmap2(data, 
                                 target_chr = "chr1",
                                 start_mbp = 120,
                                 end_mbp = 125,
                                 bed_regions = bed_regions,
                                 show_numbers = "always")

# Never show numbers even when zoomed in
plot_chromosome_matching_heatmap2(data, 
                                 target_chr = "chr1",
                                 start_mbp = 120,
                                 end_mbp = 122,
                                 bed_regions = bed_regions,
                                 show_numbers = "never")
