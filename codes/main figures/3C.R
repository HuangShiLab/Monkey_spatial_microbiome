#######################################
# Monkey Study
#
# 5 most common phylum, stacked bar plot, mocosal and content
#
# Author: HOU Shuwen
#######################################

# Load required packages
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(tidyr)

#--------------------------------------#
# Data Import (Mucosa)
#--------------------------------------#

# Set working directory for mucosa data
setwd("~/Downloads/monkey_mucosa")

# Read genus abundance data and metadata
mucosa_genus <- as.data.frame(read.table("genus_abundance.txt", header = TRUE))
metadata <- read.table("meta.txt", sep="\t", header=TRUE)
ID <- metadata$ID  # Extract sample IDs

# Aggregate data by Phylum
mucosa_genus <- mucosa_genus %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  as.data.frame()

# Format data: set row names and convert to numeric
rownames(mucosa_genus) <- mucosa_genus$Phylum
mucosa_genus <- mucosa_genus[, -1]
mucosa_genus[] <- lapply(mucosa_genus, function(x) as.numeric(as.character(x)))


#--------------------------------------#
# Data Import (Content)
#--------------------------------------#

# Set working directory for content data
setwd("~/Downloads/monkey_data")

# Read genus abundance data and metadata
content_genus <- as.data.frame(read.table("2b_genus_abundance.txt", header = TRUE))
metadata <- read.table("kraken_metadata.txt", sep="\t", header=TRUE)
ID <- metadata$ID  # Extract sample IDs

# Aggregate data by Phylum
content_genus <- content_genus %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  as.data.frame()

# Format data: set row names and convert to numeric
rownames(content_genus) <- content_genus$Phylum
content_genus <- content_genus[, -1]
content_genus[] <- lapply(content_genus, function(x) as.numeric(as.character(x)))


#--------------------------------------#
# Select Top 5 Phylum for Each Sample
#--------------------------------------#

get_top_n_rows <- function(df, n_top = 5) {
  # Ensure row names are retained by creating a Genus column
  df$Genus <- rownames(df)
  
  # Extract the top n genera for each column
  top_rows <- unique(do.call(rbind, lapply(names(df)[-ncol(df)], function(col) {
    # Sort by descending order of the column and select top n rows
    ordered_df <- df[order(-df[[col]]), ]
    ordered_df[1:n_top, ]
  })))
  
  # Restore row names and remove the temporary Genus column
  rownames(top_rows) <- top_rows$Genus
  top_rows$Genus <- NULL
  
  # Calculate the sum of the remaining genera for each column
  # Exclude the Genus column when summing
  others_sum <- colSums(df[!rownames(df) %in% rownames(top_rows), -ncol(df)])
  others_row <- as.data.frame(t(others_sum))
  rownames(others_row) <- "others"
  
  # Append the 'others' row to the top_rows data frame
  top_rows_with_others <- rbind(top_rows, others_row)
  
  return(top_rows_with_others)
}

# Process mucosa data
# Normalize values so each row sums to 1
# Check sum to confirm normalization

top_10_mucosa_genus <- get_top_n_rows(mucosa_genus, n_top = 5) %>%
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
top_10_mucosa_genus %>% summarise(across(everything(), sum, na.rm = TRUE))

# Process content data
# Normalize values so each row sums to 1
# Check sum to confirm normalization

top_10_content_genus <- get_top_n_rows(content_genus, n_top = 5) %>%
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
top_10_content_genus %>% summarise(across(everything(), sum, na.rm = TRUE))

#--------------------------------------#
# Merge Data from Both Sources
#--------------------------------------#

# Convert row names to a column for merging
top_10_content_genus <- top_10_content_genus %>% tibble::rownames_to_column("name")
top_10_mucosa_genus <- top_10_mucosa_genus %>% tibble::rownames_to_column("name")

# Merge both datasets and replace NA with 0
all_genera <- full_join(top_10_content_genus, top_10_mucosa_genus, by = "name")
all_genera[is.na(all_genera)] <- 0

# Restore row names and remove extra column
rownames(all_genera) <- all_genera[, 1]
all_genera <- all_genera[, -1]

#--------------------------------------#
# Merge sub-phylums
#--------------------------------------#

# Define the phyla you want to merge
to_merge <- c("Bacillota", "Bacillota_A", "Bacillota_C", "Firmicutes", "Firmicutes_A", "Firmicutes_C")

# Sum rows corresponding to Bacillota-related groups
firmicutes_sum <- colSums(all_genera[rownames(all_genera) %in% to_merge, , drop = FALSE])

# Remove the original rows
all_genera <- all_genera[!rownames(all_genera) %in% to_merge, , drop = FALSE]

# Add the new merged row
all_genera["Firmicutes", ] <- firmicutes_sum

#--------------------------------------#
# Data Cleaning and Normalization
#--------------------------------------#

# Transpose for plotting
all_genera <- as.data.frame(t(all_genera))

# Remove the 'human' column if it exists
all_genera <- all_genera %>% select(-human)

# Remove low-abundance genera (max < 0.05)
columns_to_remove <- all_genera %>%
  summarise(across(everything(), max, na.rm = TRUE)) %>%
  select(where(~ . < 0.05)) %>%
  colnames()
all_genera <- all_genera %>% select(-all_of(columns_to_remove))

# Normalize each row so values sum to 1
all_genera <- as.data.frame(t(all_genera)) %>%
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))

# Confirm normalization
all_genera %>% summarise(across(everything(), sum, na.rm = TRUE))

# Transpose back
all_genera <- as.data.frame(t(all_genera))

# Convert row names to a column for merging
all_genera <- all_genera %>% tibble::rownames_to_column(var = "Sample")

#--------------------------------------#
# Load Grouping Metadata
#--------------------------------------#

# Set working directory for metadata
setwd("~/Downloads/comparison_data")

# Read metadata file
group_meta <- read.table("phylum_comp_meta.txt", header = TRUE)

# Set the factor level order of 'meta' to match its appearance in the metadata
group_meta$meta <- factor(group_meta$meta, levels = unique(group_meta$meta))

# Merge abundance data with metadata
all_genera <- left_join(group_meta, all_genera, by = "Sample")

# Convert to long format for ggplot
long_all_genera <- pivot_longer(all_genera, cols = -c(Sample, meta, host),
                                names_to = "Phylum", values_to = "Abundance")

# Compute mean and standard error per group
summary_data <- long_all_genera %>%
  group_by(meta, Phylum) %>%
  summarise(
    Mean_Abundance = mean(Abundance, na.rm = TRUE),
    SE = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# grouping
summary_data <- summary_data %>%
  mutate(Sample_Type = ifelse(grepl("M_", meta, ignore.case = TRUE), "Mucosal", "Luminal"))
#--------------------------------------#
# Generate Stacked Bar Plot with Error Bars
#--------------------------------------#

p <- ggplot(summary_data, aes(x = meta, y = Mean_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack",alpha = 0.8) +
  theme_minimal() +
  labs(title = "Phylum-level Relative Abundances in Luminal and Mucosal Samples", x = "", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_rect(fill = "grey90", color = "grey50"),
        strip.text = element_text(color = "black")) +
  facet_wrap(~ Sample_Type, scales = "free_x")

ggsave("./3A.png", p, width = 11, height = 5, bg = "white")

