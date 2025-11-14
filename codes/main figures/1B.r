#######################################
# Monkey Study
#
# 5 most common phylum, stacked bar plot, monkey and human
#
#
# Author: HOU Shuwen
#######################################

# loading packages
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(phyloseq)
library(palettes)
library(vegan)
library(openxlsx)
library(tibble)
library(ggpubr)
library(patchwork)

#--------------------------------------#
# Data Import (Monkey)
#--------------------------------------#

# Set working directory for content data
setwd("~/Downloads/monkey_data")

# Read monkey phylum abundance data and metadata
monkey_phylum <- as.data.frame(read_tsv("phylum_abundance_bracken.tsv"))
metadata <- read.table("kraken_metadata.txt", sep="\t", header=TRUE)
ID <- metadata$ID  # Extract sample IDs

# Format data: set row names and convert to numeric
rownames(monkey_phylum) <- monkey_phylum$name
monkey_phylum <- monkey_phylum[, -1]
monkey_phylum[] <- lapply(monkey_phylum, function(x) as.numeric(as.character(x)))

setwd("~/Downloads/nature2023_data")

# Read human phylum abundance data and metadata
human_phylum <- as.data.frame(read_tsv("phylum_abundance_bracken.tsv"))
metadata <- read_excel("meta_WMS.xlsx")
metadata <- metadata %>% filter(subject_id %in% c(6,8,12,14))
ID <- metadata$Run  # Extract sample IDs
ID <- intersect(ID, colnames(human_phylum))

# Format data: set row names and convert to numeric
rownames(human_phylum) <- human_phylum$name
human_phylum <- human_phylum[, -1]
human_phylum[] <- lapply(human_phylum, function(x) as.numeric(as.character(x)))
human_phylum <- human_phylum[, ID, drop = FALSE]


#--------------------------------------#
# Select Top 5 Phylum for Each Sample
#--------------------------------------#

get_top_n_rows <- function(df, n_top = 5) {
  # Ensure row names are retained by creating a phylum column
  df$phylum <- rownames(df)
  
  # Extract the top n genera for each column
  top_rows <- unique(do.call(rbind, lapply(names(df)[-ncol(df)], function(col) {
    # Sort by descending order of the column and select top n rows
    ordered_df <- df[order(-df[[col]]), ]
    ordered_df[1:n_top, ]
  })))
  
  # Restore row names and remove the temporary phylum column
  rownames(top_rows) <- top_rows$phylum
  top_rows$phylum <- NULL
  
  # Calculate the sum of the remaining genera for each column
  # Exclude the phylum column when summing
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

top_5_human_phylum <- get_top_n_rows(human_phylum) %>%
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
top_5_human_phylum %>% summarise(across(everything(), sum, na.rm = TRUE))

# Process content data
# Normalize values so each row sums to 1
# Check sum to confirm normalization

top_5_monkey_phylum <- get_top_n_rows(monkey_phylum) %>%
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
top_5_monkey_phylum %>% summarise(across(everything(), sum, na.rm = TRUE))

#--------------------------------------#
# Merge Data from Both Sources
#--------------------------------------#

# Convert row names to a column for merging
top_5_human_phylum <- top_5_human_phylum %>% tibble::rownames_to_column("name")
top_5_monkey_phylum <- top_5_monkey_phylum %>% tibble::rownames_to_column("name")

# Merge both datasets and replace NA with 0
all_genera <- full_join(top_5_human_phylum, top_5_monkey_phylum, by = "name")
all_genera[is.na(all_genera)] <- 0

# Restore row names and remove extra column
rownames(all_genera) <- all_genera[, 1]
all_genera <- all_genera[, -1]


#--------------------------------------#
# Data Cleaning and Normalization
#--------------------------------------#

# Remove low-abundance genera (max < 0.05)
all_genera <- as.data.frame(t(all_genera))
columns_to_remove <- all_genera %>%
  summarise(across(everything(), max, na.rm = TRUE)) %>%
  select(where(~ . < 0.1)) %>%
  colnames()
all_genera <- all_genera %>% select(-all_of(columns_to_remove))

others_new <- 1 - rowSums(all_genera, na.rm = TRUE)
all_genera <- as.data.frame(t(all_genera))

all_genera["others",] <- all_genera["others",] + others_new

all_genera <- as.data.frame(t(all_genera))

# Convert row names to a column for merging
all_genera <- all_genera %>% tibble::rownames_to_column(var = "Run")

#--------------------------------------#
# Load Grouping Metadata
#--------------------------------------#

# Set working directory for metadata
setwd("~/Downloads/comparison_data")

# Read metadata file
group_meta <- read_excel("meta.xlsx")
group_meta <- group_meta %>% 
  mutate(group = paste(species, location, sep = "_")) %>%
  select(Run, group)

# Set the factor level order of 'meta' to match its appearance in the metadata
group_meta$group <- factor(group_meta$group, levels = unique(group_meta$group))

# Convert to long format for ggplot
long_all_genera <- pivot_longer(all_genera, cols = c(-Run),
                                names_to = "Phylum", values_to = "Abundance")

# Merge abundance data with metadata
long_all_genera <- left_join(long_all_genera, group_meta, by = "Run")

# Compute mean and standard error per group
summary_data <- long_all_genera %>%
  mutate(group = recode(
      group,
      "monkey_oral" = "monkey_upper",
      "monkey_esophagus" = "monkey_upper"
    )) %>%
  group_by(group, Phylum) %>%
  summarise(
    Abundance = mean(Abundance, na.rm = TRUE),
    SE = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# grouping
summary_data <- summary_data %>% filter(group != "monkey_kidney") %>%
  mutate(Sample_Type = ifelse(grepl("human_", group, ignore.case = TRUE), "Human", "Monkey"))
#--------------------------------------#
# Generate Stacked Bar Plot with Error Bars
#--------------------------------------#

p <- ggplot(summary_data, aes(x = group, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack",alpha = 0.8) +
  theme_minimal() +
  labs(title = "Human vs Monkey Phylum-level Relative Abundances", x = "", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey90", color = "grey50"),
        strip.text = element_text(color = "black")) +
  facet_wrap(~ Sample_Type, scales = "free_x")
p

ggsave("./1F.png", p, width = 5.5, height = 4, bg = "white")



