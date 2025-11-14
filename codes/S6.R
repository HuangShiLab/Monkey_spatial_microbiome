#######################################
# Monkey project
#
# Mucosal & Luminal -- Heatmap
#
# Author: HOU Shuwen
#######################################

## install and load necessary libraries for data analyses
#-------------------------------
p <- c("ggplot2","pheatmap", "vegan", "ANCOMBC",
       "ggpubr", "dplyr", "writexl","tidyr","tibble", "viridis")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
#-------------------------------

#function to filter
filter_features_by_prev <- function(data, prev_threshold = 0.1) {
  # Filter columns based on the proportion of non-zero entries
  filtered_data <- data[which(colSums(data != 0) / ncol(data) > prev_threshold),]
  return(filtered_data)
}

filter_features_by_abundance <- function(data, mean_abd_cutoff=0.001){
  hist(apply(data, 1, mean))
  data<-data[apply(data, 1, mean) > mean_abd_cutoff,]
  return(data)
}

#--------------------------------------#
# Data Import (Mucosa)
#--------------------------------------#

# Set working directory for mucosa data
setwd("~/Downloads/monkey_mucosa")

# Read genus abundance data and metadata
mucosa_species <- as.data.frame(read.table("species_abundance_intestine.txt", header = TRUE))
mucosa_species <- mucosa_species %>% column_to_rownames("Species")

mucosa_species <- filter_features_by_abundance(mucosa_species)

mucosa_species <- mucosa_species %>%
  rownames_to_column("species")
#--------------------------------------#
# Data Import (Content)
#--------------------------------------#

# Set working directory for content data
setwd("~/Downloads/monkey_data")

# Read genus abundance data and metadata
content_species <- as.data.frame(read.table("2b_species_abundance_intestine.txt", header = TRUE))

content_species <- content_species %>% column_to_rownames("Species")

content_species <- filter_features_by_abundance(content_species)

content_species <- content_species %>%
  rownames_to_column("species")

# combine mucosal and luminal
combined_species <- merge(content_species, mucosa_species,
                          by = "species", 
                          all = TRUE)
combined_species[is.na(combined_species)] <- 0

# Restore species names as row names and remove the merged column
rownames(combined_species) <- combined_species$species
combined_species$species <- NULL

# filter
# combined_species <- filter_features_by_abundance(combined_species)
combined_species <- sweep(combined_species, 2, colSums(combined_species), FUN = "/")


# Set working directory for metadata
setwd("~/Downloads/comparison_data")

# Read metadata file
group_meta <- read.table("phylum_comp_meta.txt", header = TRUE)
ID <- group_meta$Sample  # Extract sample IDs
combined_species <- combined_species[, ID, drop = FALSE]  # Keep only relevant columns

combined_species <- as.data.frame(t(combined_species))
rowSums(combined_species, na.rm = TRUE) # check row sums

# Ensure metadata and species table align
combined_species <- combined_species[order(match(rownames(combined_species), group_meta$Sample)), ]
group_meta <- group_meta[match(rownames(combined_species), group_meta$Sample), ]

# Prepare inputs for ANCOM
# Taxa must be in rows, samples in columns → transpose back
feature_table <- t(combined_species)
meta_data <- group_meta
rownames(meta_data) <- meta_data$Sample  # Match by Sample ID


# pseudocounts
counts <- round(feature_table * 1e6)

# Run ANCOM
ancom_out <- ancombc2(data = counts,
                   taxa_are_rows = TRUE,
                   meta_data = meta_data,
                   fix_formula = "host")

# Extract significant results
res_df <- ancom_out$res
significant_taxa <- res_df %>% filter(diff_hostmucosa == TRUE)
significant_taxa <- 
  significant_taxa[order(significant_taxa$lfc_hostmucosa, decreasing = TRUE), ]

# Subset significant taxa from data
sig_taxa_names <- significant_taxa$taxon
sig_data <- combined_species[, sig_taxa_names, drop = FALSE]
sig_data <- as.matrix(sig_data)


# Plot heatmap with clustering
p <- pheatmap(
  log(ifelse(sig_data == 0, 1e-6, sig_data)),  # Log-transformed abundance
  color = viridis(100),
  cluster_rows = FALSE,  # Enable row clustering
  cluster_cols = FALSE, # Disable column clustering
  scale = "none",       # Do not scale rows/columns
  main = "Species Abundance Heatmap",
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  angle_col = 315,
  border_color = NA
)

ggsave("./3D.png", p, width = 4.5, height = 6, bg = "white")


# Convert to data frame for writing to Excel
final_export_df <- as.data.frame(t(sig_data))
final_export_df$species <- rownames(final_export_df)


# Write to Excel
write_xlsx(final_export_df, path = "./differential_M_L_species_ANCOM.xlsx")

