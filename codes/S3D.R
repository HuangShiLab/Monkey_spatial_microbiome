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
  filtered_data <- data[, which(colSums(data != 0) / nrow(data) > prev_threshold)]
  return(filtered_data)
}

filter_features_by_abundance <- function(data, max_abd_cutoff=0.001){
  hist(apply(data, 1, max))
  data<-data[apply(data, 1, max) > max_abd_cutoff,]
  return(data)
}

#--------------------------------------#
# Data Import (Mucosa)
#--------------------------------------#

# Set working directory for mucosa data
setwd("~/Downloads/monkey_mucosa")

# Read genus abundance data and metadata
mucosa_species <- as.data.frame(read.table("species_abundance.txt", header = TRUE))

#--------------------------------------#
# Data Import (Content)
#--------------------------------------#

# Set working directory for content data
setwd("~/Downloads/monkey_data")

# Read genus abundance data and metadata
content_species <- as.data.frame(read.table("2b_species_abundance.txt", header = TRUE))


# combine mucosal and luminal
combined_species <- merge(mucosa_species, content_species, 
                          by = "Species", 
                          all = TRUE)
combined_species[is.na(combined_species)] <- 0

# Restore species names as row names and remove the merged column
rownames(combined_species) <- combined_species$Species
combined_species$Species <- NULL

# filter
combined_species <- filter_features_by_abundance(combined_species)
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

# Run ANCOM
ancom_out <- ancom(data = feature_table,
                   taxa_are_rows = TRUE,
                   meta_data = meta_data,
                   main_var = "host",          # or whatever your group variable is
                   p_adj_method = "BH",
                   prv_cut = 0.10,
                   lib_cut = 0,
                   struc_zero = TRUE,
                   neg_lb = TRUE,
                   alpha = 0.01,
                   verbose = TRUE)

# Extract significant results
res_df <- ancom_out$res
significant_taxa <- res_df %>% filter(detected_0.7 == TRUE)

# Subset significant taxa from data
sig_taxa_names <- significant_taxa$taxon_id
sig_data <- combined_species[, sig_taxa_names, drop = FALSE]

# Transpose for plotting/labeling
sig_data_t <- t(sig_data)
colnames(sig_data_t) <- group_meta$meta  # use your metadata labels

# Relabel duplicate column names
new_labels <- colnames(sig_data_t)
seen <- character()
for (i in seq_along(new_labels)) {
  if (new_labels[i] %in% seen) {
    new_labels[i] <- ""
  } else {
    seen <- c(seen, new_labels[i])
  }
}

# Plot heatmap with clustering
p <- pheatmap(
  log10(ifelse(significant_data == 0, 1e-6, significant_data)),  # Log-transformed abundance
  color = viridis(100),
  cluster_rows = TRUE,  # Enable row clustering
  cluster_cols = FALSE, # Disable column clustering
  scale = "none",       # Do not scale rows/columns
  main = "Species Abundance Heatmap",
  show_rownames = FALSE,
  labels_col = new_labels,
  angle_col = 90,
  fontsize_col = 8,
  border_color = NA
)
ggsave("./3D.png", p, width = 10, height = 5, bg = "white")



# Get the order of rows after clustering
reordered_rows <- p$tree_row$order

# Get the row names in clustered order
clustered_row_names <- rownames(significant_data)[reordered_rows]
# Reorder significant_data using clustered_row_names
final_export <- significant_data[clustered_row_names, ]

# Convert to data frame for writing to Excel
final_export_df <- as.data.frame(final_export)
final_export_df$species <- rownames(final_export_df)

# Write to Excel
write_xlsx(final_export_df, path = "./differential_M_L_species.xlsx")

