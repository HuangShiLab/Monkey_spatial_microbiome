#######################################
# Monkey project
#
# Humann result -- Heatmap
#
# Author: HOU Shuwen
#######################################

## install and load necessary libraries for data analyses
#-------------------------------
p <- c("ggplot2","pheatmap", "vegan", "readr","file2meco","RColorBrewer",
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

filter_features_by_abundance <- function(data, mean_abd_cutoff=0.001){
  hist(colMeans(data))
  data<-data[, which(colMeans(data) > mean_abd_cutoff)]
  data
}

# Species
# data input
setwd("~/Downloads/monkey_data")
bc_kraken <- as.data.frame(read_tsv("species_abundance_bracken.tsv",col_names = T))
bc_kraken$name <- gsub("'", "", bc_kraken$name)

# tidy up feature table
rownames(bc_kraken) <- bc_kraken$name
bc_kraken <- bc_kraken[,-1]
bc_kraken[] <- lapply(bc_kraken, function(x) as.numeric(as.character(x)))
# Normalize each column to sum to 1
bc_kraken <- sweep(bc_kraken, 2, colSums(bc_kraken), FUN = "/")
bc_kraken <- as.data.frame(t(bc_kraken))

rowSums(bc_kraken, na.rm = TRUE) # check row sums

metadata <- read.table("kraken_metadata.txt",sep="\t",header=TRUE)


bc_kraken <- bc_kraken[order(match(rownames(bc_kraken), rownames(metadata))), ]
# bc_kraken <- filter_features_by_prev(bc_kraken)

# Select significant features
bc_kraken$group <- metadata$group

# Remove kidney samples
bc_kraken <- bc_kraken[bc_kraken$group != "kidney", ]


# Initialize a results list
results <- data.frame(Feature = character(), p.value = numeric(), stringsAsFactors = FALSE)
# Loop through each feature (excluding SampleID and Group)
feature_cols <- setdiff(names(bc_kraken), c("SampleID", "group"))
for (feature in feature_cols) {
  kruskal_test <- kruskal.test(bc_kraken[[feature]] ~ bc_kraken$group)
  results <- rbind(results, data.frame(Feature = feature, p.value = kruskal_test$p.value))
}
# Adjust p-values using the Benjamini-Hochberg method
results$adj.p.value <- p.adjust(results$p.value, method = "BH")
results <- results[!is.na(results$adj.p.value), ]

# Filter significant features
significant_features <- results[results$adj.p.value < 0.01, ]
significant_data <- bc_kraken %>%
  select(significant_features$Feature)
significant_data <- as.matrix(t(significant_data))

# Plot heatmap with clustering
p <- pheatmap(
  log10(ifelse(significant_data == 0, 1e-6, significant_data)),  # Log-transformed abundance
  color = viridis(100),
  cluster_rows = TRUE,  # Enable row clustering
  cluster_cols = FALSE, # Disable column clustering
  scale = "none",       # Do not scale rows/columns
  main = "Species Abundance Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA
)
ggsave("./plots/S4.png", p, width = 10, height = 7)

# Get the order of rows after clustering
reordered_rows <- p$tree_row$order

# Get the row names in clustered order
clustered_row_names <- rownames(significant_data)[reordered_rows]
# Reorder significant_data using clustered_row_names
final_export <- significant_data[clustered_row_names, ]

# Convert to data frame for writing to Excel
final_export_df <- as.data.frame(final_export)
final_export_df <- tibble::rownames_to_column(final_export_df, var = "Sample")


# Pathway
# Data input
setwd("~/Downloads/monkey_data")
pathway <- read.table("merged_unstratified_pathway.tsv", 
                      sep = '\t', row.names = 1,header = TRUE)
pathway[] <- lapply(pathway, function(x) as.numeric(as.character(x)))

# Normalize each column to sum to 1
pathway <- sweep(pathway, 2, colSums(pathway), FUN = "/")
pathway <- as.data.frame(t(pathway))
rowSums(pathway, na.rm = TRUE) # check row sums

metadata <- read.table("kraken_metadata.txt",sep="\t",header=TRUE, row.names = 1)


pathway <- pathway[order(match(rownames(pathway), rownames(metadata))), ]
# pathway <- filter_features_by_prev(pathway)

# Select significant features
pathway$group <- metadata$group


# Remove kidney samples
pathway <- pathway[pathway$group != "kidney", ]
pathway <- pathway[rownames(pathway) != "Y1_SI_113cm", ]

# Initialize a results list
results <- data.frame(Feature = character(), p.value = numeric(), stringsAsFactors = FALSE)
# Loop through each feature (excluding SampleID and Group)
feature_cols <- setdiff(names(pathway), c("SampleID", "group"))
for (feature in feature_cols) {
  kruskal_test <- kruskal.test(pathway[[feature]] ~ pathway$group)
  results <- rbind(results, data.frame(Feature = feature, p.value = kruskal_test$p.value))
}
# Adjust p-values using the Benjamini-Hochberg method
results$adj.p.value <- p.adjust(results$p.value, method = "bonferroni")

# Filter significant features
significant_features <- results[results$adj.p.value < 0.01, ]
significant_data <- pathway %>%
  select(significant_features$Feature) 
significant_data <- as.matrix(t(significant_data))

# add pathway type
annotations <- read_excel("table 2.xlsx")

# Create annotation data frame
annotation_row <- data.frame(
  Type = annotations$Type,
  row.names = annotations$Pathway
)

# Match row order with heatmap data (critical!)
annotation_row <- annotation_row[rownames(significant_data), , drop = FALSE]


# Step 1: Sort the unique types
unique_types <- sort(unique(annotation_row$Type))

# Step 2: Convert Type to a factor with alphabetically ordered levels
annotation_row$Type <- factor(annotation_row$Type, levels = unique_types)

# Step 3: Create smooth plasma gradient colors in the same order
type_colors <- setNames(viridis(length(unique_types), option = "plasma"), unique_types)

# Step 4: Assign to annotation_colors
annotation_colors <- list(Type = type_colors)

# Plot heatmap with clustering
p <- pheatmap(
  log10(ifelse(significant_data == 0, 1e-6, significant_data)),  # Log-transformed abundance
  color = viridis(100),
  cluster_rows = TRUE,  # Enable row clustering
  cluster_cols = FALSE, # Disable column clustering
  scale = "none",       # Do not scale rows/columns
  main = "Pathway Abundance Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  annotation_row = annotation_row,
  annotation_names_row = FALSE,
  annotation_legend = FALSE,
  annotation_colors = list(Type = type_colors)
)
ggsave("./plots/2C.png", p, width = 10, height = 7)

# Get the order of rows after clustering
reordered_rows <- p$tree_row$order

# Get the row names in clustered order
clustered_row_names <- rownames(significant_data)[reordered_rows]
# Reorder significant_data using clustered_row_names
final_export <- significant_data[clustered_row_names, ]

# Convert to data frame for writing to Excel
final_export_df <- as.data.frame(final_export)
final_export_df <- tibble::rownames_to_column(final_export_df, var = "Sample")

colnames(df_pathways)[1] <- "pathway"

# Join with pathway map to get level 1 and level 2 annotations
df_annotated <- df_pathways %>%
  left_join(MetaCyc_pathway_map, by = c("pathway" = "pathway_id"))

# Write to Excel
write_xlsx(final_export_df, path = "./significant_pathways_clustered.xlsx")
