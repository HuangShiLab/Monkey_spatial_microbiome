#######################################
# Monkey Study
#
# Human monkey comparison --Part IV (copri Pathway)
#
#
# Author: HOU Shuwen
#######################################

# loading packages
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(cutpointr)
library(purrr)
library(reshape2)
library(palettes)
library(umap)
library(vegan)
library(openxlsx)
library(tibble)
library(ggpubr)
library(patchwork)
library(readxl)

# Reading the metadata
setwd("~/Downloads/nature2023_data/")

# Read feature table
pathway_human <- as.data.frame(read_excel("pathway_abundance_copri.xlsx"))
rownames(pathway_human) <- pathway_human[[1]]
pathway_human <- pathway_human[, -1]
pathway_human[] <- lapply(pathway_human, function(x) as.numeric(as.character(x)))
pathway_human <- pathway_human %>% rownames_to_column("name")

setwd("~/Downloads/monkey_data")
pathway_monkey <- read.table("copri_pathway_abundance.tsv", 
                             sep = '\t', row.names = 1,header = TRUE)
pathway_monkey[] <- lapply(pathway_monkey, function(x) as.numeric(as.character(x)))
pathway_monkey <- pathway_monkey %>% rownames_to_column("name")

# Combine two feature tables
pathway <- full_join(pathway_human, pathway_monkey, by = "name")
pathway <- as.data.frame(pathway)
pathway[is.na(pathway)] <- 0

rownames(pathway) <- pathway$name
pathway <- pathway[,-1]
pathway[] <- lapply(pathway, function(x) as.numeric(as.character(x)))
pathway <- t(pathway)

# Meta input
setwd("~/Downloads/comparison_data")
metadata <- as.data.frame(read_excel("meta.xlsx"))
rownames(metadata) <- metadata[[1]]
metadata[[1]] <- NULL

pathway <- pathway[order(match(rownames(pathway), rownames(metadata))), ]
pathway <- t(pathway)
pathway <- pathway[-1, ] # delete the first row
pathway <- as.matrix(pathway)

new_labels <- metadata$location
last_label <- ""
for (i in seq_along(new_labels)) {
  if (new_labels[i] == last_label) {
    new_labels[i] <- ""  # Blank out duplicates
  } else {
    last_label <- new_labels[i]
  }
}



# Plot heatmap with clustering
p <- pheatmap(
  log10(ifelse(pathway == 0, 0.01, pathway)),  # Log-transformed abundance
  color = viridis(100),
  cluster_rows = TRUE,  # Enable row clustering
  cluster_cols = FALSE, # Disable column clustering
  scale = "none",       # Do not scale rows/columns
  main = "Pathway Abundance Heatmap",
  show_rownames = FALSE,
  labels_col = new_labels,
  angle_col = 90,
  fontsize_col = 8,
  border_color = NA
)

# copri pathway comparison
pathway_filtered <- pathway[, apply(pathway, 2, median) > 0]
pathway_filtered <- pathway_filtered[, !(names(pathway_filtered) == "O1_K_1")]
pathway_filtered <- sweep(pathway_filtered, 2, colSums(pathway_filtered), FUN = "/")

# Create a new row based on column names
new_row <- ifelse(startsWith(colnames(pathway_filtered), "SRR"), "human", "monkey")

# Add the new row to the dataframe
pathway_filtered <- as.data.frame(rbind(new_row, pathway_filtered))

# Initialize a results list
results <- data.frame(Feature = character(), p.value = numeric(), stringsAsFactors = FALSE)
# Loop through each feature (excluding SampleID and Group)
feature_cols <- setdiff(rownames(pathway_filtered), c("new_row"))
for (feature in feature_cols) {
  kruskal_test <- kruskal.test(as.numeric(pathway_filtered[feature, ]) ~ 
                                 as.factor(as.character(pathway_filtered[1, ])))
  results <- rbind(results, data.frame(Feature = feature, p.value = kruskal_test$p.value))
}
# Adjust p-values using the Benjamini-Hochberg method
results$adj.p.value <- p.adjust(results$p.value, method = "BH")
results <- results[!is.na(results$adj.p.value), ]

# Filter significant features
significant_features <- results[results$adj.p.value < 0.05, ]

# Add enrichment column
enrich <- sapply(significant_features$Feature, function(feature) {
  group_labels <- as.factor(as.character(pathway_filtered["1", ]))
  values <- as.numeric(pathway_filtered[feature, ])
  group_medians <- tapply(values, group_labels, median)
  names(which.max(group_medians))  # return the group with highest median
})

# Add to the results
significant_features$enrich <- enrich
significant_data <- pathway_filtered[significant_features$Feature, ]
significant_data_clean <- as.data.frame(
  lapply(significant_data, function(x) as.numeric(as.character(x))))
rownames(significant_data_clean) <- rownames(significant_data)

# Now convert to matrix for pheatmap
significant_data_clean <- as.matrix(significant_data_clean)

p <- pheatmap(
  significant_data_clean,
  color = viridis(100),
  cluster_rows = TRUE,  # Enable row clustering
  cluster_cols = FALSE, # Disable column clustering
  scale = "none",       # Do not scale rows/columns
  main = "Species Abundance Heatmap",
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 90,
  fontsize_col = 8,
  border_color = NA
)

## gene
# Reading the metadata
setwd("~/Downloads/nature2023_data/")

# Read feature table
gene_human <- read.table("copri_gene_abundance.tsv", 
                            sep = '\t', row.names = 1,header = TRUE)
gene_human[] <- lapply(gene_human, function(x) as.numeric(as.character(x)))
gene_human <- gene_human %>% rownames_to_column("name")

setwd("~/Downloads/monkey_data")
gene_monkey <- read.table("copri_gene_abundance.tsv", 
                             sep = '\t', row.names = 1,header = TRUE)
gene_monkey[] <- lapply(gene_monkey, function(x) as.numeric(as.character(x)))
gene_monkey <- gene_monkey %>% rownames_to_column("name")

# Combine two feature tables
gene <- full_join(gene_human, gene_monkey, by = "name")
gene <- as.data.frame(gene)
gene[is.na(gene)] <- 0

rownames(gene) <- gene$name
gene <- gene[,-1]
gene[] <- lapply(gene, function(x) as.numeric(as.character(x)))

gene_filtered <- gene[, apply(gene, 2, max) > 0]
gene_filtered <- sweep(gene_filtered, 2, colSums(gene_filtered), FUN = "/")

# Create a new row based on column names
new_row <- ifelse(startsWith(colnames(gene_filtered), "SRR"), "human", "monkey")

# Add the new row to the dataframe
gene_filtered <- rbind(new_row, gene_filtered)

# Initialize a results list
results <- data.frame(Feature = character(), p.value = numeric(), stringsAsFactors = FALSE)
# Loop through each feature (excluding SampleID and Group)
feature_cols <- setdiff(rownames(gene_filtered), c("1"))
for (feature in feature_cols) {
  kruskal_test <- kruskal.test(as.numeric(gene_filtered[feature, ]) ~ as.factor(as.character(gene_filtered[1, ])))
  results <- rbind(results, data.frame(Feature = feature, p.value = kruskal_test$p.value))
}
# Adjust p-values using the Benjamini-Hochberg method
results$adj.p.value <- p.adjust(results$p.value, method = "BH")
results <- results[!is.na(results$adj.p.value), ]

# Filter significant features
significant_features <- results[results$adj.p.value < 0.01, ]

# Add enrichment column
enrich <- sapply(significant_features$Feature, function(feature) {
  group_labels <- as.factor(as.character(gene_filtered["1", ]))
  values <- as.numeric(gene_filtered[feature, ])
  group_medians <- tapply(values, group_labels, median)
  names(which.max(group_medians))  # return the group with highest median
})

# Add to the results
significant_features$enrich <- enrich
write.xlsx(significant_features, "enriched_genes.xlsx")

# calculate non zero prevalance in each group
# Prepare group assignment
group1 <- colnames(gene_filtered)[new_row == "human"]
group2 <- colnames(gene_filtered)[new_row == "monkey"]

# Initialize data frame to store prevalence
feature_rows <- setdiff(rownames(gene_filtered), "1")

prevalence_df <- data.frame(
  Feature = feature_rows,
  human_nonzero = NA_real_,
  monkey_nonzero = NA_real_,
  stringsAsFactors = FALSE
)

# Calculate non-zero prevalence
for (i in seq_along(feature_rows)) {
  feature <- feature_rows[i]
  values_human <- as.numeric(gene_filtered[feature, group1])
  values_monkey <- as.numeric(gene_filtered[feature, group2])
  
  prevalence_df$human_nonzero[i]  <- mean(values_human != 0)
  prevalence_df$monkey_nonzero[i] <- mean(values_monkey != 0)
}

filtered <- prevalence_df[
  (prevalence_df$human_nonzero >= 0.8 & prevalence_df$monkey_nonzero <= 0.2) |
    (prevalence_df$monkey_nonzero >= 0.8 & prevalence_df$human_nonzero <= 0.2),
]

