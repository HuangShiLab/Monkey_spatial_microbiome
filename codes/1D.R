#######################################
# Monkey Study
#
# Human monkey comparison --Part IV (copri Pathway)
#
#
# Author: HOU Shuwen
#######################################

# loading packages
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(vegan)
library(openxlsx)
library(tibble)
library(ggpubr)
library(patchwork)
library(clusterProfiler)
library(stringr)

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

# Initialize a results data frame
results <- data.frame(Feature = character(), p.value = numeric(), stringsAsFactors = FALSE)

# Loop through each feature (excluding SampleID / group row)
feature_cols <- setdiff(rownames(gene_filtered), "1")
for (feature in feature_cols) {
  kruskal_test <- kruskal.test(as.numeric(gene_filtered[feature, ]) ~ as.factor(as.character(gene_filtered["1", ])))
  results <- rbind(results, data.frame(Feature = feature, p.value = kruskal_test$p.value))
}

# Adjust p-values using BH method
results$adj.p.value <- p.adjust(results$p.value, method = "BH")

# Compute fold change for all features
fold_changes <- sapply(results$Feature, function(feature) {
  group_labels <- as.factor(as.character(gene_filtered["1", ]))
  values <- as.numeric(gene_filtered[feature, ])
  group_medians <- tapply(values, group_labels, median)
  max(group_medians) / max(min(group_medians), 1e-6)  # safe division
})

# Add to results
results$fold_change <- fold_changes

# Filter significant features
differential_features <- results[results$fold_change > 2, ]

# Add enrichment column
enrich <- sapply(differential_features$Feature, function(feature) {
  group_labels <- as.factor(as.character(gene_filtered["1", ]))
  values <- as.numeric(gene_filtered[feature, ])
  group_medians <- tapply(values, group_labels, median)
  names(which.max(group_medians))  # return the group with highest median
})

# Add to the results
differential_features$enrich <- enrich
write.xlsx(differential_features, "enriched_genes.xlsx")

# import GO data
universe_genes <- unique(rownames(gene_filtered))
writeLines(universe_genes, "allgenes.txt")
# copy all genes to UniPort website to get GO information
setwd("~/Downloads/")
GO_all <- read_tsv("all_GOs.tsv", quote = "")

# Use ClusterProfiler for GO enrichment analysis
gene2go <- GO_all %>%
  dplyr::select(GOID = `Gene Ontology (GO)`, GeneID = Feature) %>%
  filter(!is.na(GOID) & GOID != "") %>%
  separate_rows(GOID, sep = ";") %>%   # split rows
  mutate(GOID = trimws(GOID))

# Create list of significant genes
sig_genes <- differential_features %>%
  filter(enrich == "monkey") %>%
  pull(Feature) %>%
  unique()

# Perform enrichment
ego <- enricher(sig_genes,
                TERM2GENE = gene2go,
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.5,
                minGSSize = 2)

ego@result <- ego@result[ego@result$pvalue < 0.1, ]
ego@result$GeneRatio <- sapply(ego@result$GeneRatio, function(x) eval(parse(text = x))) 
ego@result$Description_wrapped <- str_wrap(ego@result$Description, width = 40)

dotplot <- ggplot(ego@result, aes(
  x = FoldEnrichment,
  y = reorder(Description_wrapped, FoldEnrichment),
  size = Count,
  color = pvalue
)) +
  geom_point() +
  scale_color_continuous(low = "red", high = "blue") +
  scale_size(range = c(3, 10)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.y = element_text(size = 12),   # Y-axis text size
    axis.text.x = element_text(size = 10),   # X-axis text size (optional)
    plot.title = element_text(size = 14, face = "bold")  # Title size
  ) +
  labs(
    title = "Monkey S. copri GO Enrichment",
    x = "Gene Ratio",
    y = "GO Term"
  )
ggsave("./1D.png", dotplot, width = 8, height = 4)
