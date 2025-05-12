#######################################
# Monkey Study
#
# Human monkey comparison --Part II (Pathway abundance)
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
library(phyloseq)
library(palettes)
library(umap)
library(vegan)
library(openxlsx)
library(tibble)
library(ggpubr)
library(patchwork)

class_color <- c("esophagus" = "#eacc76", 
                 "large intestine" = "#5c7272",
                 "small intestine" = "#acab4b",
                 "kidney" = "#7c99bc", 
                 "oral" = "#a5b3c1",
                 "stomach" = "#c7a8a3", 
                 "stool" = "#e9b962")

# Reading the metadata
setwd("~/Downloads/nature2023_data/")

# Read feature table
pathway_human <- read.table("unstratified_merged_pathway.tsv", 
                      sep = '\t', row.names = 1,header = TRUE)
pathway_human[] <- lapply(pathway_human, function(x) as.numeric(as.character(x)))
pathway_human <- pathway_human %>% rownames_to_column("name")

setwd("~/Downloads/monkey_data")
pathway_monkey <- read.table("merged_unstratified_pathway.tsv", 
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

OTU <- otu_table(pathway,taxa_are_rows = TRUE)

# Meta input
setwd("~/Downloads/comparison_data")
metadata <- as.data.frame(read_excel("meta.xlsx"))
rownames(metadata) <- metadata[[1]]
metadata[[1]] <- NULL
META <- sample_data(metadata)

# build phyloseq object
all <- phyloseq(OTU, META)
all <- subset_samples(all, location != "kidney")

bc_distance <- phyloseq::distance(all, method = "bray")
pcoa_results <- ordinate(all, method = "PCoA", distance = bc_distance)  # Perform PCoA

# Make Plot
pcoa_plot <- NULL
pcoa_plot <- plot_ordination(all, pcoa_results, color = 'location', shape = 'species') +
  geom_point(size = 3, alpha = 0.6) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black")) +
  labs(x = "PC1 (35.5%)", y = "PC2 (18.5%)") +
  ggtitle(paste("PCoA of HUMAnN Pathway Abundance (BCD)"))  # Create PCoA plot with title
pcoa_plot

# Extract PC1 values and metadata
pcoa_data <- as.data.frame(pcoa_results$vectors[, 1])  # Extract PC1 values
colnames(pcoa_data) <- "PC1"
pcoa_data$species <- sample_data(all)$species  # Add species information

# statistical analysis
stat_test <- pcoa_data %>%
  wilcox_test(as.formula("PC1 ~ species")) %>%  # Perform Wilcoxon test for each group
  add_xy_position(x = "species") %>% filter (xmin == "1")
#  add_significance("p")

# Create the box plot for PC1 values by species
box_plot1 <- ggplot(pcoa_data, aes(x = PC1, y = species, fill = species)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = c("#eacc76","#a76b3e")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 12),  
        axis.text.y = element_text(size = 10),  
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 5, r = 10, b = 10, l = 10)) +
  labs(x = NULL, y = NULL)

# Combine PCoA plot and box plot using patchwork
combined_plot <- pcoa_plot / box_plot1 + plot_layout(heights = c(7, 2))
combined_plot

ggsave("./1B.png", combined_plot, width = 7, height = 6)

# Extract PC2 values and metadata
pcoa_data <- as.data.frame(pcoa_results$vectors[, 2]) 
colnames(pcoa_data) <- "PC2"
pcoa_data$location <- sample_data(all)$location
pcoa_data$location <- 
  factor(pcoa_data$location, levels = c("oral", "esophagus", "stomach", "small intestine", "large intestine", "stool"))

# statistical analysis
stat_test <- pcoa_data %>%
  wilcox_test(as.formula("PC2 ~ location")) %>%  # Perform Wilcoxon test for each group
  add_xy_position(x = "location") %>% 
  add_significance("p")

# Create the box plot for PC2 values by body sites
box_plot2 <- ggplot(pcoa_data, aes(x = location, y = PC2, fill = location)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 12),  
        axis.text.x = element_text(size = 8, angle = 90),  
        axis.text.y = element_blank(), 
        plot.margin = margin(t = 5, r = 10, b = 10, l = 10)) +
  labs(x = NULL, y = NULL)
box_plot2

# umap
umap_result <- umap(as.matrix(bc_distance))
umap_result <- as.data.frame(umap_result$layout)
colnames(umap_result) <- c("UMAP1", "UMAP2")

# Add metadata information
metadata <- all@sam_data
umap_result$group <- metadata$species
umap_result$group2 <- metadata$location
umap_result$group3 <- metadata$age

# Plot umap
umap_plot <- ggplot(umap_result) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = group2, shape = group),
  size = 3, alpha = 0.8) +
  xlab("UMAP1") + ylab("UMAP2") +
  ggtitle("UMAP (BCD)") +
  theme_minimal()
umap_plot

