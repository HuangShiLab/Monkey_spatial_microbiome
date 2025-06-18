#######################################
# Monkey Study
#
# Human monkey comparison --  Part I (Taxa abundance)
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
library(vegan)
library(openxlsx)
library(patchwork)
library(ggpubr)

# Reading the metadata
setwd("~/Downloads/nature2023_data/")
md <- read_excel("meta_WMS.xlsx")

# Read feature table
feature_table <- read_tsv("species_abundance_bracken.tsv",col_names = T)

# Clean the columns
md_filter <- md %>% filter(subject_id %in% c(6,8,12,14))
SRR <- md_filter$Run

# Filter the feature table
filtered_otu <- feature_table %>% select(name,all_of(SRR))
filtered_otu$name <- gsub("'", "", filtered_otu$name)

# Combine two feature tables
setwd("~/Downloads/monkey_data")
bc_kraken <- as.data.frame(read_tsv("species_abundance_bracken.tsv",col_names = T))
bc_kraken$name <- gsub("'", "", bc_kraken$name)
bc_kraken <- full_join(filtered_otu, bc_kraken, by = "name")
bc_kraken <- as.data.frame(bc_kraken)
bc_kraken[is.na(bc_kraken)] <- 0


rownames(bc_kraken) <- bc_kraken$name
bc_kraken <- bc_kraken[,-1]
bc_kraken[] <- lapply(bc_kraken, function(x) as.numeric(as.character(x)))
bc_kraken <- bc_kraken %>%
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
bc_kraken %>% summarise(across(everything(), sum, na.rm = TRUE))

OTU <- otu_table(bc_kraken,taxa_are_rows = TRUE)

# Data input
setwd("~/Downloads/comparison_data")
metadata <- as.data.frame(read_excel("meta.xlsx"))
rownames(metadata) <- metadata[[1]]
metadata[[1]] <- NULL
META <- sample_data(metadata)

# build phyloseq object
all <- phyloseq(OTU, META)
all <- subset_samples(all, location != "kidney")
sample_data(all)$location <- factor(sample_data(all)$location, 
                                    levels = c("oral", "esophagus", "stomach", "small intestine", "large intestine"))

bc_distance <- phyloseq::distance(all, method = "bray")
pcoa_results <- ordinate(all, method = "PCoA", distance = bc_distance)  # Perform PCoA
# Make Plot
pcoa_plot <- NULL
pcoa_plot <- plot_ordination(all, pcoa_results, color = 'location', shape = 'species') +
  geom_point(size = 5, alpha = 0.7) + 
  scale_color_viridis_d(option = "plasma", begin = 0, end = 1, direction = 1) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "none"
  ) +
  labs(x = "PC1 (22.6%)", y = "PC2 (17.1%)") +
  ggtitle("PCoA of Kraken2 Species Abundance (BCD)")

# Extract PC1 values and metadata
pcoa_data <- as.data.frame(pcoa_results$vectors[, 1])  # Extract PC1 values
colnames(pcoa_data) <- "PC1"
pcoa_data$species <- sample_data(all)$species  # Add species information

# statistical analysis
#stat_test <- pcoa_data %>%
#  wilcox_test(as.formula("PC1 ~ species")) %>%  # Perform Wilcoxon test for each group
#  add_xy_position(x = "species") %>% filter (xmin == "1")

# Create the box plot for PC1 values by species
box_plot1 <- ggplot(pcoa_data, aes(x = PC1, y = species, fill = species)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = c("#eacc76","#a76b3e")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 16),
        axis.text.x = element_blank(), 
        plot.margin = margin(t = 5, r = 10, b = 10, l = 10)) +
  labs(x = NULL, y = NULL)

# Combine PCoA plot and box plot using patchwork
combined_plot <- pcoa_plot / box_plot1 + plot_layout(heights = c(4, 1))
combined_plot

ggsave("./1A.png", combined_plot, width = 6, height = 6)

# Extract PC2 values and metadata
pcoa_data <- as.data.frame(pcoa_results$vectors[, 2]) 
colnames(pcoa_data) <- "PC2"
pcoa_data$location <- sample_data(all)$location
pcoa_data$location <- 
  factor(pcoa_data$location, levels = c("oral", "esophagus", "stomach", "small intestine", "large intestine", "stool"))

# statistical analysis
#stat_test <- pcoa_data %>%
#  wilcox_test(as.formula("PC2 ~ location")) %>%  # Perform Wilcoxon test for each group
#  add_xy_position(x = "location") %>% 
#  add_significance("p")

# Create the box plot for PC2 values by body sites
box_plot2 <- ggplot(pcoa_data, aes(x = location, y = PC2, fill = location)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  scale_fill_viridis_d(option = "plasma", begin = 0, end = 1, direction = 1) + 
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 90),  
        axis.text.y = element_blank(), 
        plot.margin = margin(t = 5, r = 10, b = 10, l = 10)) +
  labs(x = NULL, y = NULL)
box_plot2
ggsave("./1AS.png", box_plot2, width = 1.5, height = 4)


# umap
#umap_result <- umap(as.matrix(bc_distance))
#umap_result <- as.data.frame(umap_result$layout)
#colnames(umap_result) <- c("UMAP1", "UMAP2")

# Add metadata information
#metadata <- all@sam_data
#umap_result$group <- metadata$species
#umap_result$group2 <- metadata$location
#umap_result$group3 <- metadata$age

# Plot umap
#umap_plot <- ggplot(umap_result) +
#  geom_point(aes(x = UMAP1, y = UMAP2, color = group2, shape = group),
#             size = 3, alpha = 0.8) +
#  scale_shape_manual(values = c(1, 2)) +
#  xlab("UMAP1") + ylab("UMAP2") +
#  ggtitle("UMAP (BCD)") +
#  theme_minimal()
#umap_plot

