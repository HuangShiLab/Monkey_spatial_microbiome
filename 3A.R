#######################################
# Monkey Study
#
# Monkey mucosal PCoA (blood)
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
library(patchwork)
library(metagMisc)

# data input
setwd("~/Downloads/monkey_mucosa")
feature_table <- as.data.frame(read.table("species_abundance.txt", header = TRUE))

metadata <- read.table("meta.txt",sep="\t",header=TRUE)
ID <- metadata$ID

# tidy up feature table
rownames(feature_table) <- feature_table$Species
feature_table <- feature_table[,-1]
feature_table[] <- lapply(feature_table, function(x) as.numeric(as.character(x)))
feature_table <- feature_table[, ID, drop = FALSE]

# make sure sum abundance equals to 1
feature_table <- feature_table %>%
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))

OTU <- otu_table(feature_table,taxa_are_rows = TRUE)

# metadata
rownames(metadata) <- metadata[[1]]
metadata[[1]] <- NULL
META <- sample_data(metadata)

# build phyloseq object
all <- phyloseq(OTU, META)

# calculate PCoA results
bc_distance <- phyloseq::distance(all, method = "bray")
pcoa_results <- ordinate(all, method = "PCoA", distance = bc_distance)  # Perform PCoA

# Make Plot
pcoa_plot1 <- NULL
pcoa_plot1 <- plot_ordination(all, pcoa_results, color = 'location',shape = 'age') +
  geom_point(size = 3, alpha = 0.6) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black"))+
  labs(x = "PC1 (27%)", y = "PC2 (12.6%)") +
  ggtitle(paste("PCoA of Species Abundance (BCD)"))  # Create PCoA plot with title
pcoa_plot1

# permanova
distance <-  phyloseq::distance(all, method = "bray")
meta <- as(sample_data(all), "data.frame")
perm_result <- adonis2(distance ~ location, data = meta)[1,]

# Separate blood sample
blood <- subset_samples(all, location == "blood")
blood <- phyloseq_filter_sample_wise_abund_trim(
  blood,  minabund = 5e-4)

bc_distance <- phyloseq::distance(blood, method = "bray")
pcoa_results <- ordinate(blood, method = "PCoA", distance = bc_distance)  # Perform PCoA

# Make Plot
pcoa_plot2 <- NULL
pcoa_plot2 <- plot_ordination(blood, pcoa_results, shape = 'age', color = 'sub.location') +
  geom_point(size = 3, alpha = 0.6) +
  geom_density_2d(aes(x = Axis.1, y = Axis.2), color = "grey", alpha = 0.8) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black")) +
  labs(x = "PC1 (46.1%)", y = "PC2 (25.3%)") +
  ggtitle("PCoA of Blood Sample Species Abundance (BCD)")
pcoa_plot2

# permanova
distance <-  phyloseq::distance(blood, method = "bray")
meta <- as(sample_data(blood), "data.frame")
perm_result <- adonis2(distance ~ sub.location, data = meta)[1,]

pcoa_plot <- ggarrange(pcoa_plot1, pcoa_plot2, ncol = 2, nrow = 1)

ggsave("./3B.png", pcoa_plot, width = 10, height = 3.5, bg = "white")



