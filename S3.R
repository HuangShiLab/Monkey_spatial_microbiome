#######################################
# Monkey project
#
# Alpha diversity -- oral and stomach Venn Diagram
#
# Author: HOU Shuwen
#######################################

# necessary packages
library(readr)
library(ggplot2)
library(phyloseq)
library(readxl)
library(dplyr)
library(ggpubr)
library(ggbeeswarm)
library(rstatix)
library(vegan)
library(metagMisc)
library(tidyr)
library(tibble)

class_color <- c(esophagus = "#eacc76", 
                 large_intestine = "#5c7272",
                 small_intestine = "#acab4b",
                 kidney = "#7c99bc", 
                 oral = "#a5b3c1",
                 stomach = "#c7a8a3", 
                 stool = "#e9b962")

### if I use 2bRAD input
# data input
setwd("~/Downloads/monkey_data")
bc_kraken <- as.data.frame(read.table("2b_species_abundance.txt", header = T))
metadata <- read.table("kraken_metadata.txt",sep="\t",header=TRUE)

# tidy up feature table
rownames(bc_kraken) <- bc_kraken$Species
### end

# Data input
setwd("~/Downloads/monkey_data")
bc_kraken <- as.data.frame(read_tsv("species_abundance_bracken.tsv",col_names = T))
rownames(bc_kraken) <- bc_kraken$name
bc_kraken <- bc_kraken[,-1]
bc_kraken[] <- lapply(bc_kraken, function(x) as.numeric(as.character(x)))
OTU <- otu_table(bc_kraken,taxa_are_rows = TRUE)

metadata <- read.table("kraken_metadata.txt",sep="\t",row.names=1,header=TRUE)
META <- sample_data(metadata)

# build phyloseq object
bc_braken <- phyloseq(OTU,META)

# filter kidney samples
oral <- subset_samples(bc_braken, location == "oral")
esophagus <- subset_samples(bc_braken, location == "esophagus")
stomach <- subset_samples(bc_braken, location == "stomach")
upper <- subset_samples(bc_braken, group == "upper")

oral <- filter_taxa(oral, function(x) sum(x) > 0, prune = TRUE)
esophagus <- filter_taxa(esophagus, function(x) sum(x) > 0, prune = TRUE)
stomach <- filter_taxa(stomach, function(x) sum(x) > 0, prune = TRUE)
upper <- filter_taxa(upper, function(x) sum(x) > 0, prune = TRUE)

library(eulerr)
# Build a named list
sets_list <- list(
  upper = taxa_names(upper),
  stomach = taxa_names(stomach)
)

# Create euler object
fit <- euler(sets_list)

# Plot
p <- plot(fit,
     fills = c("#a5b3c1", "#c7a8a3"),
     labels = list(font = 2),
     edges = TRUE,
     quantities = TRUE, 
     main = "Proportional Species Overlap")
ggsave("./plots/S3.png", p, width = 4, height = 4.5, bg = "white")





