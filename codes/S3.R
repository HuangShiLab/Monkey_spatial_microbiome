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
library(eulerr)

# --- helper: prevalence filter (≥50% by default) ---
filter_by_prevalence <- function(ps, prev = 0.50, detection = 0) {
  filter_taxa(
    ps,
    function(x) mean(x > detection, na.rm = TRUE) >= prev,
    prune = TRUE
  )
}

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
stomach <- subset_samples(bc_braken, location == "stomach")
upper <- subset_samples(bc_braken, group == "upper")
SI <- subset_samples(bc_braken, group == "small_intestine")
LI <- subset_samples(bc_braken, group == "large_intestine")

stomach   <- filter_by_prevalence(stomach,   prev = 0.50)
upper     <- filter_by_prevalence(upper,     prev = 0.50)
SI   <- filter_by_prevalence(SI,   prev = 0.50)
LI     <- filter_by_prevalence(LI,     prev = 0.50)


# Build a named list
sets_list <- list(
  upper = taxa_names(upper),
  stomach = taxa_names(stomach)
)

# Create euler object
fit <- euler(sets_list)

venn_pair <- function(ps_a, ps_b, name_a, name_b, title,
                      fills, alpha = 0.75) {
  set_a <- taxa_names(ps_a)
  set_b <- taxa_names(ps_b)
  fit <- euler(setNames(list(set_a, set_b), c(name_a, name_b)))
  plot(
    fit,
    quantities = TRUE,
    edges = TRUE,
    labels = list(font = 2),
    fills = list(fill = fills, alpha = alpha),
    main = title
  )
}

# Make the three pairwise plots (assumes your objects already filtered at 50% prev):
p_upper_stomach <- venn_pair(upper, stomach, "Upper", "Stomach",
                             "Upper ∩ Stomach",
                             fills = c("#9ecae1", "#fb9a99"))

p_stomach_SI    <- venn_pair(stomach, SI, "Stomach", "SI",
                             "Stomach ∩ SI",
                             fills = c("#fb9a99", "#fdd0a2"))

p_SI_LI         <- venn_pair(SI, LI, "SI", "LI",
                             "SI ∩ LI",
                             fills = c("#fdd0a2", "#fffa94"))


ggsave("./plots/S2-1.png", p_upper_stomach, width = 4, height = 4.5, bg = "white")
ggsave("./plots/S2-2.png", p_stomach_SI, width = 4, height = 4.5, bg = "white")
ggsave("./plots/S2-3.png", p_SI_LI, width = 4, height = 4.5, bg = "white")




