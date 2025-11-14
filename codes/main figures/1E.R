#######################################
# Monkey Study
#
# Human monkey comparison -- Part III
# Differential abundance
#
# Author: HOU Shuwen
#######################################

# loading packages
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(openxlsx)
library(metagMisc)
library(ANCOMBC)
library(ggbeeswarm)
library(ggpubr)
library(tibble)
library(microbiome)
library(patchwork)
library(writexl)

# Set working directory and read metadata
setwd("~/Downloads/nature2023_data/")
md <- read_excel("meta_WMS.xlsx")

# Read species abundance feature table
feature_table <- read_tsv("species_abundance_bracken.tsv", col_names = TRUE)

# Filter metadata for specific subject IDs
md_filter <- md %>% filter(subject_id %in% c(6, 8, 12, 14))
SRR <- md_filter$Run  # Extract selected sample run IDs

# Select relevant columns from feature table
filtered_otu <- feature_table %>% select(name, all_of(SRR))
filtered_otu$name <- gsub("'", "", filtered_otu$name)  # Remove single quotes from species names

# Set working directory to access another dataset
setwd("~/Downloads/monkey_data")
bc_kraken <- as.data.frame(read_tsv("species_abundance_bracken.tsv", col_names = TRUE))
bc_kraken$name <- gsub("'", "", bc_kraken$name)  # Remove single quotes from species names

# Merge the two feature tables by species name
bc_kraken <- full_join(filtered_otu, bc_kraken, by = "name")
bc_kraken <- as.data.frame(bc_kraken)
bc_kraken[is.na(bc_kraken)] <- 0  # Replace NA values with 0

# Convert species names into row names and normalize data
rownames(bc_kraken) <- bc_kraken$name
bc_kraken <- bc_kraken[, -1]  # Remove the 'name' column
bc_kraken[] <- lapply(bc_kraken, function(x) as.numeric(as.character(x)))
bc_kraken <- bc_kraken %>% mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
bc_kraken %>% summarise(across(everything(), sum, na.rm = TRUE))  # Summarize total abundance per sample

# Convert to Phyloseq OTU table
OTU <- otu_table(bc_kraken, taxa_are_rows = TRUE)

# Read metadata and create Phyloseq sample data
setwd("~/Downloads/comparison_data")
metadata <- as.data.frame(read_excel("meta.xlsx"))
rownames(metadata) <- metadata[[1]]  # Set first column as row names
metadata[[1]] <- NULL  # Remove the redundant first column
META <- sample_data(metadata)

# Build the Phyloseq object
all <- phyloseq(OTU, META)

# Remove samples from 'kidney' location
all <- subset_samples(all, location != "kidney")

# Create a taxonomy table
temp <- rownames(as.data.frame(all@otu_table))
TAX <- data.frame(Species = temp)
TAX$Species <- paste0("s__", TAX$Species)  # Add species prefix
TAX <- tax_table(TAX)
row.names(TAX) <- temp
colnames(TAX) <- c("Species")
all@tax_table <- TAX  # Assign taxonomy table to Phyloseq object

# Subset Phyloseq object by sample type
LI <- subset_samples(all, location == "large intestine")
SI <- subset_samples(all, location == "small intestine")
intestine <- subset_samples(all, location == "small intestine"|location == "large intestine")
oral <- subset_samples(all, location == "oral")
upper <- subset_samples(all, location == "oral" | location == "stomach")

# Function to process and analyze Phyloseq objects
analyze_phyloseq <- function(phyloseq_obj, group_var = "species") {
  # 1. Pseudocount
  mat <- as(phyloseq_obj@otu_table, "matrix")
  counts <- round(mat * 1e6)
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  storage.mode(counts) <- "integer"

  
  # 2. Run ANCOM-BC2
  out <- ancombc2(
    data       = counts,
    taxa_are_rows = TRUE,
    meta_data = as(sample_data(phyloseq_obj), "data.frame"),
    fix_formula = group_var
  )
  
  # 3. Extract results for the group variable
  # (lfc, se, W, p, q, diff for the "monkey vs human" contrast)
  res <- out$res %>%
    select(
      taxon,
      lfc_speciesmonkey, se_speciesmonkey,
      W_speciesmonkey, p_speciesmonkey, q_speciesmonkey,
      diff_robust_speciesmonkey
    ) %>%
    mutate(
      beta     = lfc_speciesmonkey,  # effect size
      direct   = case_when(
        diff_robust_speciesmonkey & beta > 0  ~ "Monkey",
        diff_robust_speciesmonkey & beta <= 0 ~ "Human",
        TRUE                           ~ "Not Significant"
      )
    ) %>%
    filter(direct != "Not Significant") %>%
    arrange(q_speciesmonkey)
  
  # 4. Return results
  list(
    ancom = as.data.frame(res)
  )
}


# List of Phyloseq objects for different locations
phyloseq_objects <- list(all = all, LI = LI, SI = SI, upper = upper, oral = oral)

# Apply analysis function to each Phyloseq object
results <- lapply(phyloseq_objects, analyze_phyloseq)

# Separate results into LEfSe and ANCOM tables
#lefse_results <- list()
ancom_results <- list()

# Iterate over results and store in lists
for (location in names(results)) {
#  lefse_table <- as.data.frame(results[[location]]$lefse)
  ancom_table <- as.data.frame(results[[location]]$ancom)
  
#  if (nrow(lefse_table) > 0) {
#    lefse_table$Location <- location
#    lefse_results[[location]] <- lefse_table
#  }
  
  if (nrow(ancom_table) > 0) {
    ancom_table$Location <- location
    ancom_results[[location]] <- ancom_table
  }
}

# Combine results into final data frames
#lefse_combined <- do.call(rbind, lefse_results)
ancom_combined <- do.call(rbind, ancom_results)

# Write results to an Excel file
write_xlsx(ancom_combined, path = "./microbiome_results.xlsx")

## box plot
# upper
upper_otu <- as.data.frame(upper@otu_table)
upper_meta <- as.data.frame(upper@sam_data)

species_names <- c("Fusobacterium polymorphum",
                   "Prevotella intermedia")

# Convert OTU table from wide to long format
otu_long <- upper_otu %>%
  rownames_to_column("Species") %>%
  filter(Species %in% species_names) %>%
  pivot_longer(-Species, names_to = "SampleID", values_to = "Abundance")

# Add SampleID as a column to metadata
upper_meta$SampleID <- rownames(upper_meta)
colnames(upper_meta)[colnames(upper_meta) == "species"] <- "host"

# Join OTU data with metadata
otu_long <- otu_long %>%
  left_join(upper_meta[, c("SampleID", "host")], by = "SampleID")

# Final data transformation
otu_long <- otu_long %>%
  mutate(
    Abundance = log10(as.numeric(Abundance) + 1),
    host = as.factor(host),
    Tract = "upper GI tract"
  )

# Plot with facet_wrap
p1 <- ggplot(otu_long, aes(x = host, y = Abundance, fill = host)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(alpha = 0.5, color = "black", size = 1, cex = 0.5) +
  scale_fill_manual(values = c("#eacc76", "#a76b3e")) +
  theme_minimal() +
  labs(
    x = "", y = "log(Relative abundance + 1)",title = "Upper GI tract",
    color = "Host"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(), 
    axis.text.y = element_text(angle = 90),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(color = "black", face = "italic"),
    legend.position = "none",
    strip.placement = "outside"
  ) +
  facet_wrap(~Species, nrow = 1, strip.position = "left")
p1


# intestine
intestine_otu <- as.data.frame(intestine@otu_table)
intestine_meta <- as.data.frame(intestine@sam_data)
species_names <- c("Segatella copri")

# Convert OTU table from wide to long format
otu_long <- intestine_otu %>%
  rownames_to_column("Species") %>%
  filter(Species %in% species_names) %>%
  pivot_longer(-Species, names_to = "SampleID", values_to = "Abundance")

# Add SampleID as a column to metadata
intestine_meta$SampleID <- rownames(intestine_meta)
colnames(intestine_meta)[colnames(intestine_meta) == "species"] <- "host"

# Join OTU data with metadata
otu_long <- otu_long %>%
  left_join(intestine_meta[, c("SampleID", "host")], by = "SampleID")

# Final data transformation
otu_long <- otu_long %>%
  mutate(
    Abundance = log10(as.numeric(Abundance) + 1),
    host = as.factor(host),
    Tract = "intestine"
  )

# Plot with facet_wrap
p2 <- ggplot(otu_long, aes(x = host, y = Abundance, fill = host)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(alpha = 0.5, color = "black", size = 1, cex = 0.5) +
  scale_fill_manual(values = c("#eacc76", "#a76b3e")) +
  theme_minimal() +
  labs(
    x = "", y = "",title = "Intestines",
    color = "Host"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(angle = 90),
    axis.ticks.y = element_line(), 
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(color = "black", face = "italic"),
    strip.placement = "outside"
  ) +
  facet_wrap(~Species, nrow = 1, strip.position = "left")
p2

species_names <- c("Treponema succinifaciens",
                   "Anaerobutyricum hallii", "Anaerostipes hadrus",
                   "Bifidobacterium adolescentis","Bifidobacterium longum")

# Convert OTU table from wide to long format
otu_long <- intestine_otu %>%
  rownames_to_column("Species") %>%
  filter(Species %in% species_names) %>%
  pivot_longer(-Species, names_to = "SampleID", values_to = "Abundance")


# Add SampleID as a column to metadata
intestine_meta$SampleID <- rownames(intestine_meta)
colnames(intestine_meta)[colnames(intestine_meta) == "species"] <- "host"

# Join OTU data with metadata
otu_long <- otu_long %>%
  left_join(intestine_meta[, c("SampleID", "host")], by = "SampleID")

# Final data transformation
otu_long <- otu_long %>%
  mutate(
    Abundance = log10(as.numeric(Abundance) + 1),
    host = as.factor(host),
    Tract = "intestine"
  )
otu_long$Species <- factor(otu_long$Species, levels = species_names)

# Plot with facet_wrap
p3 <- ggplot(otu_long, aes(x = host, y = Abundance, fill = host)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(alpha = 0.5, color = "black", size = 1, cex = 0.5) +
  scale_fill_manual(values = c("#eacc76", "#a76b3e")) +
  theme_minimal() +
  labs(
    x = "", y = "",
    color = "Host"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(), 
    axis.text.y = element_text(angle = 90),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(color = "black", face = "italic"),
    strip.placement = "outside"
  ) +
  facet_wrap(~Species, nrow = 1, strip.position = "left")+
  coord_cartesian(ylim = c(0, 0.03))
p3

combined_plot <- p1 + p2 + p3 + plot_layout(widths = c(2,0.9,6))
combined_plot

ggsave("./1C.png", combined_plot, width = 10, height = 2.5, bg = "white")



