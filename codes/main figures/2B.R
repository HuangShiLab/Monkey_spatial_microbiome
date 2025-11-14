#######################################
# Monkey project
#
# Species diversity analysis
#
# Author: HOU Shuwen
#######################################

# Load necessary R packages
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

### Optional: Load 2bRAD input instead of Bracken
# Set working directory
setwd("~/Downloads/monkey_data")

# Read 2bRAD species abundance data
bc_kraken <- as.data.frame(read.table("2b_species_abundance.txt", header = T))

# Load sample metadata
metadata <- read.table("kraken_metadata.txt", sep = "\t", header = TRUE)

# Assign species names as row names
rownames(bc_kraken) <- bc_kraken$Species
### End 2bRAD block

# Set working directory again (if not already set)
setwd("~/Downloads/monkey_data")

# Read species abundance data (Bracken output)
bc_kraken <- as.data.frame(read_tsv("species_abundance_bracken.tsv", col_names = T))

# Remove single quotes from sample names
bc_kraken$name <- gsub("'", "", bc_kraken$name)

# Load metadata
metadata <- read.table("kraken_metadata.txt", sep = "\t", header = TRUE)

# Set sample names as row names
rownames(bc_kraken) <- bc_kraken$name
bc_kraken <- bc_kraken[, -1]  # Remove name column

# Convert data to numeric
bc_kraken[] <- lapply(bc_kraken, function(x) as.numeric(as.character(x)))

# Normalize feature table by column sums
bc_kraken <- sweep(bc_kraken, 2, colSums(bc_kraken), FUN = "/")

# # Split data into old and young groups
# bc_kraken_old <- bc_kraken[1:34]
# bc_kraken_young <- bc_kraken[35:70]

# Create phyloseq object
OTU <- otu_table(bc_kraken, taxa_are_rows = TRUE)
META <- sample_data(metadata)
rownames(META) <- META$ID
bc_braken <- phyloseq(OTU, META)

# Remove kidney samples
bc_braken <- subset_samples(bc_braken, location != "kidney")

# Compute Bray-Curtis distance matrix
bc_kraken_dist <- vegdist(t(bc_kraken), method = "bray")

# Convert distance matrix to long format
dist_long <- as.data.frame(as.table(as.matrix(bc_kraken_dist))) %>%
  mutate(Sample1 = Var1, Sample2 = Var2, Distance = Freq) %>%
  filter(Sample1 != Sample2)

# Remove duplicated sample pairs
dist_long <- dist_long %>% 
  mutate(
    Sample1 = as.character(Sample1),
    Sample2 = as.character(Sample2),
    Pair = paste(pmin(Sample1, Sample2), pmax(Sample1, Sample2), sep = "-")
  ) %>% distinct(Pair, .keep_all = TRUE) %>%
  select(-Pair)

# Merge metadata for both samples
dist_long <- left_join(dist_long, metadata %>% select(Sample1 = ID, location1 = location))
dist_long <- left_join(dist_long, metadata %>% select(Sample2 = ID, location2 = location))

# Subset pairs with at least one oral sample
dist_long_oral <- dist_long %>% filter(location1 == "oral" | location2 == "oral")

# Swap sample positions to always have oral as Sample1
dist_long_oral <- dist_long_oral %>%
  mutate(
    TempSample1 = if_else(location2 == "oral", Sample2, Sample1),
    TempSample2 = if_else(location2 == "oral", Sample1, Sample2),
    TempLocation1 = if_else(location2 == "oral", location2, location1),
    TempLocation2 = if_else(location2 == "oral", location1, location2)
  ) %>%
  mutate(
    Sample1 = TempSample1,
    Sample2 = TempSample2,
    location1 = TempLocation1,
    location2 = TempLocation2
  ) %>%
  select(-TempSample1, -TempSample2, -TempLocation1, -TempLocation2)

# Remove kidney samples
dist_long_oral <- dist_long_oral %>% filter(location2 != "kidney")

# Statistical test: distance to oral vs. other sites
stat_test <- dist_long_oral %>%
  wilcox_test(as.formula("Distance ~ location2")) %>%
  add_xy_position(x = "location2") %>% 
  filter(xmin == "1") %>%
  add_significance("p") %>%
  mutate(location2 = group2)

# Set location order for plotting
dist_long_oral$location2 <- factor(dist_long_oral$location2,
                                   levels = c("oral", "esophagus", "stomach", "small_intestine", "large_intestine"))

# Boxplot: Distance to oral
p1 <- ggplot(dist_long_oral, aes(x = location2, y = Distance, fill = location2)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_viridis_d(option = "plasma", begin = 0, end = 1, direction = 1) +
  labs(title = "Community Compositional\nDissimilarity to Oral", x = "location", y = "Bray-Curtis dissimilarity") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.line = element_line(size = 0.5, color = "black")) +
  stat_pvalue_manual(stat_test, label = "{p.signif}", tip.length = 0)

# ------- Separate by age group (old vs. young) -------
# Compute Bray-Curtis distances for age groups
#bc_kraken_old_dist <- vegdist(t(bc_kraken_old), method = "bray")
#bc_kraken_young_dist <- vegdist(t(bc_kraken_young), method = "bray")

# Convert to long format
#dist_long_old <- as.data.frame(as.table(as.matrix(bc_kraken_old_dist))) %>%
#  mutate(Sample1 = Var1, Sample2 = Var2, Distance = Freq) %>%
#  filter(Sample1 != Sample2) %>% mutate(age = "old")

# dist_long_young <- as.data.frame(as.table(as.matrix(bc_kraken_young_dist))) %>%
#   mutate(Sample1 = Var1, Sample2 = Var2, Distance = Freq) %>%
#   filter(Sample1 != Sample2) %>% mutate(age = "young")

# Combine both age groups
# dist_long <- rbind(dist_long_old, dist_long_young)

# Remove duplicate pairs
# dist_long <- dist_long %>%
#   mutate(
#     Sample1 = as.character(Sample1),
#     Sample2 = as.character(Sample2),
#     Pair = paste(pmin(Sample1, Sample2), pmax(Sample1, Sample2), sep = "-")
#   ) %>% distinct(Pair, .keep_all = TRUE) %>%
#   select(-Pair)

# Add metadata
# dist_long <- left_join(dist_long, metadata %>% select(Sample1 = ID, location1 = location))
# dist_long <- left_join(dist_long, metadata %>% select(Sample2 = ID, location2 = location))

# Subset to samples involving oral site
# dist_long_oral <- dist_long %>% filter(location1 == "oral" | location2 == "oral")

# Swap positions to always have oral as Sample1
# dist_long_oral <- dist_long_oral %>%
#   mutate(
#     TempSample1 = if_else(location2 == "oral", Sample2, Sample1),
#     TempSample2 = if_else(location2 == "oral", Sample1, Sample2),
#     TempLocation1 = if_else(location2 == "oral", location2, location1),
#     TempLocation2 = if_else(location2 == "oral", location1, location2)
#   ) %>%
#   mutate(
#     Sample1 = TempSample1,
#     Sample2 = TempSample2,
#     location1 = TempLocation1,
#     location2 = TempLocation2
#   ) %>%
#   select(-TempSample1, -TempSample2, -TempLocation1, -TempLocation2)
# 
# # Remove kidney samples
# dist_long_oral <- dist_long_oral %>% filter(location2 != "kidney")
# 
# # Set factors for plotting
# dist_long_oral$location2 <- factor(dist_long_oral$location2,
#                                    levels = c("oral", "esophagus", "stomach", "small_intestine", "large_intestine"))
# dist_long_oral$age <- factor(dist_long_oral$age, levels = c("young", "old"))
# 
# # Statistical test: distance vs. age within each site
# stat_test <- dist_long_oral %>%
#   group_by(location2) %>%
#   wilcox_test(as.formula("Distance ~ age")) %>%
#   add_xy_position(x = "location2") %>%
#   add_significance("p")
# 
# # Adjust label positions for annotations
# stat_test$xmin <- stat_test$xmin + stat_test$x - 1
# stat_test$xmax <- stat_test$xmax + stat_test$x

# Boxplot: Distance by age group
#p2 <- ggplot(dist_long_oral, aes(x = interaction(age, location2), y = Distance, fill = location2)) +
#  geom_boxplot() +
#  labs(title = "Distance to oral (split by age group)", x = "location x age", y = "Bray Curtis Dissimilarity") +
#  theme(panel.background = element_rect(fill = "white", color = NA),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        plot.title = element_text(size = 16),
#        axis.title = element_text(size = 14),
#        axis.text.x = element_text(angle = 0, hjust = 0.5),
#        axis.line = element_line(size = 0.5, color = "black")) +
#  scale_x_discrete(labels = rep(levels(dist_long_oral$age), length(levels(dist_long_oral$location2)))) +
#  stat_pvalue_manual(stat_test, label = "{p.signif}")

# ----- Alpha diversity analysis -----
# Calculate Shannon index
measures <- "Shannon"
alpha_diversity <- estimate_richness(bc_braken, measures = measures)

# Merge with sample metadata
df <- merge(alpha_diversity, bc_braken@sam_data, by = "row.names")
names(df)[1] <- "Sample"
df$location <- factor(df$location, levels = c("oral", "esophagus", "stomach", 
                                              "small_intestine", "large_intestine"))

# Test: alpha diversity ~ age by location
# stat_test <- df %>%
#   group_by(location) %>%
#   wilcox_test(as.formula(paste(measures, "~ age"))) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance() %>%
#   add_xy_position(x = "location2")

# Boxplot: Shannon diversity by age
#alpha_diversity_plot <- ggplot(df, aes_string(x = "interaction(age, location)", y = measures)) +
#  geom_boxplot(aes(fill = location), alpha = 0.8) +
#  labs(y = paste(measures, "diversity"), x = "location") +
#  theme(panel.background = element_rect(fill = "white", color = NA),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        plot.title = element_text(size = 16),
#        axis.title = element_text(size = 14),
#        axis.line = element_line(size = 0.5, color = "black")) +
#  scale_x_discrete(labels = function(x) sapply(strsplit(x, "\\."), `[`, 1))

# Save alpha diversity plot
#ggsave("./plots/S2.png", alpha_diversity_plot, width = 9, height = 4.5, bg = "white")

# Test: alpha diversity ~ location (combining age groups)
stat_test <- df %>%
  wilcox_test(as.formula(paste(measures, "~ location"))) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  filter(p.adj.signif=="**") %>%
  add_xy_position(x = "location2")
stat_test[1, "xmin"] <- 3
stat_test[1, "xmax"] <- 4

# Boxplot: Shannon diversity by location
p0 <- ggplot(df, aes_string(x = "location", y = measures, fill = "location")) +
  geom_boxplot(alpha = 0.8) +
  labs(y = paste(measures, "index"), x = "location", title = "Compositional Diversity") +
  scale_fill_viridis_d(option = "plasma", begin = 0, end = 1, direction = 1) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.line = element_line(size = 0.5, color = "black"),
        legend.justification = c(0, 1)) +
  scale_x_discrete(labels = function(x) sapply(strsplit(x, "\\."), `[`, 1)) +
  stat_pvalue_manual(stat_test, label = "{p.adj.signif}", tip.length = 0, inherit.aes = FALSE)

# Combine all plots into one figure
distance_plot <- ggarrange(p0, p1, ncol = 2, nrow = 1,
                           common.legend = TRUE, legend = "right")

# Save final figure
ggsave("./plots/2B.pdf", distance_plot, width = 7.8, height = 4.5, bg = "transparent")
