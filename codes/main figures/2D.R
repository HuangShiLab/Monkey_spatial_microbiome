#######################################
# Monkey project
#
# Pathway diversity
#
# Author: HOU Shuwen
#######################################

# Load required packages
library(ade4)
library(ggplot2)
library(ggpubr)
library(vegan)
library(readr)
library(readxl)
library(dplyr)
library(rstatix)
library(phyloseq)


# ---------- Pathway-based Beta Diversity (Bray-Curtis) ----------

# Set working directory
setwd("~/Downloads/monkey_data")

# Load unstratified pathway abundance table
pathway <- read.table("merged_unstratified_pathway.tsv", 
                      sep = '\t', row.names = 1, header = TRUE)

# Convert all values to numeric
pathway[] <- lapply(pathway, function(x) as.numeric(as.character(x)))

# Normalize each sample (column) to relative abundance
pathway <- sweep(pathway, 2, colSums(pathway), FUN = "/")

# Remove one sample (excluded elsewhere)
pathway <- pathway[, colnames(pathway) != "Y1_SI_113cm"]

# Load metadata
metadata <- read.table("kraken_metadata.txt", sep = "\t", header = TRUE)

# # Split into age groups
# pathway_old <- pathway[1:34]
# pathway_young <- pathway[35:69]

# Calculate Bray-Curtis distances across all samples
pathway_dist <- vegdist(t(pathway), method = "bray")

# Convert distance matrix to long format
dist_long <- as.data.frame(as.table(as.matrix(pathway_dist))) %>%
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

# Merge metadata
dist_long <- left_join(dist_long, metadata %>% select(Sample1 = ID, location1 = location))
dist_long <- left_join(dist_long, metadata %>% select(Sample2 = ID, location2 = location))

# Keep only pairs involving oral samples
dist_long_oral <- dist_long %>%
  filter(location1 == "oral" | location2 == "oral")

# Swap columns to ensure Sample1 is oral
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

# Wilcoxon test comparing oral to other locations
stat_test <- dist_long_oral %>%
  wilcox_test(as.formula("Distance ~ location2")) %>%
  add_xy_position(x = "location2") %>% 
  filter(xmin == "1") %>%
  add_significance("p") %>%
  mutate(location2 = group2)

# Set factor levels for plotting
dist_long_oral$location2 <- factor(dist_long_oral$location2, 
                                   levels = c("oral", "esophagus", "stomach", 
                                              "small_intestine", "large_intestine"))

# Plot: Bray-Curtis distance to oral
p1 <- ggplot(dist_long_oral, aes(x = location2, y = Distance, fill = location2)) +
  geom_boxplot(alpha = 0.8) +
  labs(title = "Community Functional\nDissimilarity to Oral", x = "location", y = "Bray-Curtis dissimilarity") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.line = element_line(size = 0.5, color = "black")) +
  scale_fill_viridis_d(option = "plasma", begin = 0, end = 1, direction = 1) +
  stat_pvalue_manual(stat_test, label = "{p.signif}", tip.length = 0)

# ---------- Separate Young and Old Samples ----------

# # Calculate Bray-Curtis distances separately
# pathway_old_dist <- vegdist(t(pathway_old), method = "bray")
# pathway_young_dist <- vegdist(t(pathway_young), method = "bray")
# 
# # Convert to long format and label by age
# dist_long_old <- as.data.frame(as.table(as.matrix(pathway_old_dist))) %>%
#   mutate(Sample1 = Var1, Sample2 = Var2, Distance = Freq, age = "old") %>%
#   filter(Sample1 != Sample2)
# 
# dist_long_young <- as.data.frame(as.table(as.matrix(pathway_young_dist))) %>%
#   mutate(Sample1 = Var1, Sample2 = Var2, Distance = Freq, age = "young") %>%
#   filter(Sample1 != Sample2)
# 
# # Combine old and young
# dist_long <- rbind(dist_long_old, dist_long_young)
# 
# # Remove duplicate pairs
# dist_long <- dist_long %>%
#   mutate(
#     Sample1 = as.character(Sample1),
#     Sample2 = as.character(Sample2),
#     Pair = paste(pmin(Sample1, Sample2), pmax(Sample1, Sample2), sep = "-")
#   ) %>% distinct(Pair, .keep_all = TRUE) %>%
#   select(-Pair)
# 
# # Merge with metadata
# dist_long <- left_join(dist_long, metadata %>% select(Sample1 = ID, location1 = location))
# dist_long <- left_join(dist_long, metadata %>% select(Sample2 = ID, location2 = location))
# 
# # Filter to oral-related comparisons
# dist_long_oral <- dist_long %>% filter(location1 == "oral" | location2 == "oral")
# 
# # Reorganize to have oral in Sample1
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
#                                    levels = c("oral", "esophagus", "stomach", "small_intestine", "large_intestine", "stool"))
# dist_long_oral$age <- factor(dist_long_oral$age, levels = c("young", "old"))
# 
# # Wilcoxon test by age group
# stat_test <- dist_long_oral %>%
#   group_by(location2) %>%
#   wilcox_test(as.formula("Distance ~ age")) %>%
#   add_xy_position(x = "location2") %>%
#   add_significance("p")
# 
# # Adjust label positions
# stat_test$xmin <- stat_test$xmin + stat_test$x - 1
# stat_test$xmax <- stat_test$xmax + stat_test$x
# 
# # Plot: Distance to oral, by age group
# p2 <- ggplot(dist_long_oral, aes(x = interaction(age, location2), y = Distance, fill = location2)) +
#   geom_boxplot() +
#   labs(title = "Distance to oral (split by age group)",
#        x = "location x age", y = "Bray Curtis Dissimilarity") +
#   theme(panel.background = element_rect(fill = "white", color = NA),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.title = element_text(size = 16),
#         axis.title = element_text(size = 14),
#         axis.text.x = element_text(angle = 0, hjust = 0.5),
#         axis.line = element_line(size = 0.5, color = "black")) +
#   scale_x_discrete(labels = rep(levels(dist_long_oral$age), length(levels(dist_long_oral$location2)))) +
#   stat_pvalue_manual(stat_test, label = "{p.signif}")

# ---------- Alpha Diversity Analysis ----------

# Re-import pathway data for alpha diversity
setwd("~/Downloads/monkey_data")
pathway <- read.table("merged_unstratified_pathway.tsv", sep = '\t', row.names = 1, header = TRUE)
pathway[] <- lapply(pathway, function(x) as.numeric(as.character(x)))
pathway <- sweep(pathway, 2, colSums(pathway), FUN = "/")

# Load and filter metadata
metadata <- read.table("kraken_metadata.txt", sep = "\t", header = TRUE, row.names = 1)
metadata <- metadata[rownames(metadata) != "Y1_SI_113cm", ]

# Construct phyloseq object
OTU <- otu_table(pathway, taxa_are_rows = TRUE)
META <- sample_data(metadata)
pathway <- phyloseq(OTU, META)

# Filter kidney samples
pathway <- subset_samples(pathway, location != "kidney")

# Estimate alpha diversity (Shannon index)
measures <- "Shannon"
alpha_diversity <- estimate_richness(pathway, measures = measures)

# Merge with sample data
df <- merge(alpha_diversity, pathway@sam_data, by = "row.names")
names(df)[1] <- "Sample"
df$location <- factor(df$location, levels = c("oral", "esophagus", "stomach", 
                                              "small_intestine", "large_intestine"))

# Wilcoxon test: diversity ~ location
stat_test <- df %>%
  wilcox_test(as.formula(paste(measures, "~ location"))) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p") %>%
  add_xy_position(x = "location2") %>%
  filter(group1 == "stomach")
stat_test[1, "xmin"] <- 3
stat_test[2, "xmin"] <- 3
stat_test[1, "xmax"] <- 4
stat_test[2, "xmax"] <- 5

# Plot: Alpha diversity by location
p0 <- ggplot(df, aes_string(x = "location", y = measures, fill = "location")) +
  geom_boxplot(alpha = 0.8) +
  labs(y = paste(measures, "index"), x = "location", title = "Functional Diversity") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.line = element_line(size = 0.5, color = "black"),
        legend.justification = c(0, 1)) +
  scale_fill_viridis_d(option = "plasma", begin = 0, end = 1, direction = 1) +
  scale_x_discrete(labels = function(x) sapply(strsplit(x, "\\."), `[`, 1)) +
  stat_pvalue_manual(stat_test, label = "{p.adj.signif}", tip.length = 0, inherit.aes = FALSE)

# ---------- Combine All Plots ----------

# Arrange all plots side-by-side
distance_plot <- ggarrange(p0, p1, ncol = 2, nrow = 1,
                           common.legend = TRUE, legend = "right")

# Display combined plot
distance_plot

# Save figure
ggsave("./plots/2A.pdf", distance_plot, width = 7.8, height = 4.5)
