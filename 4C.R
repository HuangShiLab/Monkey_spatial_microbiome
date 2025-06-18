#######################################
# Monkey project
#
# Bile acids: microbial correlation with isoallo-LCA (absolute)
#
# Author: HOU Shuwen
#######################################

# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stats)
library(glmnet)
library(scales)

# 1. Read in data
setwd("~/Downloads")
rel_abund <- as.data.frame(read_tsv("monkey_data/species_abundance_bracken.tsv"))
rownames(rel_abund)=rel_abund[,1]
rel_abund=rel_abund[,-1]
metabolite <- as.data.frame(read_csv("bile acids/concentration.csv")) %>%
  filter(Metabolites != "sum")
rownames(metabolite)=metabolite[,1]
metabolite=metabolite[,-1][,-1]
qPCR_df <- read_excel("monkey_data/qPCR_copy.xlsx")

# filter out low abundance and low prevalance taxa
# Compute prevalence: proportion of samples where microbe is non-zero
prevalence <- rowSums(rel_abund > 0) / ncol(rel_abund)

# Compute mean relative abundance per microbe
mean_abundance <- rowMeans(rel_abund)

# Filter criteria
prevalence_threshold <- 0.1   # 10%
abundance_threshold <- 0.001  # Mean relative abundance > 0.001

# Filter microbes
keep <- which(prevalence > prevalence_threshold)
#              & mean_abundance > abundance_threshold)

# Filtered dataset
rel_abund <- rel_abund[keep, ]

load_vec <- setNames(qPCR_df$qPCR, qPCR_df$ID)

common_samples <- intersect(colnames(rel_abund), names(load_vec))

abs_abund <- sweep(
  rel_abund[, common_samples, drop = FALSE],
  MARGIN = 2,
  STATS  = load_vec[common_samples],
  FUN    = "*"
)
abs_abund_df <- as.data.frame(abs_abund, stringsAsFactors = FALSE)
abs_abund_df$species <- rownames(abs_abund_df)
abs_abund_df <- abs_abund_df[, c("species", setdiff(colnames(abs_abund_df), "species"))]
abs_abund_df<- abs_abund_df[, colnames(metabolite), drop = FALSE]

# 3. Choose one metabolite to test
metabolite_of_interest <- as.numeric(metabolite["Isoallolithocholic Acid", ])
names(metabolite_of_interest) <- colnames(metabolite)

# 4. Spearman correlation
cor_results <- t(
  apply(
    abs_abund_df,
    MARGIN = 1,
    FUN = function(x) {
      test <- cor.test(x, metabolite_of_interest, method = "pearson")
      c(r = as.numeric(test$estimate), p.value = test$p.value)
    }
  )
)

# Convert to data frame
cor_results <- as.data.frame(cor_results)
cor_results$p_adj <- p.adjust(cor_results$p.value, method = "fdr")

# 5. Filter significant correlations
sig_df <- cor_results %>% filter(p_adj < 0.05)

# Prepare data
species_name <- "Odoribacteraceae bacterium"
x_vals <- as.numeric(abs_abund_df[species_name, ])
y_vals <- as.numeric(metabolite_of_interest)  # Ensure sample alignment

# Data frame for plotting
data <- data.frame(
  SpeciesAbundance = x_vals,
  Metabolite = y_vals
)

# Calculate R-squared from linear model
fit <- cor.test(x_vals, y_vals, method = "pearson")
p <- fit$p.value


p1 <- ggplot(data, aes(x = SpeciesAbundance, y = Metabolite)) +
  geom_point(size = 3, color = "#bbdf27", alpha = 0.6) + 
  labs(title = "Odoribacteraceae bacterium\nCorrelation with Isoallo-LCA",
       x = "Microbial Absolute Abundance",
       y = "Isoallo-LCA Concentration") +
  geom_smooth(method = "lm", color = "#bbdf27", fill = "grey90", size = 1, alpha = 0.5) +
#  scale_x_log10() +
#  scale_y_log10() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 12),
        axis.line = element_line(size = 0.5, color = "black")) +
  geom_text(aes(x =700000, y = 2000000, label = "p = 0.002"),
            inherit.aes = FALSE,  # Avoid inheriting x/y mappings again
            color = "black", hjust = 0, vjust = 1, size = 4)

p1

species_name <- "Segatella copri"
x_vals <- as.numeric(abs_abund_df[species_name, ])
y_vals <- as.numeric(metabolite_of_interest)  # Ensure sample alignment

# Data frame for plotting
data <- data.frame(
  SpeciesAbundance = x_vals,
  Metabolite = y_vals
)

# Calculate R-squared from linear model
fit <- cor.test(x_vals, y_vals, method = "pearson")
p <- fit$p.value


p2 <- ggplot(data, aes(x = SpeciesAbundance, y = Metabolite)) +
  geom_point(size = 3, color = "skyblue", alpha = 0.6) + 
  labs(title = species_name,
       x = "Absolute Abundance",
       y = "Isoallo-LCA (µM)") +
  geom_smooth(method = "lm", color = "lightblue", fill = "grey90", size = 1, alpha = 0.5) +
  #  scale_x_log10() +
  #  scale_y_log10() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black")) +
  geom_text(aes(x =700000, y = 2000000, label = "p < 0.001"),
            inherit.aes = FALSE,  # Avoid inheriting x/y mappings again
            color = "black", hjust = 0, vjust = 1, size = 4)

p2

ggsave("./bile acids/4C.png", p1, width = 5, height = 4, bg = "white")
