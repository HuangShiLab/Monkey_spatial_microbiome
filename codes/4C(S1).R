#######################################
# Monkey project
#
# Bile acids: microbial correlation with isoallo-LCA (presence-absence)
#
# Author: HOU Shuwen
#######################################

# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stats)
# 1. Read in data
setwd("~/Downloads")
rel_abund <- as.data.frame(read_tsv("monkey_data/species_abundance_bracken.tsv"))
rownames(rel_abund) <- rel_abund[,1]
rel_abund <- rel_abund[,-1]

metabolite <- as.data.frame(read_csv("bile acids/concentration.csv")) %>%
  filter(Metabolites != "sum")
rownames(metabolite) <- metabolite[,1]
metabolite <- metabolite[,-1][,-1]

qPCR_df <- read_excel("monkey_data/qPCR_copy.xlsx")

# Convert to absence/presence (1 if > 0, else 0)
presence_absence <- (rel_abund > 0) * 1

# Compute prevalence
prevalence <- rowSums(presence_absence) / ncol(presence_absence)

# Compute mean "abundance" in binary data (same as prevalence)
# Filter criteria
prevalence_threshold <- 0.1  # 10%

# Filter microbes
keep <- which(prevalence > prevalence_threshold)

# Filtered dataset
presence_absence <- presence_absence[keep, ]

# Ensure samples match between presence_absence and metabolite
common_samples <- intersect(colnames(presence_absence), colnames(metabolite))
presence_absence <- presence_absence[, common_samples, drop = FALSE]

# 3. Choose one metabolite to test
metabolite_of_interest <- as.numeric(metabolite["Isoallolithocholic Acid", common_samples])
names(metabolite_of_interest) <- common_samples

# 4. Spearman correlation with binary presence/absence
cor_results <- t(
  apply(
    presence_absence,
    MARGIN = 1,
    FUN = function(x) {
      test <- cor.test(x, metabolite_of_interest, method = "spearman")
      c(rho = as.numeric(test$estimate), p.value = test$p.value)
    }
  )
)

# Convert to data frame
cor_results <- as.data.frame(cor_results)
cor_results$p_adj <- p.adjust(cor_results$p.value, method = "fdr")

# 5. Filter significant correlations
sig_df <- cor_results %>% filter(p_adj < 0.05)

# Add species names as a column
cor_results$species <- rownames(cor_results)

# Volcano plot
ggplot(cor_results, aes(x = rho, y = -log10(p_adj))) +
  geom_point(aes(color = p_adj < 0.05)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Microbial Correlation with Isoallolithocholic Acid",
    x = "Spearman’s rho",
    y = "-log10(FDR-adjusted p-value)",
    color = "FDR < 0.05"
  )

