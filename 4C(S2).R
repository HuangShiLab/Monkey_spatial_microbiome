#######################################
# Monkey project
#
# Bile acids: microbial correlation with isoallo-LCA (lasso)
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

# Prepare the predictor matrix (species as predictors)
X <- t(as.matrix(abs_abund_df))  # Samples x Species

# Prepare the response vector (metabolite)
y <- as.numeric(metabolite["Isoallolithocholic Acid", ])
names(y) <- colnames(metabolite)


# Ensure matching sample order
stopifnot(identical(rownames(X), names(y)))

# Run LASSO regression with cross-validation to choose lambda
set.seed(123)  # for reproducibility
cv_fit <- cv.glmnet(X, y, alpha = 1)  # alpha=1 for LASSO

# Extract the best lambda
best_lambda <- cv_fit$lambda.min

# Refit model with best lambda
lasso_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)

# Extract non-zero coefficients (selected species)
coeffs <- coef(lasso_model)
nonzero_coeffs <- coeffs[coeffs[, 1] != 0, , drop = FALSE]
selected_species <- rownames(coeffs)[-1]  # exclude intercept

# Output results
selected_df <- data.frame(
  species = selected_species,
  coefficient = as.numeric(coeffs[-1, , drop = FALSE])  # exclude intercept
)
