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
library(purrr)
library(reshape2)
library(phyloseq)
library(palettes)
library(vegan)
library(openxlsx)
library(patchwork)
library(metagMisc)

# data input
setwd("~/Downloads/monkey_data")
content_species <- as.data.frame(read.table("2b_species_abundance.txt", header = TRUE))

setwd("~/Downloads/monkey_mucosa")
feature_table <- as.data.frame(read.table("species_abundance.txt", header = TRUE))

# join
feature_table <- full_join(feature_table, content_species, by = "Species")
feature_table[is.na(feature_table)] <- 0

metadata <- read.table("meta_all.txt",sep="\t",header=TRUE)
ID <- metadata$ID

# tidy up feature table
rownames(feature_table) <- feature_table$Species
feature_table <- feature_table[,-1]
feature_table[] <- lapply(feature_table, function(x) as.numeric(as.character(x)))


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
pcoa_plot1 <- plot_ordination(all, pcoa_results, color = 'location',shape = 'part') +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = c("darkred", "#eeed89", "#f0aa73", "#d66b93")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black"))+
  labs(x = "PC1 (27%)", y = "PC2 (12.6%)") +
  ggtitle(paste("PCoA of All Samples Species Abundance (BCD)"))  # Create PCoA plot with title
pcoa_plot1

# permanova
distance <-  phyloseq::distance(all, method = "bray")
meta <- as(sample_data(all), "data.frame")
perm_result1 <- adonis2(distance ~ location, data = meta)[1,]
perm_result2 <- adonis2(distance ~ part, data = meta)[1,]

# Separate blood sample
blood <- subset_samples(all, location == "blood")
blood <- phyloseq_filter_sample_wise_abund_trim(
  blood,  minabund = 5e-4)

bc_distance <- phyloseq::distance(blood, method = "bray")
pcoa_results <- ordinate(blood, method = "PCoA", distance = bc_distance)  # Perform PCoA

# Make Plot
pcoa_plot2 <- NULL
pcoa_plot2 <- plot_ordination(blood, pcoa_results,  color = 'sub.location') +
  geom_point(size = 4, alpha = 0.7) +
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

# Supplementary
# Permanova heatmap
# 1) Build a combined grouping column on the full object
meta <- as(sample_data(all), "data.frame")

# If it's blood, keep it as its own group; otherwise combine "location+part"
meta$pair_group <- ifelse(meta$location == "blood",
                          "blood",
                          paste(meta$location, meta$part, sep = "_"))

# write back to phyloseq
sample_data(all)$pair_group <- meta$pair_group

# 2) Make the site↔site pairwise PERMANOVA heatmap
heat_pairwise <- pairwise_permanova_heatmap(
  phy = all,
  group_var = "pair_group",   # <— use the combined factor
  distance = "bray",
  permutations = 999,
  value_to_plot = "R2",       # or "p" to show adjusted p-values
  p_adjust = "BH",
  alpha = 0.05,
  title = "Pairwise PERMANOVA R2"
)
heat_pairwise
heat_blood <- pairwise_permanova_heatmap(
  phy = blood,
  group_var = "sub.location",   # <— use the combined factor
  distance = "bray",
  permutations = 999,
  value_to_plot = "R2",       # or "p" to show adjusted p-values
  p_adjust = "BH",
  alpha = 0.05,
  title = "Blood PERMANOVA R2"
)
heat_blood

heatmap_plot <- ggarrange(heat_pairwise, heat_blood, ncol = 2, nrow = 1, widths = c(3, 1))

ggsave("./3BS.png", heatmap_plot, width = 12, height = 6, bg = "white")


# function
# ---- Pairwise PERMANOVA + heatmap ----
pairwise_permanova_heatmap <- function(phy, group_var,
                                       distance = "bray",
                                       permutations = 999,
                                       value_to_plot = c("R2","p"),
                                       p_adjust = "BH",
                                       alpha = 0.05,
                                       title = NULL) {
  
  # pull groups & ensure factor with a stable order
  meta <- as(sample_data(phy), "data.frame")
  stopifnot(group_var %in% colnames(meta))
  grp <- factor(meta[[group_var]])
  levels_grp <- levels(grp)
  
  # distance on all samples, then subset per pair so indices stay aligned
  dist_all <- phyloseq::distance(phy, method = distance)
  dist_all <- as.matrix(dist_all)
  
  # storage
  R2_mat <- matrix(NA_real_, nrow = length(levels_grp), ncol = length(levels_grp),
                   dimnames = list(levels_grp, levels_grp))
  p_mat  <- R2_mat
  
  # loop pairs
  for (i in seq_along(levels_grp)) {
    for (j in seq_along(levels_grp)) {
      if (i >= j) next  # upper triangle only
      keep <- grp %in% c(levels_grp[i], levels_grp[j])
      if (sum(keep) < 4) next  # need enough samples
      
      # sub meta & dist
      sub_meta <- meta[keep, , drop = FALSE]
      sub_meta[[group_var]] <- droplevels(factor(sub_meta[[group_var]]))
      samp_ids <- rownames(sub_meta)
      sub_d <- as.dist(dist_all[samp_ids, samp_ids, drop = FALSE])
      
      # run PERMANOVA
      a2 <- adonis2(sub_d ~ sub_meta[[group_var]], permutations = permutations)
      # pick first (grouping) row
      R2_val <- a2$R2[1]
      p_val  <- a2$`Pr(>F)`[1]
      
      R2_mat[i, j] <- R2_val
      p_mat[i, j]  <- p_val
    }
  }
  
  # p-adjust across the tested pairs
  padj_vec <- p.adjust(p_mat[upper.tri(p_mat, diag = FALSE)], method = p_adjust)
  p_adj_mat <- p_mat
  p_adj_mat[upper.tri(p_adj_mat, diag = FALSE)] <- padj_vec
  
  # choose the metric to show
  show_mat <- if (value_to_plot == "R2") R2_mat else p_adj_mat
  
  # long df for plotting
  df <- as.data.frame(show_mat) |>
    mutate(row = rownames(show_mat)) |>
    pivot_longer(-row, names_to = "col", values_to = "val") |>
    mutate(row = factor(row, levels = rev(levels_grp)),
           col = factor(col, levels = levels_grp)) |>
    filter(as.integer(row) > (length(levels_grp) - as.integer(col))) # keep upper triangle
  
  # significance mask using adjusted p-values
  p_df <- as.data.frame(p_adj_mat) |>
    mutate(row = rownames(p_adj_mat)) |>
    pivot_longer(-row, names_to = "col", values_to = "p_adj") |>
    mutate(row = factor(row, levels = rev(levels_grp)),
           col = factor(col, levels = levels_grp))
  
  df <- left_join(df, p_df, by = c("row","col"))
  
  # label text
  df$label <- if (value_to_plot == "R2") {
    ifelse(is.na(df$val), "", sprintf("%.2f", df$val))
  } else {
    ifelse(is.na(df$val), "", sprintf("%.3f", df$val))
  }
  
  # choose palette & limits
  fill_lab <- if (value_to_plot == "R2") "R²" else "adj. p"
  fill_limits <- if (value_to_plot == "R2") c(0, 1) else c(0, 1)
  
  gg <- ggplot(df, aes(x = col, y = row, fill = val)) +
    geom_tile(color = NA) +
    geom_text(aes(label = label),
              size = 4, fontface = "bold", color = "black") +
    scale_fill_gradient(low = "#fde0dd", high = "#a50f15",
                        limits = fill_limits, na.value = "grey95",
                        name = fill_lab) +
    coord_equal() +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(),
          legend.position = "right") +
    ggtitle(title %||% paste0("Pairwise PERMANOVA (", distance, ") by ", group_var)) +
    guides(fill = guide_colorbar(frame.colour = "black"))
  

  gg
}
