#######################################
# Monkey Study
#
# Procruste test, mocosal and content
#
# Author: HOU Shuwen
#######################################

# Load required packages
library(vegan)
library(ggplot2)

# Set working directory and load data
setwd("~/Downloads/monkey_mucosa")
mucosa_genus <- read.table("species_abundance.txt", header = TRUE, row.names = 1)
mucosa_meta <- read.table("meta.txt", sep = "\t", header = TRUE, row.names = 1)

setwd("~/Downloads/monkey_data")
content_genus <- read.table("2b_species_abundance.txt", header = TRUE, row.names = 1)
content_meta <- read.table("kraken_metadata.txt", sep = "\t", header = TRUE, row.names = 1)

# Ensure sample IDs are consistent and match samples
common_samples <- intersect(colnames(mucosa_genus), colnames(content_genus))
mucosa_genus <- mucosa_genus[,common_samples]
content_genus <- content_genus[,common_samples]

# Normalize
mucosa_genus <- mucosa_genus / colSums(mucosa_genus)
content_genus <- content_genus / colSums(content_genus)

# Binarize data for Jaccard distance
mucosa_bin <- (mucosa_genus > 0) * 1
content_bin <- (content_genus > 0) * 1

# Jaccard distance
mucosa_dist  <- vegdist(t(mucosa_bin),  method = "jaccard")
content_dist <- vegdist(t(content_bin), method = "jaccard")

# Compute Jaccard distances
mucosa_pcoa <- cmdscale(mucosa_dist, k = 2, eig = TRUE)
content_pcoa <- cmdscale(content_dist, k = 2, eig = TRUE)

#  After computing mucosa_pcoa and content_pcoa on the correct distance matrices:
proc <- procrustes(mucosa_pcoa$points,
                   content_pcoa$points,
                   symmetric = TRUE)

proc_test <- protest(mucosa_pcoa$points,
                     content_pcoa$points,
                     permutations = 999)

print(proc_test)    # your Procrustes m² and p-value
plot(proc, kind = 1)  # base R plot of matched points

# Extract coordinates
mucosa_coords <- as.data.frame(mucosa_pcoa$points)
content_coords <- as.data.frame(content_pcoa$points)
colnames(mucosa_coords) <- colnames(content_coords) <- c("PC1", "PC2")

# Add sample info
mucosa_coords$Sample <- rownames(mucosa_coords)
content_coords$Sample <- rownames(content_coords)
mucosa_coords$Type <- "Mucosa"
content_coords$Type <- "Content"

# extract the two sets of coordinates
mucosa_df <- data.frame(proc$X,
                        Type   = "Mucosa",
                        Sample = rownames(proc$X))
content_df <- data.frame(proc$Yrot,
                         Type   = "Content",
                         Sample = rownames(proc$Yrot))
colnames(mucosa_df)[1:2]  <- colnames(content_df)[1:2]  <- c("PC1", "PC2")

# combine
df_proc <- rbind(mucosa_df, content_df)

# now extract just the sample ID (assuming Sample is e.g. "S1_mucosa")
df_proc$ID <- sub("_[^_]+$", "", df_proc$Sample)

p <- ggplot(df_proc, aes(PC1, PC2, color = Type)) +
  geom_point(aes(shape = Type), size = 4, alpha = 0.8) +
  # draw a segment _for each_ matched pair:
  geom_segment(data = df_proc %>% 
                 group_by(ID) %>% 
                 summarize(x1 = PC1[Type=="Mucosa"],
                           y1 = PC2[Type=="Mucosa"],
                           x2 = PC1[Type=="Content"],
                           y2 = PC2[Type=="Content"]),
               aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "lightgrey", alpha = 0.5) +
  scale_color_manual(values = c("Mucosa" = "#97b74e", "Content" = "#bead3a")) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black"),
        legend.text = element_text(size = 12),
        plot.caption     = element_text(size = 12)) +
  labs(title   = "Mucosa vs Content Community Congruence (Procrustes)",
       caption = paste0("Procrustes M² = ", round(proc_test$ss, 3),
                        ", p < ", format.pval(proc_test$signif, digits = 3)))


# Save plot
ggsave("plots/3D.png", p, width = 6, height = 5)
