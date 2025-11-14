#######################################
# Monkey Study
#
# GTDB-tk phylogenetic tree
#
# Author: HOU Shuwen
#######################################

# load libraries
library(metagMisc)
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)
library(viridis)
library(openxlsx)
library(ggplot2)
library("ape")
library(cowplot)
library(tidyr)
library(dtplyr)

# data import
setwd("~/Downloads/monkey_data")

# Load the GTDB-tk tree
tree <- ape::read.tree("MAG/gtdbtk.backbone.bac120.classify.tree")

# Load the GTDB-tk classification
taxonomy <- read.table("MAG/gtdbtk.bac120.summary.tsv", sep = "\t", header = TRUE)

# Extract genus information
taxonomy <- taxonomy %>%
  mutate(genus = sub(".*;g__", "", classification),  # Extract genus
         genus = sub(";s__.*", "", genus))  %>%
  mutate(phylum = sub(".*;p__", "", classification),  # Extract genus
         phylum = sub(";c__.*", "", phylum))       

# Merge metadata with the tree
tree_tips <- as.data.frame(tree$tip.label, col.names = "tip")
taxonomy <- taxonomy %>%
  filter(sample_ID %in% tree$tip.label)

# Map each genus to a representative user genome
genus_rep <- taxonomy %>%
  group_by(genus) %>%
  summarise(representative_tip = first(sample_ID))

# manipulate the metadata table
genus_meta <- taxonomy %>%
  group_by(genus) %>%
  mutate(sample = sub("cm_L.contigs.fa.metabat-bins_bin.*", "", sample_ID)) %>%
  select(phylum,genus, sample) %>%
  separate(sample, into = c("age", "location", "sub"), sep = "_") %>%
  select(-sub)

# sub-location (discard)
# mutate(sub = case_when(
#  sub < 15 ~ "CE",             
#  sub >= 16 & sub <= 30 ~ "PC", 
#  sub > 31 ~ "DC")) %>%
#  mutate(location = if_else(location == "LI", sub, location)) %>%

# Prune tree to keep only the representative genomes
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, genus_rep$representative_tip))

# Annotate genera on the tree
pruned_tree$tip.label <- genus_rep$genus[match(pruned_tree$tip.label, genus_rep$representative_tip)]

# ggtree skeleton
phylum_counts <- genus_meta %>%
  group_by(phylum) %>%
  summarise(count = n(), .groups = "drop")

find_lca <- function(node1, node2, parent_df) {
  # Initialize lists to store ancestors of node1 and node2
  ancestors_node1 <- list(node1)
  ancestors_node2 <- list(node2)
  
  # Get the parent of node1 and node2
  parent1 <- parent_df$parent[parent_df$node == node1]
  parent2 <- parent_df$parent[parent_df$node == node2]
  
  # Stop the function if any node does not have a parent
  if (length(parent1) == 0 || length(parent2) == 0) {
    return(NULL)
  }
  
  # Loop to find ancestors of node1
  while (TRUE) {
    # Find the parent of node1
    parent1 <- parent_df$parent[parent_df$node == node1]
    
    # If parent is the same as node, break the loop
    if (length(parent1) == 0 || is.na(parent1) || parent1 == node1) {
      break
    }
    
    # Add parent to ancestors
    ancestors_node1 <- c(ancestors_node1, parent1)
    
    # Update node1 to its parent for the next iteration
    node1 <- parent1
  }
  
  # Loop to find ancestors of node2
  while (TRUE) {
    # Find the parent of node2
    parent2 <- parent_df$parent[parent_df$node == node2]
    
    # If parent is the same as node, break the loop
    if (length(parent2) == 0 || is.na(parent2) || parent2 == node2) {
      break
    }
    
    # Add parent to ancestors
    ancestors_node2 <- c(ancestors_node2, parent2)
    
    # Update node2 to its parent for the next iteration
    node2 <- parent2
  }
  
  # Find the least common ancestor (LCA) by finding the first common ancestor
  lca_node <- NULL
  for (ancestor in ancestors_node1) {
    if (ancestor %in% ancestors_node2) {
      lca_node <- ancestor
      break
    }
  }
  
  return(lca_node)
}


class <- data.frame(id=c(426, 394,225,401,
                         254,386,262,233), 
                    type=c("Pseudomonadota", "Verrucomicrobiota",
                           "Spirochaetota", "Bacteroidia",
                           "Campylobacterota", "Cyanobacteriota",
                           "Bacillota", "Bacillota_I"))
class_color <- c("Pseudomonadota" = "#eacc76", 
                 "Verrucomicrobiota" = "#5c7272",
                 "Spirochaetota" = "#718c70", 
                 "Bacteroidia" = "#acab4b",
                 "Campylobacterota" = "#7c99bc", 
                 "Cyanobacteriota" = "#a5b3c1",
                 "Bacillota" = "#c7a8a3", 
                 "Bacillota_I" = "#e9b962")
gcc_plot = ggtree(pruned_tree, layout="fan", open.angle = 2, size = 0.15, alpha = 1) +
  geom_tiplab(linesize=.15, mapping = aes(label = NA)) +
  geom_highlight(data=class,mapping=aes(node=id, fill=type),alpha=0.5) +
  scale_fill_manual(values = class_color, guide = "none") + 
  coord_polar(theta = 'y', start = 0, direction = -1) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0,0,0,0))
gcc_plot

gcc_plot_data <- gcc_plot$data %>% filter(isTip)  %>%
  left_join(genus_meta, by = c("label" = "genus")) %>%
  distinct()

## layer 1: location
gcc_plot_data$type_numeric <- as.numeric(
  factor(gcc_plot_data$location, levels = c("S", "SI", "LI")))

color_type <- c("S" = "#d47090", "SI" = "#f5aa68","LI" = "#f7f06d")


layer1 <- ggplot() +
  geom_tile(
    data = gcc_plot_data,
    aes(x = y, y = type_numeric, fill = location),
    color = NA,            # no border to avoid darkening
    linewidth = 0          # (or size=0 on older ggplot2)
  ) +
  scale_fill_manual(values = color_type) +
  theme_minimal() +
  theme(
    rect = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),  # use white for consistent appearance
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(size = 0.15),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_blank()
  ) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, max(gcc_plot_data$y, na.rm = TRUE) + 1) +
  scale_y_continuous(breaks = c(1, 2), limits = c(-13, 6))

layer1

# ## layer2: age
# gcc_plot_data$type_numeric <- as.numeric(
#   factor(gcc_plot_data$age, levels = c("young", "old")))
# 
# color_type <- c("young"="#C49FA7", "old"="#6B4C54")
# 
# layer2 <- ggplot() +
#   geom_point(data = gcc_plot_data, 
#             aes(x = y, y = type_numeric, fill = age, color = age), size = 0.8) +
#   scale_fill_manual(values = color_type) +
#   scale_color_manual(values = color_type) +
#   theme_minimal() +
#   theme(rect = element_blank(),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         panel.grid = element_blank(),
#         panel.grid.major.y = element_line(size = 0.15),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         legend.position = "none",
#         plot.margin = margin(0, 0, 0, 0),
#         panel.border = element_blank()) +
#   coord_polar(start = 0, direction = -1) +
#   xlim(0, max(gcc_plot_data$y, na.rm = TRUE) + 1) +
#   scale_y_continuous(breaks = c(1,2), limits = c(-13,6))
# layer2

# combine plots
combined_plot <- ggdraw() +
#  draw_plot(layer2, x = 0, y = 0, width = 1, height = 1) + # Adjust x, y, width, height as needed
  draw_plot(layer1, x = 0, y = 0, width = 1, height = 1) + 
  draw_plot(gcc_plot, x = 0.12, y = 0.12, width = 0.76, height = 0.76) # Add the main plot
combined_plot
ggsave("./plots/tree.png", combined_plot, width = 5, height = 5)

## annoation
# Annotation plot for Location
location_annotation <- ggplot(data.frame(location = c("Stomach", "Small Intestine", "Large Intestine"), 
                                         y = 1:3), 
                              aes(x = 1, y = y, fill = location)) +
  geom_tile(width = 0.5, height = 0.5) +
  scale_fill_manual(values = c("Stomach" = "#d47090", 
                               "Small Intestine" = "#f5aa68",
                               "Large Intestine" = "#f7f06d")) +
  labs(title = "Location", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(breaks = 1:3, 
                     labels = c("Stomach", "Small Intestine", "Large Intestine"),
                     , limits = c(-1, 4)) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)))
location_annotation

# Annotation plot for Age
age_annotation <- ggplot(data.frame(age = c("Young", "Old"), 
                                    y = c(1, 2)), 
                         aes(x = 1, y = y, color = age)) +
  geom_point(size = 10) +
  scale_color_manual(values = c("Young" = "#C49FA7", "Old" = "#6B4C54")) +
  labs(title = "Age", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(breaks = c(1, 2), labels = c("Young", "Old"), limits = c(0, 3)) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)))
age_annotation
