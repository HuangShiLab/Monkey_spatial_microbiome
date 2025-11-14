#######################################
# Monkey project
#
# pathway overlap across location
#
# Author: HOU Shuwen
#######################################

library(tidyverse)
library(ComplexUpset)

# ----------------------------
# Parameters you can tweak
# ----------------------------
presence_prop_threshold <- 0.50   # >80% of samples in a location
per_sample_presence_min  <- 0     # abundance > this value counts as present in a sample
locations_order <- c("oral", "esophagus", "stomach", "small_intestine", "large_intestine")

# Set working directory
setwd("~/Downloads/monkey_data")

# Load unstratified pathway abundance table
pathway <- read.table("merged_unstratified_pathway.tsv", 
                      sep = '\t', row.names = 1, header = TRUE)

# Convert all values to numeric
pathway[] <- lapply(pathway, function(x) as.numeric(as.character(x)))


# Metadata
metadata <- read.table("kraken_metadata.txt", sep = "\t", header = TRUE)

# ----------------------------
# Align samples & locations
# ----------------------------
# Keep only samples present in both pathway (columns) and metadata$ID
common_samples <- intersect(colnames(pathway), metadata$ID)
pathway <- pathway[, common_samples, drop = FALSE]
metadata  <- metadata %>% filter(ID %in% common_samples)

# Enforce desired location order (keep all listed even if absent)
metadata$location <- factor(metadata$location, levels = locations_order)

# ----------------------------
# Presence by location (> 80%)
# ----------------------------
presence_mat <- sapply(locations_order, function(loc) {
  # Samples that claim this location, then intersect with actual columns present
  loc_samples <- metadata$ID[metadata$location == loc]
  loc_samples <- intersect(loc_samples, colnames(pathway))
  
  # If no columns remain, return zeros
  if (length(loc_samples) == 0L) {
    return(rep(0L, nrow(pathway)))
  }
  
  # Safe subset
  mat_loc <- pathway[, loc_samples, drop = FALSE]
  
  present_counts <- rowSums(mat_loc > per_sample_presence_min, na.rm = TRUE)
  denom          <- rowSums(!is.na(mat_loc))
  prop_present   <- ifelse(denom > 0, present_counts / denom, 0)
  
  as.integer(prop_present > presence_prop_threshold)
})

presence_df <- as.data.frame(presence_mat, stringsAsFactors = FALSE)
colnames(presence_df) <- locations_order
rownames(presence_df) <- rownames(pathway)

# ----------------------------
# Prep data for ComplexUpset
# ----------------------------
# Convert presence_df (pathway x locations) into a tibble with logical set columns
upset_df <- presence_df %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  # keep only locations we care about (and in correct order)
  dplyr::select(pathway, dplyr::all_of(locations_order)) %>%
  # drop pathway that are absent in all locations
  dplyr::filter(rowSums(dplyr::across(dplyr::all_of(locations_order))) > 0) %>%
  # convert 0/1 -> logical for ComplexUpset
  dplyr::mutate(dplyr::across(dplyr::all_of(locations_order), ~ .x == 1))

# ----------------------------
# Draw UpSet plot
# ----------------------------
p <- ComplexUpset::upset(
  upset_df,
  locations_order,                              # the set columns
  name = "location",
  min_size = 10,
  base_annotations = list(
    "Intersection size" = intersection_size() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  ))

print(p)

# Save final figure
ggsave("./plots/2C1.png", p, width = 6, height = 4.5, bg = "white")

# --------------------------------------------
# Make supplementary intersection table
# --------------------------------------------

# -------- helper: pathway in an exact intersection --------
get_intersection_pathways <- function(df, sets, all_sets) {
  df %>%
    filter(if_all(all_of(sets), ~ .x)) %>%                               # present in all chosen sets
    filter(if_all(all_of(setdiff(all_sets, sets)), ~ !.x)) %>%           # absent in all other sets
    pull(pathway)
}

# -------- build list of intersections -> list-column tibble --------
all_sets <- locations_order
intersection_tbl <- map_dfr(
  .x = unlist(map(1:length(all_sets), ~ combn(all_sets, .x, simplify = FALSE)), recursive = FALSE),
  .f = function(sets) {
    sp <- get_intersection_pathways(upset_df, sets, all_sets)
    tibble(
      intersection = paste(sets, collapse = "&"),
      n_pathway    = length(sp),
      pathways_list = list(sp)
    )
  }
) %>%
  # keep only non-empty intersections
  filter(n_pathway > 0)

# -------- unnest so each pathway is its own row --------
intersection_long <- intersection_tbl %>%
  unnest_longer(pathways_list, values_to = "pathway") %>%
  select(intersection, n_pathway, pathway) %>%
  arrange(desc(n_pathway), intersection, pathway)

# save
write.csv(intersection_long, "intersection_pathway_table_long.csv", row.names = FALSE)
