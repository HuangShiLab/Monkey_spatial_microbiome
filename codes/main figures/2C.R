#######################################
# Monkey project
#
# Species overlap across location
#
# Author: HOU Shuwen
#######################################
setwd("~/Downloads/monkey_data")

# Load packages
# install.packages(c("tidyverse"))
# install.packages("ComplexUpset")
# BiocManager::install("phyloseq")   # if not installed
# BiocManager::install("vegan")      # for estimate_richness
library(tidyverse)
library(ComplexUpset)
library(phyloseq)

# ----------------------------
# Parameters you can tweak
# ----------------------------
presence_prop_threshold <- 0.50   # >80% of samples in a location
per_sample_presence_min  <- 0     # abundance > this value counts as present in a sample
locations_order <- c("oral", "esophagus", "stomach", "small_intestine", "large_intestine")


# ============================
# Read abundance & metadata
# ============================
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

# Load sample metadata
metadata <- read.table("kraken_metadata.txt", sep = "\t", header = TRUE)



# ----------------------------
# Align samples & locations
# ----------------------------
# Keep only samples present in both bc_kraken (columns) and metadata$ID
common_samples <- intersect(colnames(bc_kraken), metadata$ID)
bc_kraken <- bc_kraken[, common_samples, drop = FALSE]
metadata  <- metadata %>% filter(ID %in% common_samples)

# Enforce desired location order (keep all listed even if absent)
metadata$location <- factor(metadata$location, levels = locations_order)

# ----------------------------
# Presence by location (> 80%)
# ----------------------------
presence_mat <- sapply(locations_order, function(loc) {
  # Samples that claim this location, then intersect with actual columns present
  loc_samples <- metadata$ID[metadata$location == loc]
  loc_samples <- intersect(loc_samples, colnames(bc_kraken))
  
  # If no columns remain, return zeros
  if (length(loc_samples) == 0L) {
    return(rep(0L, nrow(bc_kraken)))
  }
  
  # Safe subset
  mat_loc <- bc_kraken[, loc_samples, drop = FALSE]
  
  present_counts <- rowSums(mat_loc > per_sample_presence_min, na.rm = TRUE)
  denom          <- rowSums(!is.na(mat_loc))
  prop_present   <- ifelse(denom > 0, present_counts / denom, 0)
  
  as.integer(prop_present > presence_prop_threshold)
})

presence_df <- as.data.frame(presence_mat, stringsAsFactors = FALSE)
colnames(presence_df) <- locations_order
rownames(presence_df) <- rownames(bc_kraken)

# ----------------------------
# Prep data for ComplexUpset
# ----------------------------
# Convert presence_df (species x locations) into a tibble with logical set columns
upset_df <- presence_df %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species") %>%
  # keep only locations we care about (and in correct order)
  dplyr::select(species, dplyr::all_of(locations_order)) %>%
  # drop species that are absent in all locations
  dplyr::filter(rowSums(dplyr::across(dplyr::all_of(locations_order))) > 0) %>%
  # convert 0/1 -> logical for ComplexUpset
  dplyr::mutate(dplyr::across(dplyr::all_of(locations_order), ~ .x == 1))

# ----------------------------
# Draw UpSet plot
# ----------------------------
p <- ComplexUpset::upset(
  upset_df,
  locations_order,                              # the set columns
  name = "location",                            # label for the set dimension
  min_size = 10,
  base_annotations = list(
    "Intersection size" = intersection_size() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  ))

print(p)

# Save final figure
ggsave("./plots/2B1.png", p, width = 6, height = 4, bg = "white")

# --------------------------------------------
# Make supplementary intersection table
# --------------------------------------------

# -------- helper: species in an exact intersection --------
get_intersection_species <- function(df, sets, all_sets) {
  df %>%
    filter(if_all(all_of(sets), ~ .x)) %>%                               # present in all chosen sets
    filter(if_all(all_of(setdiff(all_sets, sets)), ~ !.x)) %>%           # absent in all other sets
    pull(species)
}

# -------- build list of intersections -> list-column tibble --------
all_sets <- locations_order
intersection_tbl <- map_dfr(
  .x = unlist(map(1:length(all_sets), ~ combn(all_sets, .x, simplify = FALSE)), recursive = FALSE),
  .f = function(sets) {
    sp <- get_intersection_species(upset_df, sets, all_sets)
    tibble(
      intersection = paste(sets, collapse = "&"),
      n_species    = length(sp),
      species_list = list(sp)
    )
  }
) %>%
  # keep only non-empty intersections
  filter(n_species > 0)

# -------- (optional) match your plot's min_size = 10 --------
intersection_tbl <- intersection_tbl %>%
  filter(n_species >= 10)

# -------- unnest so each species is its own row --------
intersection_long <- intersection_tbl %>%
  unnest_longer(species_list, values_to = "species") %>%
  select(intersection, n_species, species) %>%
  arrange(desc(n_species), intersection, species)

# save
write.csv(intersection_long, "intersection_species_table_long.csv", row.names = FALSE)

# peek
print(head(intersection_long, 20))