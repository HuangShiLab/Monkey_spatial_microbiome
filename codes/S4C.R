#######################################
# Monkey project
#
# Bile acids: Relative abundance (small vs large intestine)
#
# Author: HOU Shuwen
#######################################
# Load required packages
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(forcats)

# Read concentration data
setwd("~/Downloads/bile acids")
concentration <- read_csv("concentration.csv") %>%
  filter(Metabolites != "sum") %>%  # Remove sum row
  select(-Group)  # Remove original group column

# 2. Load and process metadata
metadata <- read_excel("meta.xlsx", sheet = "Sheet1") %>%
  rename(Sample_ID = `Sample ID`, Location = location) %>%
  mutate(
    Region = case_when(
      Location %in% c("PJ", "DJ", "IL") ~ "Small Intestine",
      Location %in% c("CE", "PC", "DC") ~ "Large intestine",
      TRUE ~ "Other"
    )
  )

# 3. Transform concentration to long format and compute relative abundances
long_data <- concentration %>%
  pivot_longer(-Metabolites, names_to = "Sample_ID", values_to = "Concentration") %>%
  inner_join(metadata, by = "Sample_ID") %>%
  group_by(Sample_ID) %>%
  mutate(
    Sample_Total = sum(Concentration, na.rm = TRUE),
    Relative_Abundance = Concentration / Sample_Total
  ) %>%
  ungroup()

# 1. find significant metabolites (same as before)
significant_metabolites <- long_data %>%
  group_by(Metabolites) %>%
  summarise(Max_Abs = max(Relative_Abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(Max_Abs >= 0.05) %>%
  pull(Metabolites)

# 2A. keep only the significant ones
sig_df <- long_data %>%
  filter(Metabolites %in% significant_metabolites)

# 2B. everything else (those that were “<5% everywhere”)
tiny_df <- long_data %>%
  filter(!(Metabolites %in% significant_metabolites))

# 3. for each sample, collapse all “tiny” metabolites into a single “Others” row
others_df <- tiny_df %>%
  group_by(Sample_ID, Location, Region) %>% 
  summarise(
    Concentration       = sum(Concentration, na.rm = TRUE),        # (optional if you ever need raw concentration)
    Sample_Total        = unique(Sample_Total),                     # each sample’s total is the same
    Relative_Abundance  = sum(Relative_Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Metabolites = "Others") %>%
  # If you need to preserve any other columns (e.g. Region), keep them in summarise() above.
  select(Metabolites, everything())

# 4. bind the significant + “Others”
filtered_all <- bind_rows(
  sig_df,
  others_df
) %>%
  # now set the desired order for Location and Sample_ID as before:
  mutate(
    Location  = factor(Location, levels = c("PJ","DJ","IL","CE","PC","DC")),
    Sample_ID = factor(Sample_ID, levels = unique(metadata$Sample_ID)),
    # ensure that “Others” appears last in the fill legend
    Metabolites = fct_relevel(Metabolites, c(sort(setdiff(unique(significant_metabolites), "Others")), "Others"))
  )

# 5. plot exactly as before (no further changes needed)
p_all <- ggplot(filtered_all, aes(x = Sample_ID, y = Relative_Abundance, fill = Metabolites)) +
  geom_col(position = "stack") +
  facet_wrap(~ Location, scales = "free_x", ncol = 3) +
  labs(
    title = "Major Metabolites in Small & Large Intestine",
    x = "Sample ID", y = NULL
  ) +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    strip.background= element_rect(fill = "grey90", color = "grey50"),
    strip.text      = element_text(color = "black"),
    legend.text     = element_text(size = 8),
    legend.position = "right",
    panel.spacing   = unit(0.2, "lines"),
    plot.title      = element_text(hjust = .5)
  )


# print or save
print(p_all)
ggsave("4C.png", p_all, width=12, height=8)


