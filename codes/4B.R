#######################################
# Monkey project
#
# Bile acids: Classification plot
#
# Author: HOU Shuwen
#######################################
# Load required packages
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Monkey
# Read and process concentration data
setwd("~/Downloads/bile acids")
concentration_data <- read_csv("concentration.csv") %>%
  filter(Metabolites != "sum") 

# Read metadata
metadata <- read_excel("meta.xlsx", sheet = "Sheet1") %>%
  rename(Sample_ID = "Sample ID", Location = location)

# 2. Tidy concentration data
#    Columns: Metabolites, Group, then sample columns
conc_long <- concentration_data %>%
  pivot_longer(
    cols = -c(Metabolites, Group),
    names_to  = "Sample",
    values_to = "Concentration"
  ) %>%
  rename(
    BileAcid = Metabolites,
    Classification = Group
  )

# 3. Join with metadata
#    Metadata columns: Sample ID and location
df <- conc_long %>%
  left_join(
    metadata %>% 
      rename(Sample = `Sample_ID`, LocationPart = Location) %>% 
      select(Sample, LocationPart),
    by = "Sample"
  ) %>%
  mutate(
    LocationPart = factor(LocationPart, levels = c("PJ","DJ","IL","CE","PC","DC"))
  ) %>%
  filter(!is.na(LocationPart))

df_avg <- df %>%
  group_by(LocationPart, Classification) %>%
  summarize(
    TotalConcentration = mean(Concentration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(LocationPart) %>%
  mutate(
    RelativeAbundance = TotalConcentration / sum(TotalConcentration)
  ) %>%
  ungroup()

# 3. Plot: one stacked bar per LocationPart, filled by Classification
p2<-ggplot(df_avg, aes(x = LocationPart, y = RelativeAbundance, fill = Classification)) +
  geom_col() +
  labs(
    title = "Average Abundance Split by\nBile Acid Composition (Monkey)",
    x     = "Location",
    y     = "Mean Relative Abundance",
    fill  = "Acid Class"
  ) +
  scale_fill_viridis_d(option = "viridis", direction = 1, end = 0.9) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA), 
    axis.line = element_line(size = 0.5, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12),
    strip.text = element_text(face = "bold")
  )
p2

# Human
# Read and process concentration data
setwd("~/Downloads/nature2023_data")
concentration_data <- read_excel("BA_concentration.xlsx", sheet = 2, skip = 1)


# 2. Tidy concentration data
#    Columns: Metabolites, Group, then sample columns
conc_long <- concentration_data %>%
  pivot_longer(
    cols = -c(Metabolite, Group),
    names_to  = "Site",
    values_to = "Concentration"
  ) %>%
  rename(
    BileAcid = Metabolite,
    Classification = Group
  ) %>%
  filter(Site %in% c("Small Intestine 1", "SI 2",
                     "Large Intestine 1", "LI 2", "Stool")) %>%
  mutate(
    Site = factor(Site, levels = c("Small Intestine 1", "SI 2",
                                   "Large Intestine 1", "LI 2", "Stool"))
  )



# 3. Plot: one stacked bar per LocationPart, filled by Classification
p3<-ggplot(conc_long, aes(x = Site, y = Concentration, fill = Classification)) +
  geom_col() +
  labs(
    title = "Average Abundance Split by\nBile Acid Composition (Human)",
    x     = "Location",
    y     = "Mean Relative Abundance",
    fill  = "Acid Class"
  ) +
  scale_fill_viridis_d(option = "viridis", direction = 1, end = 0.9) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA), 
    axis.line = element_line(size = 0.5, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12),
    strip.text = element_text(face = "bold")
  )
p3
abundance_plot <- ggarrange(p2, p3,
                                  common.legend = TRUE, legend = "right")
ggsave("./4B.png", abundance_plot, width = 10, height = 4, bg = "white")
