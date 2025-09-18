#######################################
# Monkey project
#
# Bile acids: isoalloLCA concentration
#
# Author: HOU Shuwen
#######################################
# Load required libraries
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Read in the data
setwd("~/Downloads/bile acids")
concentration <- read_csv("concentration.csv")
meta <- read_excel("meta.xlsx")

ealca <- concentration %>%
  filter(grepl("Isoallolithocholic Acid", Metabolites, ignore.case = TRUE))

# 3. Identify only real sample columns (exclude QC)
sample_cols <- intersect(names(ealca), meta$`Sample ID`)

# 4. Pivot to long format
ealca_long <- ealca %>%
  select(all_of(sample_cols)) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "SampleID",
    values_to = "Concentration"
  )

# 5. Join metadata to bring in Location
ealca_long <- ealca_long %>%
  left_join(meta, by = c("SampleID" = "Sample ID")) %>%
  mutate(
    Location = factor(.data$location,
                      levels = c("PJ","DJ","IL","CE","PC","DC")))

# 6. Plot boxplots of EALCA by Location
p4 <- ggplot(ealca_long, aes(x = Location, y = Concentration, fill = Location)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1500000)) +
  geom_beeswarm(alpha = 0.5, color = "black", size = 1, cex = 0.5) +
  labs(
    title = "Isoallo-LCA Concentrations Across Locations",
    y     = "Isoallo-LCA concentration"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA), 
    axis.line = element_line(size = 0.5, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = ifelse(levels(ealca_long$Location) %in% c("PJ", "DJ", "IL"), 
                                    "#31688EFF", "#3E4A89FF"))

ggsave("./4D.png", p4, width = 5.5, height = 4, bg = "white")

