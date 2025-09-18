#######################################
# Monkey project
#
# Microbial load
#
# Author: HOU Shuwen
#######################################

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# Step 1: Read the Excel file
# Set working directory for mucosa data
setwd("~/Downloads/monkey_data")
df <- read_excel("qPCR_copy.xlsx")

site_order <- rev(c("oral", "esophagus", "stomach","DJ", "PJ",  "IL", "CE", "PC", "DC"))

# Calculate mean qPCR per site
df_summary <- df %>%
  group_by(site) %>%
  summarise(mean_qPCR = mean(qPCR, na.rm = TRUE)) %>%
  mutate(site = factor(site, levels = site_order)) %>%
  filter(site %in% site_order)

p <- ggplot(df_summary, aes(x = site, y = mean_qPCR)) +
  geom_bar(stat = "identity",  fill = "#3E4A89FF") +
  scale_y_log10() + 
  labs(
    title = "Microbial load (log scale)",
    x = "",
    y = "16S copy number/ g sample"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.text = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_line(color = "black", size = 0.5), 
    legend.position = "none"
  ) +
  coord_flip() 
p
ggsave("./plots/2D.png", p, width = 3, height = 4, bg = "white")

