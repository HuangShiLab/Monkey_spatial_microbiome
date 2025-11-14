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
library(ggbreak)

# Step 1: Read the Excel file
setwd("~/Downloads/monkey_data")
df <- read_excel("qPCR_copy.xlsx")

site_order <- rev(c("oral", "esophagus", "stomach","DJ", "PJ",  "IL", "CE", "PC", "DC"))

# Calculate mean, SD, and SE per site (choose either SD or SE for error bars)
df_summary <- df %>%
  group_by(site) %>%
  summarise(
    mean_qPCR = mean(qPCR, na.rm = TRUE),
    sd_qPCR = sd(qPCR, na.rm = TRUE),  # Standard deviation
    n = n(),                           # Sample size (for SE)
    se_qPCR = sd_qPCR / sqrt(n)        # Standard error
  ) %>%
  mutate(site = factor(site, levels = site_order)) %>%
  filter(site %in% site_order)

# Step 2: Plot with error bars (use either sd_qPCR or se_qPCR)
p <- ggplot(df_summary, aes(x = site, y = mean_qPCR)) +
  geom_bar(stat = "identity", fill = "#3E4A89FF", width = 0.7) +  # Slightly narrower bars for clarity
  # Add error bars: ymin = mean - error, ymax = mean + error
  geom_errorbar(
    aes(ymin = mean_qPCR - se_qPCR, ymax = mean_qPCR + se_qPCR),  # Replace se_qPCR with sd_qPCR if preferred
    width = 0.2,  # Width of the error bar "caps"
    color = "black", 
    linewidth = 0.5
  ) +
  labs(
    title = "Microbial load\n(Square-root scale)",
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
  scale_y_sqrt() + 
  coord_flip() 

p
ggsave("./plots/2D.png", p, width = 3, height = 4, bg = "white")
