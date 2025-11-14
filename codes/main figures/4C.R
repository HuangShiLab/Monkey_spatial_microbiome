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
    values_to = "IsoalloLCA"
  )

# 5. Join metadata to bring in Location
ealca_long <- ealca_long %>%
  left_join(meta, by = c("SampleID" = "Sample ID")) %>%
  mutate(
    location = factor(.data$location,
                      levels = c("DJ","PJ","IL","CE","PC","DC")))

# 6. Other
totals_long <- concentration %>%
  filter(grepl("^Sum$", Metabolites, ignore.case = TRUE)) %>%
  select(all_of(sample_cols)) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "SampleID",
    values_to = "Total"
  )

plot_df <- ealca_long %>%
  left_join(totals_long, by = "SampleID") %>%
  mutate(
    Other = pmax(0, Total - IsoalloLCA)  # guard against negatives
  )


plot_df_mean <- plot_df %>%
  group_by(location) %>%
  summarise(
    IsoalloLCA = mean(IsoalloLCA, na.rm = TRUE),
    Other      = mean(Other,      na.rm = TRUE),
    .groups = "drop"
  )

# 6. Plot boxplots of EALCA by Location
p4 <- ggplot(plot_df_mean, aes(x = location)) +
  geom_col(aes(y = IsoalloLCA, fill = "IsoalloLCA"), position = "stack") +
  geom_col(aes(y = Other,      fill = "Other"),      position = "stack") +
  scale_fill_manual(
    name   = "Bile Acid",
    values = c("IsoalloLCA" = "#bdd048", "Other" = "grey80"),
    breaks = c("IsoalloLCA", "Other"),
    labels = c("IsoalloLCA", "Other")
  ) +
  scale_y_continuous(limits = c(0, 800000), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "IsoalloLCA and Other Concentration by Location", y = "Concentration") +
  theme(
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(size = 0.5, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 16),
    axis.title       = element_text(size = 14),
    strip.text       = element_text(face = "bold")
  )


p4

ggsave("./4D.png", p4, width = 6, height = 4, bg = "white")

