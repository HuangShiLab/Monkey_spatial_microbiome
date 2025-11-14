#######################################
# Monkey project
#
# Bile acids: Absolute abundance box plot
#
# Author: HOU Shuwen
#######################################
# Load required packages
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

# Monkey
# Read concentration data
setwd("~/Downloads/bile acids")
concentration <- read_csv("concentration.csv")

# Extract sum row and transpose
sum_data <- concentration %>% 
  filter(Metabolites == "Sum") %>%
  select(-Metabolites, -Group) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample_ID") %>%
  rename(Sum_Concentration = V1)

# Read metadata
metadata <- read_excel("meta.xlsx") %>%
  rename(Sample_ID = "Sample ID", Location = location)

# Merge data
merged_data <- sum_data %>%
  inner_join(metadata, by = "Sample_ID") %>%
  mutate(Sum_Concentration = as.numeric(Sum_Concentration),
         Location = factor(Location, levels = c("DJ", "PJ", "IL", "CE", "PC", "DC")))

# Create boxplot
p1<-ggplot(merged_data, aes(x = Location, y = Sum_Concentration, fill = Location)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 1800000)) +
  geom_beeswarm(alpha = 0.5, color = "black", size = 1, cex = 0.5) +
  labs(title = "Sum Concentration across Locations (Monkey)",
       y = "Sum concentration", x = "Location") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA), 
        axis.line = element_line(size = 0.5, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = ifelse(levels(merged_data$Location) %in% c("PJ", "DJ", "IL"), 
                                    "#f0aa73", "#eeed89"))
p1

# Human
# Read concentration data
setwd("~/Downloads/nature2023_data")
concentration <- read_csv("bile_salts.csv")

# Extract sum row and transpose
sum_data <- concentration %>% 
  filter(`Subject Number` == "sum" |`Subject Number` == "Metabolites") %>%
  t() %>%
  as.data.frame() %>%
  rename(Metabolites = V1) %>%
  rename(Sum_Concentration = V2) %>%
  slice(-1) %>%
  mutate(Sum_Concentration = as.numeric(Sum_Concentration),
         Metabolites = factor(Metabolites, levels = c("DJ", "PJ", "IL", "PC","DC")))


# Create boxplot
p2<-ggplot(sum_data, aes(x = Metabolites, y = Sum_Concentration, fill = Metabolites)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 6e+06)) +
  geom_beeswarm(alpha = 0.5, color = "black", size = 1, cex = 0.5) +
  labs(title = "Sum Concentration across Locations (Human)",
       y = "Sum concentration", x = "Location") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA), 
        axis.line = element_line(size = 0.5, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = ifelse(levels(sum_data$Metabolites) %in% c("PJ", "DJ", "IL"), 
                                    "#f0aa73", "#eeed89"))
p2
concentration_plot <- ggarrange(p1, p2)
ggsave("./4A.png", concentration_plot, width = 12, height = 4, bg = "white")

