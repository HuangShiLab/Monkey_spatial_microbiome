#######################################
# Monkey project
#
# Supplementary: CheckM QC
#
# Author: HOU Shuwen
#######################################

# ---- Packages (install once if needed) ----
# install.packages(c("tidyverse", "readxl", "cowplot"))
library(tidyverse)
library(readxl)
library(cowplot)

# ---- Read & prep ----
path <- "~/Downloads/monkey_data/MAG/CheckM.xlsx" 
# Make sure column names match your Excel ("Completeness", "Contamination")
# If they have spaces or different case, adjust here:
df <- read_excel(path)
names(df)

# Example: if your Excel columns are "Completeness" and "Contamination"
df <- df %>%
  mutate(
    Quality = case_when(
      Completeness >= 95  ~ ">95%",
      Completeness >= 90  ~ ">90%",
      Completeness >= 80  ~ ">80%",
      TRUE                ~ "<80%"
    ),
    Quality = factor(Quality, levels = c(">95%", ">90%", ">80%", "<80%"))
  )



# ---- Main scatter ----
p_scatter <-
  ggplot(df, aes(Contamination, Completeness, fill = Quality)) +
  geom_point(shape = 22, size = 2.4, alpha = 0.9) +
  labs(x = "Contamination (%)", y = "Completeness (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

# ---- Marginal histograms ----
p_top <-
  ggplot(df, aes(Contamination)) +
  geom_histogram(binwidth = 0.1, fill = "grey80", color = "grey30") +
  labs(y = "Count", x = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_right <-
  ggplot(df, aes(Completeness)) +
  geom_histogram(binwidth = 0.1, fill = "grey80", color = "grey30") +
  labs(x = NULL, y = "Count") +
  coord_flip() +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_empty <- ggplot() + theme_void()

# ---- Assemble ----
final_plot <- plot_grid(
  plot_grid(p_top, p_empty, ncol = 2, rel_widths = c(4, 1.2)),
  plot_grid(p_scatter, p_right, ncol = 2, rel_widths = c(4, 1.2), align = "h"),
  nrow = 2, rel_heights = c(1.2, 4), align = "v"
)

# ---- Save ----
ggsave("s1.pdf", final_plot, width = 8, height = 6)
