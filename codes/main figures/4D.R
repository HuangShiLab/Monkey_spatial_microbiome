#######################################
# Monkey project
#
# Gene cluster plot
#
# Author: HOU Shuwen
#######################################

# ---- packages ----
library(tidyverse)
library(gggenes)

setwd("~/Downloads/monkey_data")

# ============================================================
# 1) INTERPRETED tidy data from your annotations (no parsing)
# ============================================================
genes <- tribble(
  ~molecule,                                   ~contig,        ~gene_id,            ~gene_name,   ~role,       ~start, ~end,   ~strand,
  # species A → Parabacteroides sp900549585
  "O1_LI_10cm_L_bin_18", "k141_428973", "k141_428973_3",      "5AR",        "5AR",        2192,   2986,    "+",
  "O1_LI_10cm_L_bin_18", "k141_428973", "k141_428973_4",      "5BR",        "5BR",        2962,   4224,    "+",
  "O1_LI_10cm_L_bin_18", "k141_428973", "k141_428973_5",      "3β-HSDH",    "3β-HSDH",    4221,   5066,    "+",
  
  "O1_LI_19cm_L_bin_144", "k141_30244", "k141_30244_4",      "5AR",        "5AR",        3143,   3937,    "-",
  "O1_LI_19cm_L_bin_144", "k141_30244", "k141_30244_3",      "5BR",        "5BR",        1905,   3167,    "-",
  "O1_LI_19cm_L_bin_144", "k141_30244", "k141_30244_2",      "3β-HSDH",    "3β-HSDH",    1063,   1908,    "-",
  
  "O1_LI_28cm_L_bin_42", "k141_403932", "k141_403932_8",      "5AR",        "5AR",        9023,   9817,    "+",
  "O1_LI_28cm_L_bin_42", "k141_403932", "k141_403932_9",      "5BR",        "5BR",        9793,   11055,    "+",
  "O1_LI_28cm_L_bin_42", "k141_403932", "k141_403932_10",      "3β-HSDH",    "3β-HSDH",    11052,   11897,    "+",

  "O1_LI_46cm_L_bin_112", "k141_9009", "k141_9009_4",      "5AR",        "5AR",        3009, 3803,    "-",
  "O1_LI_46cm_L_bin_112", "k141_9009", "k141_9009_3",      "5BR",        "5BR",        1771, 3033,    "-",
  "O1_LI_46cm_L_bin_112", "k141_9009", "k141_9009_2",      "3β-HSDH",    "3β-HSDH",    929, 1774,    "-",

  "O1_LI_55cm_L_bin_93", "k141_525680", "k141_525680_3",      "5AR",        "5AR",        1052, 1846,    "+",
  "O1_LI_55cm_L_bin_93", "k141_525680", "k141_525680_4",      "5BR",        "5BR",        1822, 3084,    "+",
  "O1_LI_55cm_L_bin_93", "k141_525680", "k141_525680_5",      "3β-HSDH",    "3β-HSDH",    3081, 3926,    "+",
  
  "O2_LI_10cm_L_bin_9", "k141_10624", "k141_10624_7",      "5AR",        "5AR",        6433, 7227,    "-",
  "O2_LI_10cm_L_bin_9", "k141_10624", "k141_10624_6",      "5BR",        "5BR",        5195, 6457,    "-",
  "O2_LI_10cm_L_bin_9", "k141_10624", "k141_10624_5",      "3β-HSDH",    "3β-HSDH",    4353, 5198,    "-",
  
  # species B → Egerieousia sp900544815
  "Y2_LI_19cm_L_bin_41",     "k141_106483", "k141_106483_11",      "5AR",        "5AR",        11596, 12327,    "-",
  "Y2_LI_19cm_L_bin_41",     "k141_106483", "k141_106483_10",      "5BR",        "5BR",        10348, 11571,    "-",
  "Y2_LI_19cm_L_bin_41",     "k141_106483", "k141_106483_9",      "other",        "other",     9300, 10124,    "-",
  "Y2_LI_19cm_L_bin_41",     "k141_106483", "k141_106483_8",     "3β-HSDH",    "3β-HSDH",    8288, 9193,   "-",

  "Y2_LI_28cm_L_bin_48",     "k141_329263", "k141_329263_7",      "5AR",        "5AR",        6333,   7064,    "+",
  "Y2_LI_28cm_L_bin_48",     "k141_329263", "k141_329263_8",      "5BR",        "5BR",        7089,   8312,    "+",
  "Y2_LI_28cm_L_bin_48",     "k141_329263", "k141_329263_9",      "other",        "other",     8536,   9360,    "+",
  "Y2_LI_28cm_L_bin_48",     "k141_329263", "k141_329263_10",     "3β-HSDH",    "3β-HSDH",    9467,   10372,   "+"
  ) %>%
  mutate(
    start = pmin(start, end),
    end   = pmax(start, end),
    forward = strand == "+"
  )

# --- 1.5) Flip all '-' genes to '+' coordinates using a big per-molecule base ---
# Pick a per-molecule BIG value just larger than any coordinate on that molecule
flip_base <- genes %>%
  group_by(molecule) %>%
  summarise(
    BIG = 10^(ceiling(log10(max(c(start, end)) + 1))),  # next power of 10
    .groups = "drop"
  )

genes <- genes %>%
  left_join(flip_base, by = "molecule") %>%
  mutate(
    # normalise inputs
    s0 = pmin(start, end),
    e0 = pmax(start, end),
    
    # flip only the '-' strand: new_start = BIG - old_end; new_end = BIG - old_start
    start = if_else(strand == "-", BIG - e0, s0),
    end   = if_else(strand == "-", BIG - s0, e0),
    
    # after flipping, treat everything as '+' so arrows point right
    strand  = "+",
    forward = TRUE
  ) %>%
  select(-s0, -e0, -BIG)

# ============================================================
# 2) ALIGN ON ANCHOR (5AR) & set window (optional)
# ============================================================
anchor_role <- "5AR"
# If you later add neighbors, set a wider window (e.g., 20000 for ±10 kb)
window_bp <- NA_integer_  # keep all genes (no trimming); set to a number to clip

anchors <- genes %>%
  filter(role == anchor_role) %>%
  group_by(molecule) %>%
  summarize(anchor = start, .groups = "drop")

genes_rel <- genes %>%
  left_join(anchors, by = "molecule") %>%
  mutate(
    xmin = start - anchor,
    xmax = end   - anchor,
    mid  = (xmin + xmax)/2
  )

if (!is.na(window_bp)) {
  genes_rel <- genes_rel %>%
    filter(xmax >= -window_bp, xmin <= window_bp)
}

# Order panels (Parabacteroides first)
mol_levels <- genes_rel %>%
  distinct(molecule) %>%
  arrange(desc(str_detect(molecule, regex("^Parabacteroides", ignore_case = TRUE)))) %>%
  pull(molecule)
genes_rel$molecule <- factor(genes_rel$molecule, levels = mol_levels)

# ============================================================
# 3) PLOT
# ============================================================
role_cols <- c("5AR"="#4e5d88","5BR"="#83659d","3β-HSDH"="#d66b93","other"="grey80")
genes_rel$role <- factor(genes_rel$role, levels = c("5AR","5BR","3β-HSDH","other"))


p <- ggplot(
  genes_rel,
  aes(xmin = xmin, xmax = xmax, y = molecule, forward = forward, fill = role)
) +
  geom_gene_arrow(size = 0.25, color = "grey20",
                  arrowhead_height = grid::unit(3, "mm"),
                  arrowhead_width  = grid::unit(3, "mm")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35") +
  geom_gene_label(aes(x = mid, label = gene_name),
                  fill = "white",
                  size = 3,
                  label.padding = grid::unit(0.15, "lines"),
                  label.size = 0.15,
                  label.r = grid::unit(0.05, "lines")) +
  scale_fill_manual(values = role_cols, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 4039),
    breaks = scales::pretty_breaks(4),
    labels = ~ scales::number(.x/1000, accuracy = 0.1, suffix = " kb")
  ) +
  labs(
    x = sprintf("Position relative to %s (kb)", anchor_role),
    y = NULL, fill = "Function",
    title = "5AR and 3β-HSDH Gene Clusters in Monkey MAGs"
  ) +
  gggenes::theme_genes() +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

print(p)

# Export (vector + print)
ggsave("./plots/4D.png", p, width = 6, height = 3.5, bg = "white")

