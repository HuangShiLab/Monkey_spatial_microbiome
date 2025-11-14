# ============================
# Enriched CAZymes & PULs (R)
# ============================

# ===== Setup =====
library(tidyverse)
library(janitor)



# ---------- Load dbCAN outputs ----------
# Expecting long format with columns: Family, Sample, Count
setwd("~/Downloads/monkey_copri")


# ----------------------------
# Parameters you can adjust
# ----------------------------
min_prev <- 0.75    # minimum prevalence
top_k    <- 30      # always keep at least top_k by combined rank
pul_mode <- "mag"        # "mag" = % MAGs with ≥1 PUL for a category; "pul" = total PUL counts

# ----------------------------
# Load inputs (robust to headers)
# ----------------------------
cazy <- readr::read_tsv("CAZy_counts_long.tsv", col_types = readr::cols()) %>%
  clean_names()  # expect columns: family, sample, count
stopifnot(all(c("family","sample","count") %in% names(cazy)))

pul_ids_raw <- readr::read_tsv("PUL_ids_per_sample.tsv", col_types = readr::cols(), col_names = TRUE)
# Rename first two cols -> sample, pul_id (handles files without headers)
pul_ids <- pul_ids_raw %>%
  transmute(sample = .[[1]], pul_id = .[[2]]) %>%
  filter(!is.na(pul_id), pul_id != "") %>%
  mutate(sample = as.character(sample), pul_id = as.character(pul_id))

# substrate mapping
substrate <- try(readr::read_tsv("PUL_substrates_long.tsv", col_types = readr::cols()) %>% clean_names(), silent = TRUE)



# ==================================================
# 1) CAZy enrichment (prevalence)
# ==================================================

# Total MAGs (samples)
n_mags <- cazy %>% distinct(sample) %>% nrow()

# Prevalence per family (Fehlner-Peach logic = presence across MAGs)
cazy_prev <- cazy %>%
  filter(count > 0) %>%
  filter(!str_starts(family, "GT")) %>%
  filter(family != "CBM0") %>%
  distinct(sample, family) %>%
  dplyr::count(family, name = "count") %>%
  mutate(prevalence = count / n_mags) %>%
  filter(prevalence >= min_prev)

# Abundance per family (Blanco-Míguez logic)
cazy_abund <- cazy %>%
  group_by(family) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),
    mean_count  = mean(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rel_total = total_count / sum(total_count)
  )


# ==================================================
# 2) PUL enrichment (prevalence)
# ==================================================

n_mags_pul <- pul_ids %>% distinct(sample) %>% nrow()

# Prevalence of each PUL across samples (presence/absence)
pul_prev <- pul_ids %>%
  distinct(sample, pul_id) %>%
  dplyr::count(pul_id, name = "counts") %>%
  mutate(prevalence = counts / n_mags_pul) %>%
  arrange(desc(counts)) 

n_mags_sub <- substrate %>% distinct(sample) %>% nrow()

# Prevalence of each substrate across samples (presence/absence)
substrate_prev <- substrate %>%
  distinct(sample, substrate) %>%
  dplyr::count(substrate, name = "counts") %>%
  mutate(
    prevalence = counts / n_mags_pul
  ) %>%
  arrange(desc(counts))

map_group <- function(x){
  case_when(
    x %in% c("GH13","GH57","GH77","GH133","CBM20","CBM48") ~ "Glycogen/Starch",
    x %in% c("CE1","CE4","GH39","GH43","CBM32")            ~ "Hemicellulose",
    x %in% c("GH28","GH105","CE8","CE12","CE20")           ~ "Pectin",
    x %in% c("GH23","GH25","CE11","CBM50")          ~ "Peptidoglycan",
    x %in% c("GH2","GH20","GH29","GH31","GH92","GH95","GH97","GH109","GH127") ~ "Oligosaccharide"
  )
}

cazy_prev <- cazy_prev %>% mutate(Functional_Group = map_chr(family, map_group))

grp_order <- c("Glycogen/Starch",
               "Hemicellulose",
               "Pectin",
               "Peptidoglycan",
               "Oligosaccharide")

cazy_prev <- cazy_prev %>%
  mutate(Functional_Group = factor(Functional_Group, levels = grp_order)) %>%
  arrange(Functional_Group, desc(prevalence)) %>%
  mutate(family = factor(family, levels = unique(family))) 

p <- ggplot(cazy_prev,
       aes(x = family,
           y = "•",
           size = prevalence,
           fill = Functional_Group)) +
  geom_point(shape = 21) +
  scale_fill_viridis_d(name = "Functional group", guide = "none") +
  labs(x = NULL, y = NULL, fill = NULL, size = NULL,
       title = "CAZyme functional categories of Monkey S. copri MAGs") +
  facet_grid(~ Functional_Group, scales = "free_x", space = "free_x") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle=90),
        panel.grid = element_blank(),
        label = NULL,
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA)) +
  scale_size(range = c(4, 8), limits = c(0.5, 1),
             breaks = seq(0.5, 1, 0.25), name = "Prevalence")

ggsave("./1E.png", p, width = 10, height = 1.8, bg = "transparent")

