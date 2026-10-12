#!/usr/bin/env Rscript
# =============================================================================
# Fig. 5 - number of species detected per biological replicate, by amplicon
# (8.0 and 9.5 kb). Box plots group the three replicates at each sampling
# location; points show the individual replicates.
#   A) field sampling locations   B) aquarium (AQ)
# Colours match the Fig. 5 in the submitted manuscript.
#
# Input (written by 07_stats_figures.R):
#   results/stats/species_per_sample.csv   one row per sample x fragment
#     (species-level detections after the barcode-bleed filter)
#
# Output:
#   results/figures/fig5_species_per_sample.pdf / .png   (180 mm wide)
#
# Usage (from the lrpcr_blast folder):  Rscript 11_fig5_species_per_sample.R
#   or give the folder:                 Rscript 11_fig5_species_per_sample.R /path/to/lrpcr_blast
# Packages: dplyr, tidyr, readr, stringr, ggplot2 (tidyverse) and patchwork
#   (install.packages("patchwork") if needed).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
proj <- if (length(args) >= 1) args[1] else "."
in_file <- file.path(proj, "results", "stats", "species_per_sample.csv")
fig_dir <- file.path(proj, "results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(in_file)) stop("Not found: ", in_file, " - run 07_stats_figures.R first.")

FRAGS  <- c("8.0kb", "9.5kb")
LABELS <- c("8.0kb" = "8.0 kb", "9.5kb" = "9.5 kb")
COLS   <- c("8.0kb" = "#483D8B", "9.5kb" = "#C1CDC1")   # as in the submitted Fig. 5
DODGE  <- 0.75

d <- read_csv(in_file, show_col_types = FALSE) %>%
  filter(source %in% c("Field", "Aquarium"), fragment %in% FRAGS) %>%
  select(sample_id, source, station, replicate, fragment, n_species) %>%
  # every replicate gets a value for both amplicons (0 if none detected)
  complete(nesting(sample_id, source, station, replicate), fragment = FRAGS,
           fill = list(n_species = 0L)) %>%
  mutate(fragment = factor(fragment, levels = FRAGS),
         stn_num  = suppressWarnings(as.numeric(str_remove(station, "^P"))))

# Field locations in station-number order (P1, P3, P4, P7, P28, P38).
field_levels <- d %>% filter(source == "Field") %>% distinct(station, stn_num) %>%
  arrange(stn_num) %>% pull(station)

panel <- function(dd, x_levels, tag, show_legend) {
  dd <- dd %>% mutate(station = factor(station, levels = x_levels))
  ymax <- max(dd$n_species)
  ggplot(dd, aes(x = station, y = n_species, fill = fragment, colour = fragment)) +
    geom_boxplot(position = position_dodge(width = DODGE), width = 0.6,
                 alpha = 0.45, outlier.shape = NA, linewidth = 0.5) +
    geom_point(position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0,
                                               dodge.width = DODGE, seed = 1),
               size = 1.8) +
    scale_fill_manual(values = COLS, labels = LABELS, name = "Fragment", drop = FALSE) +
    scale_colour_manual(values = COLS, labels = LABELS, name = "Fragment", drop = FALSE) +
    scale_y_continuous(limits = c(0, ymax + 1), breaks = scales::breaks_pretty(n = 5),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = NULL, y = NULL, tag = tag) +
    theme_classic(base_size = 10) +
    theme(axis.text = element_text(colour = "grey30"),
          plot.tag = element_text(face = "bold", size = 13),
          legend.position = if (show_legend) "right" else "none")
}

pA <- panel(filter(d, source == "Field"), field_levels, "A)", FALSE) +
  labs(y = "Species (#)", x = "Sampling location")
pB <- panel(filter(d, source == "Aquarium") %>% mutate(station = "AQ"), "AQ", "B)", TRUE)

p <- pA + pB + plot_layout(widths = c(6, 1.3), guides = "collect")

ggsave(file.path(fig_dir, "fig5_species_per_sample.pdf"), p, width = 180, height = 90, units = "mm")
ggsave(file.path(fig_dir, "fig5_species_per_sample.png"), p, width = 180, height = 90, units = "mm", dpi = 300)

# Per-group summary (should match species_per_sample_summary.csv from step 7).
se <- function(x) if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
print(d %>% group_by(source, fragment) %>%
        summarise(n = n(), mean = mean(n_species), se = se(n_species),
                  min = min(n_species), max = max(n_species), .groups = "drop"))
message("Wrote ", file.path(fig_dir, "fig5_species_per_sample.(pdf|png)"))
