#!/usr/bin/env Rscript
# =============================================================================
# Step 7 - detection statistics and figures for the revision.
#
# Inputs (from lrpcr_assign.py and lists/):
#   results/assignment/taxon_counts.tsv     reads per sample x fragment x taxon
#   results/assignment/sample_summary.tsv   reads per sample x fragment
#   lists/sample_metadata.csv               barcode -> sample, station, source
#   lists/pcr_success_tableS1.csv           gel bands for all 72 field samples (Table S1)
#   lists/aquarium_species.tsv              tank census (Table S2)
#
# Outputs:
#   results/stats/*.csv, results/stats/stats_summary.txt
#   results/figures/fig6_heatmap_field.(pdf|png)      Fig. 6: reads per species x field sample
#   results/figures/fig7_heatmap_aquarium.(pdf|png)   Fig. 7: reads per species x aquarium sample
#
# Reviewer comments addressed: R3.1, R3.7, R3.10, R3.14, R3.15, R2.7, R4.2, R4.27.
#
# Usage (from the lrpcr_blast folder):  Rscript 07_stats_figures.R
#   or give the folder:                 Rscript 07_stats_figures.R /path/to/lrpcr_blast
# Packages: dplyr, tidyr, readr, stringr, ggplot2 (all in the tidyverse).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
proj <- if (length(args) >= 1) args[1] else "."
res_dir <- file.path(proj, "results", "assignment")
out_dir <- file.path(proj, "results", "stats")
fig_dir <- file.path(proj, "results", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

FRAGS <- c("8.0kb", "9.5kb")
se <- function(x) if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
log_lines <- character(0)
say <- function(...) { line <- paste0(...); log_lines <<- c(log_lines, line); message(line) }

# ---------------------------------------------------------------- inputs -----
meta <- read_csv(file.path(proj, "lists", "sample_metadata.csv"), show_col_types = FALSE) %>%
  transmute(barcode   = demulti,
            sample_id = SampleID,
            source    = Source,                                   # Field / Aquarium / NC
            station   = case_when(Source == "Field" ~ paste0("P", Station_num),
                                  Source == "Aquarium" ~ "AQ",
                                  TRUE ~ as.character(Station_num)),
            station_name = str_trim(Station_name),
            replicate = str_extract(SampleID, "[ABC]$"))

counts <- read_tsv(file.path(res_dir, "taxon_counts.tsv"), show_col_types = FALSE,
                   col_types = cols(kept = col_character(), .default = col_guess())) %>%
  mutate(kept = tolower(kept) == "true",
         barcode = str_extract(sample, "barcode\\d+"))

summ <- read_tsv(file.path(res_dir, "sample_summary.tsv"), show_col_types = FALSE) %>%
  mutate(barcode = str_extract(sample, "barcode\\d+"))

pcr  <- read_csv(file.path(proj, "lists", "pcr_success_tableS1.csv"), show_col_types = FALSE)
tank <- read_tsv(file.path(proj, "lists", "aquarium_species.tsv"), show_col_types = FALSE)

sp <- counts %>% filter(kept, rank == "species") %>% inner_join(meta, by = "barcode")
samples <- meta %>% filter(source %in% c("Field", "Aquarium"))

# ------------------------------------------- 1. species per sample (R3.14) ---
per_frag <- sp %>%
  group_by(barcode, fragment) %>%
  summarise(n_species = n_distinct(taxon), species_reads = sum(reads), .groups = "drop")
per_any <- sp %>%
  group_by(barcode) %>%
  summarise(n_species = n_distinct(taxon), species_reads = sum(reads), .groups = "drop") %>%
  mutate(fragment = "either")

assigned_reads <- summ %>%
  filter(fragment %in% FRAGS) %>%
  group_by(barcode, fragment) %>% summarise(reads_after_q = sum(reads), .groups = "drop") %>%
  bind_rows(summ %>% filter(fragment %in% FRAGS) %>% group_by(barcode) %>%
              summarise(reads_after_q = sum(reads), .groups = "drop") %>% mutate(fragment = "either"))

per_sample <- samples %>%
  crossing(fragment = c(FRAGS, "either")) %>%
  left_join(bind_rows(per_frag, per_any), by = c("barcode", "fragment")) %>%
  left_join(assigned_reads, by = c("barcode", "fragment")) %>%
  mutate(across(c(n_species, species_reads, reads_after_q), ~ replace_na(.x, 0))) %>%
  arrange(source, station, replicate, fragment)
write_csv(per_sample, file.path(out_dir, "species_per_sample.csv"))

by_group <- per_sample %>%
  group_by(source, fragment) %>%
  summarise(n_samples = n(), mean_species = mean(n_species), se_species = se(n_species),
            min_species = min(n_species), max_species = max(n_species),
            samples_with_zero = sum(n_species == 0), .groups = "drop")
write_csv(by_group, file.path(out_dir, "species_per_sample_summary.csv"))
say("Species per sample (mean +/- SE, zeros included):")
for (i in seq_len(nrow(by_group))) with(by_group[i, ], say(sprintf(
  "  %-8s %-6s n=%d  %.2f +/- %.2f  (range %d-%d; %d with zero)",
  source, fragment, n_samples, mean_species, se_species, min_species, max_species, samples_with_zero)))

# ------------------------------------------------- 2. per station (R2.7) ----
per_station <- sp %>%
  group_by(source, station, fragment) %>%
  summarise(n_species = n_distinct(taxon), .groups = "drop") %>%
  bind_rows(sp %>% group_by(source, station) %>%
              summarise(n_species = n_distinct(taxon), .groups = "drop") %>% mutate(fragment = "either"))
write_csv(per_station, file.path(out_dir, "species_per_station.csv"))

# ------------------------------------- 3. 8.0 vs 9.5 kb overlap (R3.15) -----
overlap <- sp %>%
  distinct(source, taxon, fragment) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = fragment, values_from = present, values_fill = FALSE) %>%
  { if (!"8.0kb" %in% names(.)) mutate(., `8.0kb` = FALSE) else . } %>%
  { if (!"9.5kb" %in% names(.)) mutate(., `9.5kb` = FALSE) else . } %>%
  mutate(found_in = case_when(`8.0kb` & `9.5kb` ~ "both", `8.0kb` ~ "8.0kb only", TRUE ~ "9.5kb only"))
write_csv(overlap, file.path(out_dir, "fragment_overlap.csv"))
say("Species by fragment (8.0 vs 9.5 kb):")
ov <- overlap %>% count(source, found_in)
for (i in seq_len(nrow(ov))) say(sprintf("  %-8s %-11s %d", ov$source[i], ov$found_in[i], ov$n[i]))

# ------------------------------------------------- 4. PCR dropout (R3.7, R3.10)
pcr_long <- pcr %>%
  pivot_longer(c(`band_9.5kb`, `band_8.0kb`), names_to = "fragment", values_to = "band") %>%
  mutate(fragment = sub("band_", "", fragment), amplified = band == "y")
dropout <- pcr_long %>%
  group_by(fragment) %>%
  summarise(samples = n(), dropout = 1 - mean(amplified), amplified = sum(amplified), .groups = "drop")
by_station_pcr <- pcr_long %>%
  group_by(station, fragment) %>%
  summarise(reps_amplified = sum(amplified), .groups = "drop") %>%
  pivot_wider(names_from = fragment, values_from = reps_amplified, names_prefix = "reps_amplified_")
write_csv(dropout, file.path(out_dir, "pcr_dropout.csv"))
write_csv(by_station_pcr, file.path(out_dir, "pcr_success_by_station.csv"))
say("PCR dropout across 72 field samples (single PCR per sample, gel):")
for (i in seq_len(nrow(dropout))) with(dropout[i, ], say(sprintf(
  "  %s: %d of %d amplified (dropout %.0f%%)", fragment, amplified, samples, 100 * dropout)))
any_station <- pcr_long %>% group_by(station) %>% summarise(any = any(amplified))
say(sprintf("  Stations with >=1 band in any replicate or fragment: %d of %d",
            sum(any_station$any), nrow(any_station)))

# Gel band vs sequencing result for the sequenced field samples.
gel_vs_seq <- per_sample %>%
  filter(source == "Field", fragment %in% FRAGS) %>%
  left_join(pcr_long %>% select(station, replicate, fragment, amplified),
            by = c("station", "replicate", "fragment")) %>%
  mutate(species_detected = n_species > 0) %>%
  count(fragment, gel_band = amplified, species_detected)
write_csv(gel_vs_seq, file.path(out_dir, "gel_band_vs_sequencing.csv"))

# ------------------------------------------------------ 5. controls (R4.2) --
controls <- meta %>% filter(source == "NC") %>%
  left_join(summ %>% group_by(barcode) %>% summarise(reads_after_q = sum(reads), .groups = "drop"),
            by = "barcode") %>%
  left_join(counts %>% group_by(barcode) %>%
              summarise(taxa_any = n_distinct(taxon), assigned_reads = sum(reads), .groups = "drop"),
            by = "barcode") %>%
  mutate(across(c(reads_after_q, taxa_any, assigned_reads), ~ replace_na(.x, 0)))
write_csv(controls, file.path(out_dir, "controls.csv"))
say("Controls: ", paste(sprintf("%s %d reads, %d taxa", controls$sample_id,
                                 controls$reads_after_q, controls$taxa_any), collapse = "; "))

# ------------------------------------------ 6. genus-level reads (Sebastes) --
genus_level <- counts %>% filter(kept, rank == "genus") %>% inner_join(meta, by = "barcode") %>%
  group_by(source, sample_id, fragment, taxon) %>% summarise(reads = sum(reads), .groups = "drop")
write_csv(genus_level, file.path(out_dir, "genus_level_reads.csv"))

# ------------------------------------------- 7. aquarium vs census (R4.27) --
aq <- sp %>% filter(source == "Aquarium")
aq_det <- aq %>%
  group_by(taxon) %>%
  summarise(replicates_detected = n_distinct(barcode),
            reads_8.0kb = sum(reads[fragment == "8.0kb"]),
            reads_9.5kb = sum(reads[fragment == "9.5kb"]), .groups = "drop") %>%
  rename(species = taxon)
aq_table <- tank %>% full_join(aq_det, by = "species") %>%
  mutate(across(c(replicates_detected, reads_8.0kb, reads_9.5kb), ~ replace_na(.x, 0)),
         status = case_when(is.na(tank_census) ~ "not on Aquarium list",
                            tank_census == 0 ~ "on list, census 0",
                            TRUE ~ "in tank")) %>%
  arrange(status, desc(tank_census), desc(reads_8.0kb + reads_9.5kb))
write_csv(aq_table, file.path(out_dir, "aquarium_detection_vs_census.csv"))
in_tank <- aq_table %>% filter(status == "in tank")
say(sprintf("Aquarium: %d of %d tank species detected", sum(in_tank$replicates_detected > 0), nrow(in_tank)))
for (lo in list(c(4, Inf), c(1, 3))) {
  s <- in_tank %>% filter(tank_census >= lo[1], tank_census <= lo[2])
  say(sprintf("  census %s: %d of %d detected", ifelse(is.infinite(lo[2]), ">=4", "1-3"),
              sum(s$replicates_detected > 0), nrow(s)))
}
say("  Off-list species: ", paste(aq_table$species[aq_table$status == "not on Aquarium list"], collapse = ", "))

writeLines(log_lines, file.path(out_dir, "stats_summary.txt"))

# ================================================================ FIGURES ====
theme_ms <- theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = "italic"), legend.position = "right")

sample_order <- meta %>% filter(source %in% c("Field", "Aquarium")) %>%
  mutate(o = match(barcode, meta$barcode)) %>% arrange(o) %>% pull(sample_id)
sp_fig <- sp %>% filter(source %in% c("Field", "Aquarium"))   # controls never plotted
species_order <- sp_fig %>% group_by(taxon) %>% summarise(r = sum(reads)) %>% arrange(r) %>% pull(taxon)

# Fig. 6 (field) and Fig. 7 (aquarium): species x sample, fill = log10 reads, one panel per
# fragment. Each figure lists only the species detected in that sample set, ordered by total
# reads. Every sample gets a column, including those with no species detected (shown blank).
# The two figures share one colour scale so they can be compared.
reads_max <- max(sp_fig$reads)
heatmap_for <- function(src) {
  smp <- samples %>% filter(source == src) %>% select(sample_id)
  sp_src <- sp_fig %>% filter(source == src)
  sp_order <- sp_src %>% group_by(taxon) %>% summarise(r = sum(reads)) %>% arrange(r) %>% pull(taxon)
  heat <- smp %>% crossing(fragment = FRAGS, taxon = sp_order) %>%
    left_join(sp_src %>% group_by(sample_id, fragment, taxon) %>%
                summarise(reads = sum(reads), .groups = "drop"),
              by = c("sample_id", "fragment", "taxon")) %>%
    mutate(sample_id = factor(sample_id, levels = intersect(sample_order, smp$sample_id)),
           taxon = factor(taxon, levels = sp_order),
           fragment = recode(fragment, "8.0kb" = "8.0 kb fragment", "9.5kb" = "9.5 kb fragment"))
  p <- ggplot(heat, aes(sample_id, taxon, fill = reads)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    facet_wrap(~ fragment, nrow = 1) +
    scale_fill_gradient(name = "Reads\n(log10 scale)", trans = "log10",
                        low = "#9ecae1", high = "#08306b", na.value = "grey95",
                        limits = c(1, reads_max),
                        breaks = c(1, 10, 100, 1000, 10000),
                        labels = c("1", "10", "100", "1,000", "10,000")) +
    labs(x = "Sample", y = NULL) + theme_ms +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          strip.text = element_text(size = 8))
  list(plot = p, n_sp = length(sp_order))
}
for (f in list(list(src = "Field", name = "fig6_heatmap_field", w = 180),
               list(src = "Aquarium", name = "fig7_heatmap_aquarium", w = 140))) {
  h <- heatmap_for(f$src)
  ht <- 35 + 3.6 * h$n_sp
  ggsave(file.path(fig_dir, paste0(f$name, ".pdf")), h$plot, width = f$w, height = ht, units = "mm")
  ggsave(file.path(fig_dir, paste0(f$name, ".png")), h$plot, width = f$w, height = ht, units = "mm", dpi = 300)
}
# The occupancy figure (former Fig. 7) was dropped in favour of the two heatmaps.

message("\nWrote tables to ", out_dir, " and figures to ", fig_dir)
message("R ", R.version$major, ".", R.version$minor, "; ",
        paste(sapply(c("dplyr", "tidyr", "readr", "stringr", "ggplot2"),
                     function(p) paste(p, as.character(packageVersion(p)))), collapse = "; "))
