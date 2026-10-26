# Supplementary figures: replicate occupancy vs. mean reads per positive replicate
#
#   Fig A  LR-PCR, field and aquarium: replicates with detection (of 3) vs.
#          log10 mean reads per positive replicate, with a binomial (3-trial) fit.
#   Fig B  Aquarium only: MiFish vs LR-PCR, same axes and model, fitted per assay,
#          plus a test of whether the occupancy-reads relationship differs by assay.
#
# Inputs: long-format tables, one row per sample x (fragment) x species with reads > 0
#   LR-PCR: results/assignment/species_by_sample.tsv   (columns: sample, fragment, species, reads)
#   MiFish: MiFish species_by_sample.tsv for the aquarium (columns: sample, species, reads)
# Edit CONFIG if your column names or paths differ.
#
# Usage: Rscript figS_occurrence_vs_reads.R [lrpcr.tsv] [mifish.tsv] [output_prefix]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---- CONFIG -----------------------------------------------------------------
args       <- commandArgs(trailingOnly = TRUE)
lr_file    <- if (length(args) >= 1) args[1] else "results/assignment/species_by_sample.tsv"
mf_file    <- if (length(args) >= 2) args[2] else "results_mifish/assignment/species_by_sample.tsv"
out_pref   <- if (length(args) >= 3) args[3] else "FigS_occupancy"

col_sample   <- "sample"
col_fragment <- "fragment"   # LR-PCR only
col_species  <- "species"
col_reads    <- "reads"

drop_samples <- c("EB-4/29", "NTC-8kb", "NTC-9.5kb")  # negative controls
meta_file    <- "reference/sample_metadata.csv"         # only used if samples are barcode names
# Fallback if meta_file isn't found: LR-PCR barcode -> sample ID (from sample_metadata.csv)
barcode_key <- c(
  barcode01 = "P28-A", barcode02 = "P28-B", barcode03 = "P28-C",
  barcode04 = "P1-A",  barcode05 = "P1-B",  barcode06 = "P1-C",
  barcode07 = "P3-A",  barcode08 = "P3-B",  barcode09 = "P3-C",
  barcode10 = "P4-A",  barcode11 = "P4-B",  barcode12 = "P4-C",
  barcode13 = "P7-A",  barcode14 = "P7-B",  barcode15 = "P7-C",
  barcode16 = "P38-A", barcode17 = "P38-B", barcode18 = "P38-C",
  barcode19 = "EB-4/29",
  barcode20 = "AQ202408-A", barcode21 = "AQ202408-B", barcode22 = "AQ202408-C",
  barcode23 = "NTC-8kb", barcode24 = "NTC-9.5kb")
aq_pattern   <- "^AQ"         # LR-PCR samples from the aquarium

# MiFish samples to use (NULL = all samples in the MiFish file).
# Set to the aquarium replicate IDs if the file also holds other samples
# (either the PCR-replicate names or the collapsed biological-sample names work).
mf_samples   <- NULL

# MiFish PCR replicates: each biological replicate (18_mifish_bt1, 19_mifish_bt2,
# 20_mifish_bt3) was sequenced in triplicate (_a_all, _b_all, _c_all). The suffix
# matched by this pattern is stripped, so the PCR replicates collapse into one
# biological sample. A biological replicate is positive if any of its PCR
# replicates detected the species. (The trailing "_all" is optional, so names
# ending in just _a/_b/_c also work.)
mf_pcr_suffix   <- "_[abc](_all)?$"
# How to combine reads across PCR replicates: "sum" (pooled reads) or "mean"
# (closer to the single-PCR depth of each LR-PCR replicate)
mf_pcr_collapse <- "sum"

# In Fig B, LR-PCR is scored per species: a replicate is positive if either
# fragment detected the species, and its reads are the sum over both fragments.
# Set FALSE to show the 8.0 and 9.5 kb fragments as separate assays instead.
lr_combine_fragments <- TRUE

# Optional: text file with one tank-census species name per line. If present,
# Fig B uses filled points for tank species and open points for off-list species.
tank_list_file <- "tank_species.txt"

frag_cols  <- c("8.0 kb" = "#1f3b73", "9.5 kb" = "#8cc084")
assay_cols <- c("MiFish" = "#c2185b", "LR-PCR" = "#1f3b73",
                "LR-PCR 8.0 kb" = "#1f3b73", "LR-PCR 9.5 kb" = "#8cc084")
# -----------------------------------------------------------------------------

# Column lookup: the CONFIG name first, then common alternatives (case-insensitive)
alt_names <- list(
  sample   = c("sample", "sample_id", "sampleid", "SampleID", "barcode", "library"),
  species  = c("species", "taxon", "taxonomy", "lca", "assignment", "sci_name", "scientific_name"),
  reads    = c("reads", "n_reads", "read_count", "count", "nreads", "n"),
  fragment = c("fragment", "amplicon", "frag", "target")
)
find_col <- function(d, preferred, role, f) {
  cands <- c(preferred, alt_names[[role]])
  hit <- names(d)[match(tolower(cands), tolower(names(d)))]
  hit <- hit[!is.na(hit)][1]
  if (is.na(hit)) stop("Could not find the '", role, "' column in ", f,
                       ".\n  Columns present: ", paste(names(d), collapse = ", "),
                       "\n  Set col_", role, " in the CONFIG block to the right name.", call. = FALSE)
  hit
}

read_long <- function(f, has_fragment, map_barcodes = FALSE) {
  d <- read.delim(f, check.names = FALSE)
  cp <- find_col(d, col_species, "species", f)

  # Wide format (one column of read counts per sample, e.g. barcode01, barcode02 ...):
  # if no sample column exists, every numeric column other than species/fragment
  # is treated as a sample and the table is reshaped to long format.
  has_sample_col <- any(tolower(c(col_sample, alt_names$sample)) %in% tolower(names(d)))
  if (!has_sample_col) {
    id_cols <- c(cp, if (has_fragment) find_col(d, col_fragment, "fragment", f))
    smp_cols <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], id_cols)
    if (length(smp_cols) == 0)
      stop("No sample column and no numeric per-sample columns in ", f,
           ".\n  Columns present: ", paste(names(d), collapse = ", "), call. = FALSE)
    message(basename(f), ": wide format, reshaping ", length(smp_cols), " sample columns (",
            paste(head(smp_cols, 3), collapse = ", "), if (length(smp_cols) > 3) ", ..." else "", ")")
    d <- pivot_longer(d, all_of(smp_cols), names_to = ".sample", values_to = ".reads")
    cs <- ".sample"; cr <- ".reads"
  } else {
    cs <- find_col(d, col_sample, "sample", f)
    cr <- find_col(d, col_reads,  "reads",  f)
  }
  message(basename(f), ": using columns sample='", cs, "', species='", cp, "', reads='", cr, "'",
          if (has_fragment) paste0(", fragment='", find_col(d, col_fragment, "fragment", f), "'") else "")
  out <- data.frame(sample  = as.character(d[[cs]]),
                    species = as.character(d[[cp]]),
                    reads   = as.numeric(d[[cr]]))
  if (has_fragment) {
    fr <- as.character(d[[find_col(d, col_fragment, "fragment", f)]])
    # normalise "8.0kb", "8kb", "8.0 kb", "8.0" -> "8.0 kb"; same for 9.5
    fr <- ifelse(grepl("^\\s*8(\\.0)?\\s*(kb)?\\s*$", fr, ignore.case = TRUE), "8.0 kb",
          ifelse(grepl("^\\s*9\\.5\\s*(kb)?\\s*$",     fr, ignore.case = TRUE), "9.5 kb", fr))
    unknown <- setdiff(unique(fr), names(frag_cols))
    if (length(unknown)) warning("Unrecognised fragment labels: ", paste(unknown, collapse = ", "))
    out$fragment <- fr
  } else {
    out$fragment <- NA_character_
  }
  # If samples are barcode names (barcode01...), translate to sample IDs (P28-A...)
  # (LR-PCR run only; MiFish barcode names are a different run and are left as-is)
  if (map_barcodes && all(grepl("^barcode\\d+$", out$sample, ignore.case = TRUE))) {
    if (file.exists(meta_file)) {
      meta <- read.csv(meta_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
      key  <- setNames(trimws(meta$SampleID), tolower(trimws(meta$demulti)))
      src  <- meta_file
    } else {
      key <- barcode_key; src <- "built-in table (from sample_metadata.csv)"
    }
    unmapped <- setdiff(unique(tolower(out$sample)), names(key))
    if (length(unmapped)) stop("No sample ID for: ", paste(unmapped, collapse = ", "), call. = FALSE)
    out$sample <- unname(key[tolower(out$sample)])
    message("  translated barcode names to sample IDs using ", src)
  }
  out <- filter(out, !sample %in% drop_samples, !is.na(sample), !is.na(reads), reads > 0)
  if (nrow(out) == 0) stop("No rows with reads > 0 left in ", f, call. = FALSE)
  out
}

# Occupancy per group: n_occ = positive replicates, n_trials = replicates sampled,
# mean_reads = mean reads over the positive replicates only.
occupancy <- function(d, n_trials_tbl, ...) {
  d |>
    group_by(..., station, sample, species) |>
    summarise(reads = sum(reads), .groups = "drop") |>          # sums fragments if not grouped by them
    group_by(..., station, species) |>
    summarise(n_occ = n_distinct(sample), mean_reads = mean(reads), .groups = "drop") |>
    left_join(n_trials_tbl, by = "station") |>
    mutate(log_reads = log10(mean_reads))
}

# Binomial GLM with n trials; returns model, LRT vs null, and fitted curve + 95% CI
fit_binom <- function(d) {
  m  <- glm(cbind(n_occ, n_trials - n_occ) ~ log_reads, family = binomial, data = d)
  m0 <- update(m, . ~ 1)
  grid <- data.frame(log_reads = seq(min(d$log_reads), max(d$log_reads), length.out = 200))
  pr   <- predict(m, newdata = grid, type = "link", se.fit = TRUE)
  n    <- unique(d$n_trials)
  if (length(n) != 1) stop("Groups with different replicate numbers must be fit separately.")
  grid$fit <- n * plogis(pr$fit)
  grid$lo  <- n * plogis(pr$fit - 1.96 * pr$se.fit)
  grid$hi  <- n * plogis(pr$fit + 1.96 * pr$se.fit)
  grid$n_trials <- n
  list(model = m, lrt = anova(m0, m, test = "Chisq"), curve = grid)
}

report <- function(label, f, d) {
  cat("====", label, "====\n")
  cat("n =", nrow(d), "species x station (x fragment) combinations; replicates per station =",
      paste(unique(d$n_trials), collapse = ","), "\n")
  print(summary(f$model)$coefficients)
  cat("Residual deviance:", round(deviance(f$model), 2), "on", df.residual(f$model), "df\n")
  cat("Likelihood-ratio test vs intercept-only:\n"); print(f$lrt); cat("\n")
}

reads_axis <- scale_x_continuous(labels = function(x) scales::label_comma()(10^x))
base_theme <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"))

# ===== Fig A: LR-PCR, field and aquarium ======================================
lr <- read_long(lr_file, has_fragment = TRUE, map_barcodes = TRUE) |>
  mutate(station = sub("-[A-Z]$", "", sample),
         source  = ifelse(grepl(aq_pattern, station), "Seattle Aquarium", "Puget Sound (field)"))

# Three biological replicates per station, including replicates with no reads at all
lr_trials <- data.frame(station = unique(lr$station), n_trials = 3)

occA <- occupancy(lr, lr_trials, source, fragment)
fitsA <- split(occA, occA$source) |> lapply(fit_binom)
curvesA <- bind_rows(lapply(names(fitsA), \(s) mutate(fitsA[[s]]$curve, source = s)))

figA <- ggplot(occA, aes(log_reads, n_occ)) +
  geom_ribbon(data = curvesA, aes(y = fit, ymin = lo, ymax = hi), fill = "grey80", alpha = 0.6) +
  geom_line(data = curvesA, aes(y = fit), linewidth = 0.9, colour = "grey20") +
  geom_point(aes(colour = fragment), alpha = 0.75, size = 2,
             position = position_jitter(width = 0, height = 0.08, seed = 1)) +
  facet_wrap(~ source, scales = "free_x") +
  scale_y_continuous(breaks = 0:3, limits = c(0, 3.2)) +
  reads_axis +
  scale_colour_manual(values = frag_cols, name = "LR-PCR fragment") +
  labs(x = "Mean reads per positive replicate (log10 scale)",
       y = "Biological replicates with detection (of 3)") +
  base_theme

ggsave(paste0(out_pref, "_A_LRPCR.pdf"), figA, width = 7.5, height = 4)
ggsave(paste0(out_pref, "_A_LRPCR.png"), figA, width = 7.5, height = 4, dpi = 300)

# ===== Fig B: aquarium, MiFish vs LR-PCR ======================================
mf <- read_long(mf_file, has_fragment = FALSE)
mf <- mutate(mf, bio_sample = sub(mf_pcr_suffix, "", sample))
if (!is.null(mf_samples)) mf <- filter(mf, sample %in% mf_samples | bio_sample %in% mf_samples)

# Collapse PCR replicates within each biological sample
pcr_per_bio <- mf |> distinct(bio_sample, sample) |> count(bio_sample, name = "n_pcr")
if (all(pcr_per_bio$n_pcr == 1))
  warning("No MiFish PCR replicates were collapsed: check mf_pcr_suffix against the sample names (",
          paste(head(unique(mf$sample), 3), collapse = ", "), ", ...)", call. = FALSE)
message("MiFish PCR replicates collapsed (", mf_pcr_collapse, " of reads): ",
        paste0(pcr_per_bio$bio_sample, " [", pcr_per_bio$n_pcr, " PCR]", collapse = ", "))
n_pcr_total <- max(pcr_per_bio$n_pcr)   # replicates with zero reads for a species count as 0 in the mean
mf <- mf |>
  group_by(bio_sample, species) |>
  summarise(reads = if (mf_pcr_collapse == "mean") sum(reads) / n_pcr_total else sum(reads),
            .groups = "drop") |>
  transmute(sample = bio_sample, species, reads, fragment = NA_character_,
            station = "AQ", assay = "MiFish")

lr_aq <- lr |> filter(source == "Seattle Aquarium") |> mutate(station = "AQ")
lr_aq <- if (lr_combine_fragments) {
  mutate(lr_aq, assay = "LR-PCR")
} else {
  mutate(lr_aq, assay = paste("LR-PCR", fragment))
}

n_mf <- n_distinct(mf$sample)
n_lr <- 3
message("MiFish replicates found: ", n_mf, " (", paste(sort(unique(mf$sample)), collapse = ", "), ")")

occB <- bind_rows(
  occupancy(select(mf,    -fragment), data.frame(station = "AQ", n_trials = n_mf), assay),
  occupancy(select(lr_aq, -fragment), data.frame(station = "AQ", n_trials = n_lr), assay)
) |>
  mutate(assay = factor(assay, levels = intersect(names(assay_cols), unique(assay))),
         prop  = n_occ / n_trials)

if (file.exists(tank_list_file)) {
  tank <- trimws(readLines(tank_list_file)); tank <- tank[tank != ""]
  occB$tank <- ifelse(occB$species %in% tank, "Tank census", "Not on census")
} else {
  occB$tank <- "All species"
}

fitsB   <- split(occB, occB$assay, drop = TRUE) |> lapply(fit_binom)
curvesB <- bind_rows(lapply(names(fitsB), \(a) mutate(fitsB[[a]]$curve, assay = a))) |>
  mutate(assay = factor(assay, levels = levels(occB$assay)))

# Does the occupancy-reads relationship differ between assays? (proportion scale,
# so valid even if MiFish has a different number of replicates)
mB_add <- glm(cbind(n_occ, n_trials - n_occ) ~ log_reads + assay, family = binomial, data = occB)
mB_int <- glm(cbind(n_occ, n_trials - n_occ) ~ log_reads * assay, family = binomial, data = occB)
mB_one <- glm(cbind(n_occ, n_trials - n_occ) ~ log_reads,         family = binomial, data = occB)

# y-axis: replicate count if both assays have the same number of replicates, else proportion
same_n <- length(unique(occB$n_trials)) == 1
if (same_n) {
  occB$y <- occB$n_occ
  curvesB <- mutate(curvesB, y = fit, ylo = lo, yhi = hi)
  yscale <- scale_y_continuous(breaks = 0:n_lr, limits = c(0, n_lr + 0.2))
  ylab   <- paste0("Biological replicates with detection (of ", n_lr, ")")
  jit    <- 0.08
} else {
  occB$y <- occB$prop
  curvesB <- mutate(curvesB, y = fit / n_trials, ylo = lo / n_trials, yhi = hi / n_trials)
  yscale <- scale_y_continuous(limits = c(0, 1.07), labels = scales::label_percent())
  ylab   <- "Replicates with detection"
  jit    <- 0.025
}

figB <- ggplot(occB, aes(log_reads, y, colour = assay)) +
  geom_ribbon(data = curvesB, aes(y = y, ymin = ylo, ymax = yhi, fill = assay),
              colour = NA, alpha = 0.15) +
  geom_line(data = curvesB, aes(y = y), linewidth = 0.9) +
  geom_point(aes(shape = tank), alpha = 0.8, size = 2.2, stroke = 0.8,
             position = position_jitter(width = 0, height = jit, seed = 1)) +
  yscale + reads_axis +
  scale_colour_manual(values = assay_cols, name = NULL) +
  scale_fill_manual(values = assay_cols, guide = "none") +
  scale_shape_manual(values = c("Tank census" = 16, "Not on census" = 1, "All species" = 16),
                     name = NULL, guide = if (file.exists(tank_list_file)) "legend" else "none") +
  labs(x = "Mean reads per positive replicate (log10 scale)", y = ylab) +
  base_theme

ggsave(paste0(out_pref, "_B_MiFish_vs_LRPCR.pdf"), figB, width = 5.5, height = 4.5)
ggsave(paste0(out_pref, "_B_MiFish_vs_LRPCR.png"), figB, width = 5.5, height = 4.5, dpi = 300)

# ===== Model summaries ========================================================
sink(paste0(out_pref, "_model_summary.txt"))
cat("Binomial GLMs: cbind(replicates detected, replicates not detected) ~ log10(mean reads per positive replicate)\n\n")
cat("######## Fig A: LR-PCR ########\n\n")
for (s in names(fitsA)) report(s, fitsA[[s]], filter(occA, source == s))
cat("######## Fig B: aquarium, MiFish vs LR-PCR ########\n\n")
for (a in names(fitsB)) report(a, fitsB[[a]], filter(occB, assay == a))
cat("Assay effect on intercept (log_reads vs log_reads + assay):\n")
print(anova(mB_one, mB_add, test = "Chisq"))
cat("\nAssay effect on slope (log_reads + assay vs log_reads * assay):\n")
print(anova(mB_add, mB_int, test = "Chisq"))
cat("\nInteraction model coefficients:\n")
print(summary(mB_int)$coefficients)
sink()

write.table(occA, paste0(out_pref, "_A_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(occB, paste0(out_pref, "_B_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

message("Wrote ", out_pref, "_A_LRPCR.pdf/.png, ", out_pref, "_B_MiFish_vs_LRPCR.pdf/.png, ",
        out_pref, "_model_summary.txt and the plotted data tables")
