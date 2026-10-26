#!/usr/bin/env Rscript
# =============================================================================
# Fig. 1 - map of the WOAC stations sampled for eDNA, with a West Coast inset.
#   Main map: central and southern Puget Sound. Each station is shaded with
#   topography and bathymetry, and its symbol shows whether its samples were
#   sequenced on ONT or only tested by PCR.
#   Inset: the North American West Coast, with a box around the main map.
#   Both panels use the same depth colour scale.
#
# Inputs (lists/):
#   lists/pcr_success_tableS1.csv   the 72 field samples (24 stations) tested by
#                                   LR-PCR (Table S1); these are the stations mapped
#   lists/sample_metadata.csv       the sequenced samples; field stations here are
#                                   "ONT sequencing", the rest "PCR only"
#   lists/WOAC_cruise_plan.csv      station coordinates (columns StationID, Lat, Lon).
#                                   Copy it from data/WOAC_cruise_plan.csv in the
#                                   manuscript folder. Sequenced stations fall back to
#                                   the lat/lon in sample_metadata.csv.
#
# Bathymetry: NOAA ETOPO, downloaded by marmap::getNOAA.bathy. The first run needs
#   an internet connection; the grids are saved in data/bathy/ and reused after that.
#
# Outputs:
#   results/figures/fig1_sampling_map.pdf / .png / .tiff   (6 x 6 in; TIFF is LZW, 300 dpi)
#
# Reviewer comments addressed: R4.30 (legend: ONT sequencing vs PCR only).
#
# Usage (from the lrpcr_blast folder):  Rscript 12_fig1_sampling_map.R
#   or give the folder:                 Rscript 12_fig1_sampling_map.R /path/to/lrpcr_blast
# Packages: dplyr, readr, stringr, tibble, ggplot2 (tidyverse), sf, marmap, ggrepel
#   (install.packages(c("sf", "marmap", "ggrepel")) if needed).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tibble); library(ggplot2)
  library(sf); library(marmap); library(ggrepel)
})
set.seed(110011)                      # fixes the ggrepel label positions

args <- commandArgs(trailingOnly = TRUE)
proj <- if (length(args) >= 1) args[1] else "."
list_dir  <- file.path(proj, "lists")
fig_dir   <- file.path(proj, "results", "figures")
bathy_dir <- file.path(proj, "data", "bathy")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(bathy_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Map extents and settings ------------------------------------------------
MAIN_X <- c(-123.7, -122.1); MAIN_Y <- c(47.1, 48.4)   # main map
WC_X   <- c(-129, -122);     WC_Y   <- c(39, 53)       # inset
RES_MAIN <- 0.1    # bathymetry grid resolution, arc-minutes (main map)
RES_WC   <- 3      # (inset)
BUF      <- 0.3    # degrees of extra bathymetry around the main map, so contours reach the edge
DEPTH_COLS <- c(low = "#0D0D47", mid = "gainsboro", high = "palegreen4")
SEQUENCED  <- "ONT sequencing"; PCR_ONLY <- "PCR only"
SEATTLE <- tibble(landmark = "Seattle, WA", latitude = 47.609153, longitude = -122.341500)

# ---- Stations ----------------------------------------------------------------
need <- function(f) { if (!file.exists(f)) stop("Not found: ", f, call. = FALSE); f }

pcr <- read_csv(need(file.path(list_dir, "pcr_success_tableS1.csv")), show_col_types = FALSE)
meta <- read_csv(need(file.path(list_dir, "sample_metadata.csv")), show_col_types = FALSE,
                 locale = locale(encoding = "UTF-8")) %>%
  rename_with(~ str_remove(.x, "^﻿"))
seq_st <- meta %>% filter(Source == "Field") %>%
  transmute(StationID = paste0("P", Station_num), Lat = lat, Lon = lon) %>% distinct()

cruise_file <- file.path(list_dir, "WOAC_cruise_plan.csv")
cruise <- if (file.exists(cruise_file)) {
  read_csv(cruise_file, show_col_types = FALSE) %>%
    transmute(StationID = str_trim(StationID), Lat = as.numeric(Lat), Lon = as.numeric(Lon))
} else {
  message("Note: ", cruise_file, " not found; only the sequenced stations have coordinates.")
  tibble(StationID = character(), Lat = numeric(), Lon = numeric())
}

stations <- tibble(StationID = unique(str_trim(pcr$station))) %>%
  left_join(cruise, by = "StationID") %>%
  left_join(seq_st, by = "StationID", suffix = c("", ".seq")) %>%
  mutate(Lat = coalesce(Lat, Lat.seq), Lon = coalesce(Lon, Lon.seq),
         Test = factor(if_else(StationID %in% seq_st$StationID, SEQUENCED, PCR_ONLY),
                       levels = c(SEQUENCED, PCR_ONLY))) %>%
  select(StationID, Lat, Lon, Test)

missing <- stations %>% filter(is.na(Lat) | is.na(Lon))
if (nrow(missing) > 0)
  warning("No coordinates for ", nrow(missing), " station(s), left off the map: ",
          paste(missing$StationID, collapse = ", "), call. = FALSE)
stations <- stations %>% filter(!is.na(Lat), !is.na(Lon))
message(nrow(stations), " stations mapped: ", sum(stations$Test == SEQUENCED), " sequenced, ",
        sum(stations$Test == PCR_ONLY), " PCR only.")

# ---- Bathymetry --------------------------------------------------------------
# Land cells are set to +|deepest point in the inset|, so land is one green in both
# panels, and both fill scales share the limits [deepest, -deepest], centred on 0.
bathy <- function(x, y, res) {
  b <- getNOAA.bathy(x[1], x[2], y[1], y[2], resolution = res, antimeridian = FALSE,
                     keep = TRUE, path = bathy_dir)
  as.xyz(b) %>% setNames(c("lon", "lat", "z"))
}
wc   <- bathy(WC_X, WC_Y, RES_WC)
main <- bathy(MAIN_X + c(-BUF, BUF), MAIN_Y + c(-BUF, BUF), RES_MAIN)
deepest <- min(wc$z, main$z)
wc$z[wc$z > -1]     <- -deepest
main$z[main$z > -1] <- -deepest
fill_scale <- scale_fill_gradient2(low = DEPTH_COLS[["low"]], mid = DEPTH_COLS[["mid"]],
                                   high = DEPTH_COLS[["high"]], midpoint = 0,
                                   limits = c(deepest, -deepest), guide = "none")

map_theme <- theme_classic() +
  theme(axis.line = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank())

base_map <- function(d, contour_step) {
  ggplot() +
    geom_tile(data = d, aes(x = lon, y = lat, fill = z)) +
    fill_scale +
    geom_contour(data = d, aes(x = lon, y = lat, z = z), breaks = 0,
                 colour = "grey5", linewidth = 0.3) +
    geom_contour(data = d, aes(x = lon, y = lat, z = z), binwidth = contour_step,
                 colour = "grey65", linewidth = 0.1)
}

# ---- Inset: West Coast -------------------------------------------------------
p_wc <- base_map(wc, 1000) +
  annotate("rect", xmin = MAIN_X[1], xmax = MAIN_X[2], ymin = MAIN_Y[1], ymax = MAIN_Y[2],
           alpha = 0.1, colour = "black", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(-129, -122, by = 3)) +
  scale_y_continuous(breaks = seq(40, 52, by = 4)) +
  coord_sf(xlim = WC_X, ylim = WC_Y, crs = 4326) +
  labs(x = NULL, y = NULL) +
  map_theme +
  theme(panel.background = element_rect(fill = "grey90"),
        plot.background  = element_rect(fill = "grey90"))

# ---- Main map: Puget Sound stations ------------------------------------------
p_main <- base_map(main, 100) +
  geom_point(data = stations, aes(x = Lon, y = Lat, shape = Test), size = 2) +
  geom_label_repel(data = stations, aes(x = Lon, y = Lat, label = StationID),
                   point.padding = 0.25, alpha = 0.8, min.segment.length = 0, size = 3,
                   label.padding = 0.1, box.padding = 0.1) +
  geom_point(data = SEATTLE, aes(x = longitude, y = latitude), shape = 1, size = 5) +
  geom_label_repel(data = SEATTLE, aes(x = longitude, y = latitude, label = landmark),
                   point.padding = 5, min.segment.length = 0, size = 3,
                   label.padding = 0.2, box.padding = 0.2, nudge_x = 1, nudge_y = 0) +
  scale_shape_manual(values = setNames(c(16, 17), c(SEQUENCED, PCR_ONLY)), drop = FALSE) +
  coord_sf(xlim = MAIN_X, ylim = MAIN_Y, crs = 4326) +
  labs(x = "Longitude", y = "Latitude", shape = "Method") +
  map_theme +
  # legend inside the map, top left ("inside" needs ggplot2 >= 3.5.0)
  (if (packageVersion("ggplot2") >= "3.5.0")
     theme(legend.position = "inside", legend.position.inside = c(0.17, 0.82))
   else theme(legend.position = c(0.17, 0.82))) +
  theme(legend.background = element_rect(fill = "grey99", linewidth = 0.2, linetype = "solid"),
        legend.text  = element_text(size = rel(0.7)),
        legend.title = element_text(size = rel(0.7)),
        legend.margin = margin(3, 3, 3, 3))

p <- p_main + annotation_custom(ggplotGrob(p_wc), xmin = -123.7, xmax = -123.25,
                                ymin = 47.15, ymax = 47.9)

# ---- Save --------------------------------------------------------------------
out <- file.path(fig_dir, "fig1_sampling_map")
ggsave(paste0(out, ".pdf"),  p, width = 6, height = 6, units = "in")
ggsave(paste0(out, ".png"),  p, width = 6, height = 6, units = "in", dpi = 300)
ggsave(paste0(out, ".tiff"), p, width = 6, height = 6, units = "in", dpi = 300, compression = "lzw")
message("Wrote ", out, ".(pdf|png|tiff)")
