# =============================================================================
# figureS2_gse_band_by_year.R
# Figure S2: values of the GSE band A07 at the seagrass presence samples
# (PA = 1) in each survey year, with the points coloured by region. A07 is
# shown as an example band.
#
# Inputs  (result/)  dataPA_<year>.csv for 2018-2023 (R/workflow/01_prepare_point_data.R)
# Outputs (result/paper_tables/)  FigureS2_gse_A07_by_year.png / .tiff / .pdf / .svg
#
# The page size (10.42 x 5.21 inches) and the font size (11 pt) are those of the
# figure in the Supplementary Material, which was saved at 1000 x 500 pixels
# (96 dpi). The points are jittered with a fixed seed, so their positions are
# the same in every run. For a higher resolution, set fig_dpi before running, e.g.
#   fig_dpi <- 600; source(here::here("R", "figures_tables", "figureS2_gse_band_by_year.R"))
# =============================================================================
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)

in_dir  <- if (exists("in_dir"))  in_dir  else here::here("result")
out_dir <- if (exists("out_dir")) out_dir else here::here("result", "paper_tables")
fig_dpi <- if (exists("fig_dpi")) fig_dpi else 300
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

years <- 2018:2023
band  <- "GSE_A07"

# -----------------------------------------------------------------------------
# Presence samples of each survey year
# -----------------------------------------------------------------------------
pa1 <- map_dfr(years, function(yr) {
  path <- file.path(in_dir, paste0("dataPA_", yr, ".csv"))
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE) %>%
    filter(PA == 1, !is.na(.data[[band]]), !is.na(loc)) %>%
    transmute(year = factor(yr, levels = years), loc, value = .data[[band]])
})

cat(sprintf("Presence samples: %d\n", nrow(pa1)))
print(table(pa1$year, pa1$loc))

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
set.seed(42)
p <- ggplot(pa1, aes(x = year, y = value, colour = loc)) +
  geom_boxplot(fill = "#66CCFF", outlier.shape = NA, alpha = 0.6, colour = "grey30") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 1.5) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Year", y = "Pixel Value", colour = "Location") +
  theme_minimal(base_size = 11, base_family = "Arial")

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
w <- 1000 / 96; h <- 500 / 96
stem <- file.path(out_dir, "FigureS2_gse_A07_by_year")
ggsave(paste0(stem, ".png"),  p, width = w, height = h, dpi = fig_dpi, bg = "white",
       device = ragg::agg_png)
ggsave(paste0(stem, ".tiff"), p, width = w, height = h, dpi = fig_dpi, bg = "white",
       device = ragg::agg_tiff, compression = "lzw")
ggsave(paste0(stem, ".pdf"),  p, width = w, height = h, bg = "white", device = cairo_pdf)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(stem, ".svg"), p, width = w, height = h, bg = "white")
}
cat("Figure S2 saved in", out_dir, "\n")
