# =============================================================================
# figure6_accuracy_map.R
# Figure 6: predicted AGC and its expected accuracy for one region and survey
# year (Rote 2021 in the paper).
#   (a) Expected accuracy class of the mapped seagrass extent, from the distance
#       to the nearest AGC sample of that year. Shown at 100 m, with the best
#       class of each 100 m cell, so that thin coastal meadows stay visible.
#   (b) Zoom on the sampled coast: expected accuracy class at 10 m.
#   (c) Zoom: AGC predicted by the proxy workflow (10 m), with the samples.
# The MEC in the legend is the MEC of the GSE-only model in this region
# (Table S8) for the design whose median distance falls in each class.
# The layout in the paper was arranged in QGIS from the same layers.
#
# Inputs (result/gee/, exports of gee/12_figure6_layers.js)
#   Fig6_<REGION><YEAR>_layers.tif, Fig6_<REGION><YEAR>_S2_RGB.tif,
#   Fig6_<REGION><YEAR>_samples.csv
# Outputs (result/paper_tables/)  Figure6_<REGION><YEAR>_map.png / .tiff / .pdf / .svg
# =============================================================================
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))   # degree sign, en dash
suppressMessages({
  library(terra)
  library(ggplot2)
  library(patchwork)
})

in_dir  <- if (exists("in_dir"))  in_dir  else here::here("result", "gee")
out_dir <- if (exists("out_dir")) out_dir else here::here("result", "paper_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- CONFIGURATION (only this block changes between regions) -----------------
cfg <- list(
  name   = "Rote2021", year = 2021, crs = "EPSG:32751", lat0 = -10.80,
  ov_ext = c(458000, 551000, 8780000, 8850000),        # all mapped seagrass of Rote
  zm_ext = c(477500, 486500, 8789500, 8799500),        # west coast, 1,457 of 1,461 samples
  ov_bar = 20, zm_bar = 2,                             # scale bars (km)
  ov_xb = seq(122.6, 123.4, 0.2), ov_yb = seq(-11.0, -10.4, 0.2),
  zm_xb = seq(122.78, 122.90, 0.02), zm_yb = seq(-10.94, -10.88, 0.02),
  # Rote, GSE-only (Table S8); median distances in Rote: pixel 6 m, 100 m blocks
  # 41 m, 250 m blocks 89 m, 500 m blocks 153 m, 1000 m blocks 452 m
  cls_labs = c("1" = "≤ 10 m (MEC 0.43)",
               "2" = "10–100 m (MEC 0.20–0.33)",
               "3" = "100–250 m (MEC 0.10)",
               "4" = "250 m–1 km (MEC −0.11)",
               "5" = "> 1 km (untested within the survey)")
)
lay <- rast(file.path(in_dir, paste0("Fig6_", cfg$name, "_layers.tif")))
s2  <- rast(file.path(in_dir, paste0("Fig6_", cfg$name, "_S2_RGB.tif")))
pts <- read.csv(file.path(in_dir, paste0("Fig6_", cfg$name, "_samples.csv")))

ov_ext <- ext(cfg$ov_ext)
zm_ext <- ext(cfg$zm_ext)
LAT0 <- cfg$lat0
deg_per_km <- 1 / (111.32 * cos(LAT0 * pi / 180))

# ---- colours ------------------------------------------------------------------
ink <- "#2b2b29"; muted <- "#6f6e69"
cls_cols <- c("1" = "#08306b", "2" = "#2563a8", "3" = "#4a90c8", "4" = "#8fc1e3",
              "5" = "#d9711c")
agc_cols <- c("#e5f5e0", "#a1d99b", "#41ab5d", "#006d2c", "#00441b")

# ---- helpers ------------------------------------------------------------------
to_ll <- function(x, res_m, method) {
  project(x, "EPSG:4326", res = res_m / 1000 * deg_per_km, method = method)
}
grey_base <- function(rgb_ll) {
  lum <- 0.299 * rgb_ll[[1]] + 0.587 * rgb_ll[[2]] + 0.114 * rgb_ll[[3]]
  g <- 0.60 + 0.40 * (lum / 255)                    # light grey basemap
  m <- as.matrix(g, wide = TRUE)
  m[is.na(m) | m <= 0.6001] <- 1                    # outside the image: white
  list(raster = as.raster(matrix(grey(m), nrow = nrow(m))), e = as.vector(ext(g)))
}
add_base <- function(b) {
  annotation_raster(b$raster, xmin = b$e[1], xmax = b$e[2], ymin = b$e[3], ymax = b$e[4],
                    interpolate = TRUE)
}
scale_bar <- function(xr, yr, km, label, right = FALSE) {
  x0 <- if (right) xr[2] - 0.05 * diff(xr) - km * deg_per_km else xr[1] + 0.05 * diff(xr)
  y0 <- yr[1] + 0.05 * diff(yr)
  x1 <- x0 + km * deg_per_km
  list(annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y0 + 0.012 * diff(yr),
                fill = ink, colour = NA),
       annotate("text", x = (x0 + x1) / 2, y = y0 + 0.03 * diff(yr), label = label,
                size = 2.5, colour = ink, vjust = 0))
}
lon_lab <- function(x) sprintf("%.2f°E", x)
lat_lab <- function(y) sprintf("%.2f°S", -y)

theme_map <- theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "#52514e", size = 6.5),
        axis.title = element_blank(),
        panel.border = element_rect(fill = NA, colour = muted, linewidth = 0.3),
        legend.title = element_text(colour = ink, size = 7.5),
        legend.text = element_text(colour = ink, size = 7),
        plot.tag = element_text(face = "bold", colour = ink, size = 9),
        plot.background = element_rect(fill = "white", colour = NA))

v_pts <- vect(pts, geom = c("lon", "lat"), crs = "EPSG:4326")

# ---- (a) overview --------------------------------------------------------------
cls_ov <- aggregate(crop(lay[["acc_class"]], ov_ext), fact = 10, fun = "min", na.rm = TRUE)
cls_ov <- as.data.frame(to_ll(cls_ov, 100, "near"), xy = TRUE)
names(cls_ov)[3] <- "cls"
cls_ov$cls <- factor(round(cls_ov$cls), levels = 1:5)
base_ov <- grey_base(to_ll(aggregate(crop(s2, ov_ext), fact = 5, fun = "mean"), 100, "bilinear"))
xr_ov <- base_ov$e[1:2]; yr_ov <- base_ov$e[3:4]

zm_ll <- as.vector(ext(project(as.polygons(zm_ext, crs = cfg$crs), "EPSG:4326")))

pa <- ggplot() + add_base(base_ov) +
  geom_raster(data = cls_ov, aes(x, y, fill = cls)) +
  geom_point(data = as.data.frame(geom(v_pts)), aes(x, y), size = 0.15, colour = "black") +
  annotate("rect", xmin = zm_ll[1], xmax = zm_ll[2], ymin = zm_ll[3], ymax = zm_ll[4],
           fill = NA, colour = ink, linewidth = 0.5) +
  scale_bar(xr_ov, yr_ov, cfg$ov_bar, paste(cfg$ov_bar, "km")) +
  scale_fill_manual(values = cls_cols, labels = cfg$cls_labs, drop = FALSE,
                    name = paste0("Expected accuracy class\n(distance to an AGC sample of ",
                                  cfg$year, ")")) +
  scale_x_continuous(labels = lon_lab, expand = c(0, 0), breaks = cfg$ov_xb) +
  scale_y_continuous(labels = lat_lab, expand = c(0, 0), breaks = cfg$ov_yb) +
  coord_fixed(ratio = 1 / cos(LAT0 * pi / 180), xlim = xr_ov, ylim = yr_ov) +
  labs(tag = "(a)") + theme_map +
  theme(legend.position = c(0.012, 0.985), legend.justification = c(0, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.9), colour = NA),
        legend.margin = margin(3, 4, 3, 4),
        legend.key.size = unit(3.2, "mm")) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(colour = NA)))

# ---- (b) and (c) zoom ------------------------------------------------------------
zm <- crop(lay, zm_ext)
agc_zm <- as.data.frame(to_ll(zm[["AGC"]], 10, "near"), xy = TRUE); names(agc_zm)[3] <- "AGC"
cls_zm <- as.data.frame(to_ll(zm[["acc_class"]], 10, "near"), xy = TRUE); names(cls_zm)[3] <- "cls"
cls_zm$cls <- factor(round(cls_zm$cls), levels = 1:5)
base_zm <- grey_base(to_ll(crop(s2, zm_ext), 20, "bilinear"))
xr_zm <- base_zm$e[1:2]; yr_zm <- base_zm$e[3:4]
pts_zm <- as.data.frame(geom(v_pts))
pts_zm <- pts_zm[pts_zm$x > xr_zm[1] & pts_zm$x < xr_zm[2] & pts_zm$y > yr_zm[1] & pts_zm$y < yr_zm[2], ]

zoom_axes <- list(
  scale_x_continuous(labels = lon_lab, expand = c(0, 0), breaks = cfg$zm_xb),
  scale_y_continuous(labels = lat_lab, expand = c(0, 0), breaks = cfg$zm_yb),
  coord_fixed(ratio = 1 / cos(LAT0 * pi / 180), xlim = xr_zm, ylim = yr_zm))

pb <- ggplot() + add_base(base_zm) +
  geom_raster(data = agc_zm, aes(x, y, fill = AGC)) +
  geom_point(data = pts_zm, aes(x, y), size = 0.06, colour = "black", alpha = 0.7) +
  scale_bar(xr_zm, yr_zm, cfg$zm_bar, paste(cfg$zm_bar, "km"), right = TRUE) +
  scale_fill_gradientn(colours = agc_cols, limits = c(10, 55), oob = scales::squish,
                       name = expression("AGC (gC m"^-2*")")) +
  zoom_axes + labs(tag = "(c)") + theme_map +
  theme(legend.position = c(0.03, 0.975), legend.justification = c(0, 1),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = scales::alpha("white", 0.9), colour = NA),
        legend.margin = margin(3, 5, 3, 5)) +
  guides(fill = guide_colourbar(barwidth = unit(28, "mm"), barheight = unit(2.2, "mm"),
                                title.position = "top"))

pc <- ggplot() + add_base(base_zm) +
  geom_raster(data = cls_zm, aes(x, y, fill = cls)) +
  scale_bar(xr_zm, yr_zm, cfg$zm_bar, paste(cfg$zm_bar, "km"), right = TRUE) +
  scale_fill_manual(values = cls_cols, drop = FALSE, guide = "none") +
  zoom_axes + labs(tag = "(b)") + theme_map

fig <- pa / (pc | pb) + plot_layout(heights = c(1.12, 1))

save_fig <- function(p, stem, w_mm, h_mm) {
  w <- w_mm / 25.4; h <- h_mm / 25.4
  ggsave(file.path(out_dir, paste0(stem, ".png")),  p, width = w, height = h, dpi = 300, bg = "white",
         device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(stem, ".tiff")), p, width = w, height = h, dpi = 300, bg = "white",
         device = ragg::agg_tiff, compression = "lzw")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")),  p, width = w, height = h, bg = "white",
         device = cairo_pdf)
  if (requireNamespace("svglite", quietly = TRUE))
    ggsave(file.path(out_dir, paste0(stem, ".svg")), p, width = w, height = h, bg = "white")
}
fig_h <- if (exists("fig_h")) fig_h else 214
save_fig(fig, paste0("Figure6_", cfg$name, "_map"), 174, fig_h)
cat("Figure 6 saved in", out_dir, "\n")
