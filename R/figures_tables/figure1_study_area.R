# =============================================================================
# figure1_study_area.R
# Figure 1: the six study regions on a map of Indonesia, with six insets that
# show the field samples (seagrass present or absent) on satellite imagery.
#
# The layout follows the figure in the manuscript, which was made in QGIS:
# the main map fills the page, the insets sit in two rows at the top and the
# bottom, and the legend is at the top right. Positions are in the pixel units
# of that figure (1378 x 974) on a 6.26 x 4.43 inch page, so the layout is the
# same at every resolution.
#
# Basemaps are downloaded with the maptiles package: Esri World Light Gray
# Canvas for the main map and Esri World Imagery for the insets. The QGIS
# figure used Bing VirtualEarth imagery, which maptiles does not provide, so
# the imagery in the insets looks slightly different. The basemap name is
# printed in the legend; give the imagery source in the figure caption.
#
# Inputs (data/)
#   field_points.csv      xcoord, ycoord (decimal degrees), loc (region), PA (0/1)
#   eez_indonesia.shp     Indonesian Exclusive Economic Zone (Marine Regions,
#                         World EEZ v12, https://www.marineregions.org)
# Outputs (result/paper_tables/)  Figure1_study_area.png / .tiff / .pdf
#
# To match the QGIS figure, set the centre of each inset (centre_lon,
# centre_lat) in the table `insets` below. When they are NA, the inset is
# centred on the field samples of the region.
# For a higher resolution, set fig_dpi before running, e.g.
#   fig_dpi <- 600; source(here::here("R", "figures_tables", "figure1_study_area.R"))
# =============================================================================
library(sf)
library(dplyr)
library(readr)
library(terra)
library(maptiles)
library(grid)

out_dir    <- if (exists("out_dir"))    out_dir    else here::here("result", "paper_tables")
data_dir   <- if (exists("data_dir"))   data_dir   else here::here("data")
fig_dpi    <- if (exists("fig_dpi"))    fig_dpi    else 300
fix_labels <- if (exists("fix_labels")) fix_labels else FALSE
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

FONT <- "Tahoma"
PX_W <- 1378; PX_H <- 974                       # layout units (pixels of the QGIS figure)
W <- PX_W / 220; H <- PX_H / 220                # page size (inches)
fs  <- function(px) px * 72 / 220               # font size in layout pixels to points
inch <- function(px) unit(px / 220, "inches")

# colours of the QGIS figure
SEA <- "#D0CFD4"; ABSENT <- "#E31A1C"; PRESENT <- "#33A02C"
MARKER <- "#E56769"; EEZ_COL <- "#6E6E73"; LABEL <- "#2B2B2B"; COUNTRY <- "#808080"

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
main_lon     <- c(89.3, 141.7)                  # west and east edge of the main map
main_lat_top <- 16.6                            # north edge; the south edge follows from the page shape
main_zoom    <- 6                               # tile zoom of the main map
main_tiles   <- "Esri.WorldGrayCanvas"
inset_tiles  <- "Esri.WorldImagery"
basemap_name <- "Esri World Imagery"

# One row per inset. col/row: position in the layout. width_km: width of the
# inset. label_*, bar_*: position of the region name and of the 1 km scale bar,
# in pixels from the top-left corner of the inset. ink: colour of the text and
# the scale bar. zoom: tile zoom.
insets <- data.frame(
  loc        = c("Bintan", "Manado", "Ayau", "Karimunjawa", "Komodo", "Rote"),
  col        = c(1, 2, 3, 1, 2, 3),
  row        = c(1, 1, 1, 2, 2, 2),
  width_km   = c(5.35, 2.42, 2.40, 5.56, 2.62, 2.67),
  centre_lon = NA_real_,
  centre_lat = NA_real_,
  label_x    = c(157, 142, 9, 5, 180, 14),
  label_y    = c(192, 192, 192, 20, 19, 167),
  bar_x      = c(181, 116, 11, 15, 123, 11),
  bar_y      = c(218, 218, 218, 41, 53, 189),
  ink        = c("white", "white", "black", "white", "white", "white"),
  zoom       = c(15, 16, 16, 15, 16, 16),
  stringsAsFactors = FALSE)
inset_x0 <- c(0, 285, 570); inset_w <- 278      # columns of the insets
inset_y0 <- c(0, 695);      inset_h <- 279      # rows of the insets

if (!fix_labels) {
  leg <- list(eez = "Indonesia Exclusive Economic Zone", title = "Seagrass present-absent",
              absent = "false", present = "true", basemap = "basemap")
} else {
  leg <- list(eez = "Indonesian Exclusive Economic Zone", title = "Field samples",
              absent = "Seagrass absent", present = "Seagrass present", basemap = "Basemap of the insets")
}

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
pts <- read_csv(file.path(data_dir, "field_points.csv"), show_col_types = FALSE) %>%
  filter(!is.na(xcoord), !is.na(ycoord), !is.na(PA), loc %in% insets$loc) %>%
  mutate(PA = as.integer(as.logical(PA))) %>%
  st_as_sf(coords = c("xcoord", "ycoord"), crs = 4326) %>%
  st_transform(3857)
xy <- st_coordinates(pts)
pts$X <- xy[, 1]; pts$Y <- xy[, 2]

eez <- st_read(file.path(data_dir, "eez_indonesia.shp"), quiet = TRUE) %>% st_transform(3857)

merc <- function(lon, lat) {
  p <- st_coordinates(st_transform(st_sfc(st_point(c(lon, lat)), crs = 4326), 3857))
  c(p[1], p[2])
}

# main map window (EPSG:3857), with the shape of the page
mx <- c(merc(main_lon[1], 0)[1], merc(main_lon[2], 0)[1])
m_per_px <- diff(mx) / PX_W
my <- merc(0, main_lat_top)[2] - c(PX_H * m_per_px, 0)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
# tiles of a window (EPSG:3857) as a raster image and its extent
get_image <- function(xlim, ylim, provider, zoom) {
  win <- st_as_sfc(st_bbox(c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
                           crs = st_crs(3857)))
  r <- get_tiles(win, provider = provider, zoom = zoom, crop = TRUE)
  a <- as.array(r)[, , 1:3]
  a[is.na(a)] <- 255
  img <- rgb(a[, , 1], a[, , 2], a[, , 3], maxColorValue = 255)
  dim(img) <- dim(a)[1:2]
  list(img = as.raster(img), ext = as.vector(ext(r)))
}

# viewport of a map window placed at (x0, y0, w, h) layout pixels from the top left
map_vp <- function(x0, y0, w, h, xlim, ylim) {
  viewport(x = inch(x0), y = unit(H, "inches") - inch(y0), width = inch(w), height = inch(h),
           just = c("left", "top"), xscale = xlim, yscale = ylim, clip = "on")
}

draw_image <- function(im) {
  e <- im$ext
  grid.raster(im$img, x = e[1], y = e[3], width = e[2] - e[1], height = e[4] - e[3],
              just = c("left", "bottom"), default.units = "native", interpolate = TRUE)
}

# text in layout pixels (y from the top), with an optional white halo
txt <- function(x, y, label, px, face = "plain", col = LABEL, just = "left", halo = NA) {
  pos <- function(dx, dy) list(x = inch(x + dx), y = unit(H, "inches") - inch(y + dy))
  gp  <- function(colour) gpar(fontfamily = FONT, fontsize = fs(px), fontface = face, col = colour)
  if (!is.na(halo)) {
    for (a in seq(0, 2 * pi, length.out = 17)[-17]) {
      p <- pos(2.2 * cos(a), 2.2 * sin(a))
      grid.text(label, p$x, p$y, just = c(just, "centre"), gp = gp(halo))
    }
  }
  p <- pos(0, 0)
  grid.text(label, p$x, p$y, just = c(just, "centre"), gp = gp(col))
}

px_rect <- function(x0, y0, x1, y1, fill, col = NA, lty = "solid", lwd = 1) {
  grid.rect(inch(x0), unit(H, "inches") - inch(y0), inch(x1 - x0), inch(y1 - y0),
            just = c("left", "top"), gp = gpar(fill = fill, col = col, lty = lty, lwd = lwd))
}
px_line <- function(xs, ys, col = "black", lwd = 1.6, lty = "solid", arrow = NULL) {
  grid.lines(inch(xs), unit(H, "inches") - inch(ys), arrow = arrow,
             gp = gpar(col = col, fill = col, lwd = lwd, lty = lty))
}

# -----------------------------------------------------------------------------
# Inset windows
# -----------------------------------------------------------------------------
for (k in seq_len(nrow(insets))) {
  if (is.na(insets$centre_lon[k])) {
    p <- pts[pts$loc == insets$loc[k], ]
    b <- st_bbox(st_transform(p, 4326))
    insets$centre_lon[k] <- mean(b[c("xmin", "xmax")])
    insets$centre_lat[k] <- mean(b[c("ymin", "ymax")])
  }
}
insets <- insets %>%
  rowwise() %>%
  mutate(cx = merc(centre_lon, centre_lat)[1], cy = merc(centre_lon, centre_lat)[2],
         scale = 1 / cos(centre_lat * pi / 180),            # 3857 metres per ground metre
         half_w = width_km * 1000 / 2 * scale,
         half_h = half_w * inset_h / inset_w,
         bar_px = inset_w / width_km) %>%                   # 1 km in layout pixels
  ungroup()

cat("Inset centres (lon, lat) used:\n")
print(as.data.frame(insets[, c("loc", "centre_lon", "centre_lat", "width_km")]))

# -----------------------------------------------------------------------------
# Download the basemaps once
# -----------------------------------------------------------------------------
main_img  <- get_image(mx, my, main_tiles, main_zoom)
inset_img <- lapply(seq_len(nrow(insets)), function(k) {
  with(insets[k, ], get_image(cx + c(-half_w, half_w), cy + c(-half_h, half_h), inset_tiles, zoom))
})

# region markers: mean position of the field samples of each region
markers <- pts %>% st_drop_geometry() %>% group_by(loc) %>%
  summarise(X = mean(X), Y = mean(Y), .groups = "drop") %>%
  mutate(px = (X - mx[1]) / m_per_px, py = (my[2] - Y) / m_per_px)

# -----------------------------------------------------------------------------
# The figure
# -----------------------------------------------------------------------------
draw <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = SEA, col = NA))

  # main map
  pushViewport(map_vp(0, 0, PX_W, PX_H, mx, my))
  draw_image(main_img)
  for (g in st_geometry(st_cast(st_boundary(eez), "MULTILINESTRING"))) {
    for (ln in g) grid.lines(ln[, 1], ln[, 2], default.units = "native",
                             gp = gpar(col = EEZ_COL, lwd = 0.5))
  }
  popViewport()
  ll <- function(lon, lat) { m <- merc(lon, lat); c((m[1] - mx[1]) / m_per_px, (my[2] - m[2]) / m_per_px) }
  a <- ll(113.6, -0.45); txt(a[1], a[2], "Indonesia", 27, "bold.italic", COUNTRY, just = "centre")
  a <- ll(132.2, -17.4); txt(a[1], a[2], "Australia", 27, "bold.italic", COUNTRY, just = "centre")

  # region markers and names
  for (k in seq_len(nrow(markers))) {
    grid.circle(inch(markers$px[k]), unit(H, "inches") - inch(markers$py[k]), inch(11),
                gp = gpar(fill = adjustcolor(MARKER, 0.55), col = MARKER, lwd = 2.2))
    txt(markers$px[k] + 9, markers$py[k], markers$loc[k], 19, "bold", LABEL, halo = "white")
  }

  # north arrow and 500 km scale bar
  px_line(c(27, 27), c(532, 392))
  star <- seq(pi / 2, by = pi / 5, length.out = 10)          # white star at the top
  rad  <- rep(c(10, 4), 5)
  grid.polygon(inch(27 + rad * cos(star)), unit(H, "inches") - inch(381 - rad * sin(star)),
               gp = gpar(fill = "white", col = NA))
  px_line(c(27, 68), c(532, 418), arrow = arrow(length = inch(9), angle = 22, type = "open"))
  arc <- seq(pi / 2, atan2(532 - 418, 68 - 27), length.out = 20)   # dotted arc from north to the arrow
  px_line(27 + 124 * cos(arc), 532 - 124 * sin(arc), lwd = 1.2, lty = "22")
  bar <- 500e3 / (m_per_px * cos(4 * pi / 180))          # 500 km at about 4 degrees S
  px_rect(22, 538, 22 + bar, 552, fill = "black")
  txt(22, 579, "0", 20, just = "centre")
  txt(22 + bar, 579, "500 km", 20, just = "centre")

  # insets
  for (k in seq_len(nrow(insets))) {
    i <- insets[k, ]
    x0 <- inset_x0[i$col]; y0 <- inset_y0[i$row]
    pushViewport(map_vp(x0, y0, inset_w, inset_h, i$cx + c(-i$half_w, i$half_w), i$cy + c(-i$half_h, i$half_h)))
    draw_image(inset_img[[k]])
    p <- pts[pts$loc == i$loc, ]
    for (v in c(0L, 1L)) {
      q <- p[p$PA == v, ]
      if (nrow(q) > 0) grid.points(q$X, q$Y, pch = 16, size = inch(9 / 1.05), default.units = "native",
                                   gp = gpar(col = if (v == 1) PRESENT else ABSENT))
    }
    popViewport()
    txt(x0 + i$label_x, y0 + i$label_y, i$loc, 22, "bold", i$ink)
    px_rect(x0 + i$bar_x, y0 + i$bar_y - 7, x0 + i$bar_x + i$bar_px, y0 + i$bar_y + 7, fill = i$ink)
    txt(x0 + i$bar_x, y0 + i$bar_y + 34, "0", 20, col = i$ink, just = "centre")
    txt(x0 + i$bar_x + i$bar_px - 5, y0 + i$bar_y + 34, "1 km", 20, col = i$ink)
  }

  # legend
  px_rect(972, 26, 1005, 43, fill = NA, col = "#7F7F7F", lty = "22", lwd = 0.8)
  txt(1015, 35, leg$eez, 22)
  txt(972, 75, leg$title, 27, "bold")
  grid.circle(inch(988), unit(H, "inches") - inch(105), inch(4.5), gp = gpar(fill = ABSENT, col = NA))
  txt(1013, 105, leg$absent, 22)
  grid.circle(inch(988), unit(H, "inches") - inch(135), inch(4.5), gp = gpar(fill = PRESENT, col = NA))
  txt(1013, 135, leg$present, 22)
  txt(972, 173, leg$basemap, 22)
  txt(972, 211, basemap_name, 27, "bold")
}

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
stem <- file.path(out_dir, "Figure1_study_area")
ragg::agg_png(paste0(stem, ".png"), width = W, height = H, units = "in", res = fig_dpi,
              background = "white"); draw(); invisible(dev.off())
ragg::agg_tiff(paste0(stem, ".tiff"), width = W, height = H, units = "in", res = fig_dpi,
               background = "white", compression = "lzw"); draw(); invisible(dev.off())
cairo_pdf(paste0(stem, ".pdf"), width = W, height = H, bg = "white"); draw(); invisible(dev.off())
cat("Figure 1 saved in", out_dir, "\n")
