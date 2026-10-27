# =============================================================================
# figure2_workflow.R
# Figure 2: modelling workflow for annual seagrass above-ground carbon (AGC).
# Data input -> five Random Forest proxy models -> final AGC model, inside the
# fixed seagrass extent (presence probability >= 0.8 in any year, 2017-2024).
#
# Schematic only; no data are read. Positions are in the pixel units of the
# original drawing (1089 x 410), placed on a 7.0 inch wide page, so the layout
# is the same at every size and resolution. Font sizes are given in the same
# pixel units and converted to points.
#
# fix_labels
#   FALSE  the text of the figure in the submitted manuscript
#   TRUE   corrected text: no trailing commas, "Water depth", "Seagrass
#          presence", "Above-ground", "Fixed seagrass extent", and the arrow
#          from leaf morphology to percent cover drawn dashed, as in the legend
#
# Font: Arial. On Linux, fontconfig uses Liberation Sans, which has the same
# letter widths, so the layout is the same on every system.
#
# Outputs (result/paper_tables/)
#   Figure2_workflow.png / .tiff  raster, 300 dpi by default
#   Figure2_workflow.pdf / .svg   vector, editable text
# For a higher resolution, set fig_dpi before running, e.g.
#   fig_dpi <- 600; source(here::here("R", "figures_tables", "figure2_workflow.R"))
# =============================================================================
library(grid)

out_dir    <- if (exists("out_dir")) out_dir else here::here("result", "paper_tables")
fig_dpi    <- if (exists("fig_dpi")) fig_dpi else 300
fix_labels <- if (exists("fix_labels")) fix_labels else FALSE
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

FONT <- "Arial"
PX_W <- 1089; PX_H <- 410                       # drawing units (pixels of the original)
W <- 7.0; H <- W * PX_H / PX_W                  # page size (inches)
px_pt  <- 72 * W / PX_W                         # one drawing unit in points
fs     <- function(px) px * px_pt               # font size in points
lw_px  <- function(px) px * px_pt * 96 / 72     # line width in drawing units to grid lwd

# colours of the original drawing
GREY  <- "#666666"; BLACK <- "#000000"; TXT <- "#000000"
INPUT <- "#EEEEEE"; NONSG <- "#FAD7AC"; SG <- "#B9E0A5"
AGC_FILL <- "#009900"; AGC_EDGE <- "#007A00"; NOTE <- "#4D4D4D"

# -----------------------------------------------------------------------------
# Labels. Each line is a list of (text, face) pieces; faces are plain, bold,
# italic and bold.italic.
# -----------------------------------------------------------------------------
p <- function(...) list(...)
seg <- function(text, face = "plain") c(text, face)

if (!fix_labels) {
  L <- list(
    legend1 = ": Generate raster",
    legend2 = ": Additional variable",
    input = list(
      p(seg("1. Time-series Google's Satellite Embedding ("), seg("predictor variable", "italic"), seg("),")),
      p(seg("2. Field data set ("), seg("response variables", "italic"), seg("),")),
      p(seg("3. Bathymetry ("), seg("additional predictor variable", "italic"), seg("),")),
      p(seg("4. Distance to land ("), seg("additional predictor variable", "italic"), seg("),"))),
    m1 = list(p(seg("(1) Present - absent", "bold")),
              p(seg("(Random Forest probability)"))),
    m2 = list(p(seg("(2) Leaf morphology (", "bold"), seg("short-leaf,", "bold.italic")),
              p(seg("long-leaf, mixed-leaf", "bold.italic"), seg(")", "bold")),
              p(seg("(Random Forest classification)"))),
    m3 = list(p(seg("(3) Percent Cover (SPC)", "bold")),
              p(seg("(Random Forest regression)"))),
    m4 = list(p(seg("(4) Above ground biomass", "bold")),
              p(seg("(AGB)", "bold")),
              p(seg("(Random Forest regression)"))),
    m5 = list(p(seg("(5) Carbon index (C Index)", "bold")),
              p(seg("(Random Forest regression)"))),
    agc = "Seagrass Above ground carbon (AGC)",
    nonsg = c("Non", "seagrass", "area"),
    extent = "Consistent seagrass area (presence probability ≥ 0.8 in any year from 2017 to 2024)",
    dashed_morph_arrow = FALSE)
} else {
  L <- list(
    legend1 = ": Predicted layer",
    legend2 = ": Additional predictor",
    input = list(
      p(seg("1. Annual Google Satellite Embeddings, GSE ("), seg("predictors", "italic"), seg(")")),
      p(seg("2. Field data ("), seg("response variables", "italic"), seg(")")),
      p(seg("3. Water depth ("), seg("additional predictor", "italic"), seg(")")),
      p(seg("4. Distance to land ("), seg("additional predictor", "italic"), seg(")"))),
    m1 = list(p(seg("(1) Seagrass presence", "bold")),
              p(seg("(Random Forest probability)"))),
    m2 = list(p(seg("(2) Leaf morphology (", "bold"), seg("short leaf,", "bold.italic")),
              p(seg("long leaf, mixed leaf", "bold.italic"), seg(")", "bold")),
              p(seg("(Random Forest classification)"))),
    m3 = list(p(seg("(3) Percent cover (SPC)", "bold")),
              p(seg("(Random Forest regression)"))),
    m4 = list(p(seg("(4) Above-ground biomass", "bold")),
              p(seg("(AGB)", "bold")),
              p(seg("(Random Forest regression)"))),
    m5 = list(p(seg("(5) Carbon index", "bold")),
              p(seg("(Random Forest regression)"))),
    agc = "Seagrass above-ground carbon (AGC)",
    nonsg = c("Non-", "seagrass", "area"),
    extent = "Fixed seagrass extent (presence probability ≥ 0.8 in any year from 2017 to 2024)",
    dashed_morph_arrow = TRUE)
}

# -----------------------------------------------------------------------------
# Drawing helpers (all positions in drawing units; y runs downwards)
# -----------------------------------------------------------------------------
rect <- function(x0, y0, x1, y1, fill, col = NA, lwd = 1, r = 0) {
  # centre in drawing units, size in inches (the y axis of the viewport runs downwards)
  xc <- unit((x0 + x1) / 2, "native"); yc <- unit((y0 + y1) / 2, "native")
  w  <- unit((x1 - x0) * W / PX_W, "inches"); h <- unit((y1 - y0) * H / PX_H, "inches")
  gp <- gpar(fill = fill, col = col, lwd = lw_px(lwd))
  if (r > 0) {
    grid.roundrect(xc, yc, w, h, r = unit(r * W / PX_W, "inches"), gp = gp)
  } else {
    grid.rect(xc, yc, w, h, gp = gp)
  }
}

# arrowhead in the style of the original ("classic": filled, with a notch);
# (x, y) is the tip, (dx, dy) the direction of the line
head <- function(x, y, dx, dy, col, len = 15, half = 6.5, notch = 0.30) {
  n <- sqrt(dx^2 + dy^2); ux <- dx / n; uy <- dy / n
  px <- -uy; py <- ux
  bx <- x - ux * len; by <- y - uy * len
  nx <- x - ux * len * (1 - notch); ny <- y - uy * len * (1 - notch)
  grid.polygon(c(x, bx + px * half, nx, bx - px * half),
               c(y, by + py * half, ny, by - py * half),
               default.units = "native", gp = gpar(fill = col, col = col, lwd = 0.5))
}

# polyline with an arrowhead at the last point; the line stops inside the head
path <- function(xs, ys, col, lwd = 3, arrow = TRUE, lty = "solid", len = 15) {
  n <- length(xs)
  dx <- xs[n] - xs[n - 1]; dy <- ys[n] - ys[n - 1]
  xe <- xs; ye <- ys
  if (arrow) {
    d <- sqrt(dx^2 + dy^2)
    xe[n] <- xs[n] - dx / d * len * 0.7; ye[n] <- ys[n] - dy / d * len * 0.7
  }
  grid.lines(xe, ye, default.units = "native",
             gp = gpar(col = col, lwd = lw_px(lwd), lty = lty, lineend = "butt", linejoin = "mitre"))
  if (arrow) head(xs[n], ys[n], dx, dy, col, len = len)
}

# one line of mixed-face text, centred at x (just = "centre") or starting at x
rich <- function(x, y, pieces, size, col = TXT, just = "centre") {
  gps <- lapply(pieces, function(s) gpar(fontfamily = FONT, fontsize = fs(size), fontface = s[2], col = col))
  ws  <- mapply(function(s, g) convertWidth(grobWidth(textGrob(s[1], gp = g)), "native", valueOnly = TRUE),
                pieces, gps)
  x0 <- if (just == "centre") x - sum(ws) / 2 else x
  for (k in seq_along(pieces)) {
    grid.text(pieces[[k]][1], x0 + sum(ws[seq_len(k - 1)]), y, just = c("left", "centre"),
              default.units = "native", gp = gps[[k]])
  }
}

txt <- function(x, y, label, size, face = "plain", col = TXT, just = "centre") {
  grid.text(label, x, y, just = c(just, "centre"), default.units = "native",
            gp = gpar(fontfamily = FONT, fontsize = fs(size), fontface = face, col = col))
}

# box of a model, with its lines of text centred
model_box <- function(x0, x1, y0, y1, lines, ys, border = NA, size = 12.1) {
  rect(x0, y0, x1, y1, fill = "white", col = border, lwd = 1)
  for (k in seq_along(lines)) rich((x0 + x1) / 2, ys[k], lines[[k]], size)
}

# -----------------------------------------------------------------------------
# The figure
# -----------------------------------------------------------------------------
draw <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, PX_W), yscale = c(PX_H, 0)))

  # areas
  rect(18, 200, 138, PX_H, fill = NONSG)
  rect(158, 200, PX_W, PX_H, fill = SG)
  for (k in seq_along(L$nonsg)) txt(78, c(277, 307, 337)[k], L$nonsg[k], 25.3, "italic", NOTE)
  txt(181, 386, L$extent, 19.8, "italic", NOTE, just = "left")

  # legend
  path(c(8, 48), c(20, 20), BLACK)
  txt(53, 20, L$legend1, 11.9, just = "left")
  path(c(8, 48), c(55, 55), BLACK, lty = "22")
  txt(56, 60, L$legend2, 11.9, just = "left")

  # data input
  rect(383, 9, 763, 110, fill = INPUT, r = 12)
  txt(389, 19, "Data input", 14.3, "bold", just = "left")
  for (k in seq_along(L$input)) rich(415, c(46.5, 61.5, 76, 91)[k], L$input[[k]], 12, just = "left")

  # model boxes
  model_box(68, 228, 120, 160, L$m1, c(133, 148), border = BLACK)
  model_box(258, 448, 221, 270, L$m2, c(231, 245.5, 260))
  model_box(478, 668, 221, 270, L$m3, c(237.5, 253))
  model_box(678, 868, 221, 270, L$m4, c(231, 245.5, 260))
  model_box(878, 1068, 221, 270, L$m5, c(237.5, 253))

  # AGC model
  rect(458, 320, 688, 360, fill = AGC_FILL, col = AGC_EDGE, lwd = 1, r = 5)
  txt(573, 332.5, L$agc, 12.5, "bold", "white")
  txt(573, 347, "(Random Forest regression)", 12, "plain", TXT)

  # grey arrows: inputs to the models
  path(c(383, 148, 148), c(85, 85, 120), GREY)
  path(c(148, 148), c(160, 180), GREY, arrow = FALSE)
  path(c(78, 229), c(180, 180), GREY, arrow = FALSE)
  path(c(78, 78), c(180, 200), GREY)
  path(c(229, 229), c(180, 200), GREY)
  path(c(478, 478, 353, 353), c(110, 165, 165, 221), GREY)
  path(c(573, 573), c(110, 221), GREY)
  path(c(668, 668, 773, 773), c(110, 165, 165, 221), GREY)
  path(c(763, 973, 973), c(85, 85, 221), GREY)
  path(c(448, 478), c(245, 245), GREY, lty = if (L$dashed_morph_arrow) "22" else "solid")

  # black arrows: predicted layers to the AGC model
  path(c(228, 248, 248, 458), c(150, 150, 341, 341), BLACK)
  path(c(353, 353, 516, 516), c(270, 290, 290, 320), BLACK)
  path(c(573, 573), c(270, 320), BLACK)
  path(c(773, 773, 630, 630), c(270, 290, 290, 320), BLACK)
  path(c(973, 973, 688), c(270, 331, 331), BLACK)
  path(c(763, 1079, 1079, 688), c(35, 35, 350, 350), BLACK)

  popViewport()
}

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
stem <- file.path(out_dir, "Figure2_workflow")
ragg::agg_png(paste0(stem, ".png"), width = W, height = H, units = "in", res = fig_dpi,
              background = "white"); draw(); invisible(dev.off())
ragg::agg_tiff(paste0(stem, ".tiff"), width = W, height = H, units = "in", res = fig_dpi,
               background = "white", compression = "lzw"); draw(); invisible(dev.off())
cairo_pdf(paste0(stem, ".pdf"), width = W, height = H, bg = "white"); draw(); invisible(dev.off())
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(paste0(stem, ".svg"), width = W, height = H, bg = "white"); draw(); invisible(dev.off())
}
cat("Figure 2 saved in", out_dir, "\n")
