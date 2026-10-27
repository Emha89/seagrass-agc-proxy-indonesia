# =============================================================================
# figure3_validation_designs.R
# Figure 3: leak-free validation designs, drawn as Input / Process / Output.
# Five designs ordered by distance (pixel, spatial blocks, new survey, unseen
# region, single-region transfer), with the median distances of Table 3, the
# cross-fitted proxy run and the sensitivity run with proxy models trained on
# all samples. The quadrat sketch shows one transect with 3-5 quadrats per
# 10 m pixel.
#
# Schematic only; no data are read. Region sample sizes are typed in from
# Table 2. All positions are in inches on a 7.0 x 6.45 inch page (full page
# width), and font sizes are in points.
#
# Font: Arial. On Linux, fontconfig uses Liberation Sans, which has the same
# letter widths, so the layout is the same on every system.
#
# Outputs (result/paper_tables/)
#   Figure3_validation_design.png / .tiff  raster, 300 dpi by default
#   Figure3_validation_design.pdf / .svg   vector, editable text
# For a higher resolution, set fig_dpi before running, e.g.
#   fig_dpi <- 600; source(here::here("R", "figures_tables", "figure3_validation_designs.R"))
# =============================================================================
library(grid)

out_dir <- if (exists("out_dir")) out_dir else here::here("result", "paper_tables")
fig_dpi <- if (exists("fig_dpi")) fig_dpi else 300
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

FONT <- "Arial"
W <- 7.0; H <- 6.45                              # page size (inches)

# colours (Okabe-Ito training and test colours, colour-blind safe)
TRAIN <- "#0072B2"; TEST <- "#E69F00"
EDGE <- "#7F8794"; FILL <- "#F3F4F6"; BAND <- "#3B4A5A"; SHARED <- "#E6EEF6"
TXT <- "#1F2328"; GREY <- "#5A6270"; LINE <- "#B8BEC7"

regions <- c("Ayau", "Bintan", "Karimunjawa", "Komodo", "Manado", "Rote")
n_reg   <- c(114, 59, 263, 1398, 292, 1550)

# -----------------------------------------------------------------------------
# Drawing helpers. Every element is stored with a layer number (z) and drawn
# from the lowest to the highest layer: boxes 1-2, cells 3, arrows 4, text 5.
# -----------------------------------------------------------------------------
grobs <- list(); zs <- numeric(0)
add <- function(g, z) {
  grobs[[length(grobs) + 1]] <<- g
  zs <<- c(zs, z)
}

pt  <- function(size) size / 72                 # points to inches
lwd <- function(lw) lw * 96 / 72                 # line width in points to grid lwd
EM_LINE <- (1491 + 431) / 2048                   # font ascender + descender (share of the font size)
DASH <- "31"                                     # dashed line pattern (in 1/96 inch)

# rounded box; (x, y) is the lower-left corner
box <- function(x, y, w, h, fc = FILL, ec = EDGE, lw = 0.7, r = 0.04, z = 1, ls = "-") {
  add(roundrectGrob(x, y, w, h, just = c("left", "bottom"), default.units = "native",
                    r = unit(r, "inches"),
                    gp = gpar(fill = fc, col = ec, lwd = lwd(lw),
                              lty = if (ls == "--") DASH else "solid")), z)
}

# plain rectangle (bars and cells)
rect <- function(x, y, w, h, fc, ec = NA, lw = 0.6, z = 3) {
  add(rectGrob(x, y, w, h, just = c("left", "bottom"), default.units = "native",
               gp = gpar(fill = fc, col = ec, lwd = lwd(lw))), z)
}

# one training or test cell
cell <- function(x, y, w, h, test) rect(x, y, w, h, fc = if (test) TEST else TRAIN, ec = "white")

# arrow with a filled head (head length 0.4 x ms, half width 0.2 x ms, in points).
# The line is shortened a little, so that the outline of the head ends at (x1, y1).
arrow_line <- function(x0, y0, x1, y1, color = GREY, lw = 0.9, ms = 7) {
  angle <- atan(0.2 / 0.4)
  back  <- pt(0.5 * lw / sin(angle))
  len   <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  x1 <- x1 - back * (x1 - x0) / len
  y1 <- y1 - back * (y1 - y0) / len
  add(segmentsGrob(x0, y0, x1, y1, default.units = "native",
                   arrow = arrow(angle = angle * 180 / pi, length = unit(sqrt(0.4^2 + 0.2^2) * ms, "points"),
                                 type = "closed"),
                   gp = gpar(col = color, fill = color, lwd = lwd(lw), linejoin = "mitre",
                             lineend = "butt")), 4)
}

# text; "\n" starts a new line. ha = left / center / right, va = center / top.
# Each line takes a fixed height (linespacing x font height), and the letters
# of each line are centred vertically within it.
txt <- function(x, y, s, size = 7, weight = "normal", color = TXT, ha = "left",
                va = "center", style = "normal", rot = 0, linespacing = 1.2, z = 5) {
  face <- if (weight == "bold") 2 else if (style == "italic") 3 else 1
  gp <- gpar(fontsize = size, fontface = face, col = color, fontfamily = FONT)
  L <- linespacing * EM_LINE * pt(size)          # height of one line (inches)
  hj <- c(left = 0, center = 0.5, right = 1)[[ha]]
  # baseline offset from the middle of a line: half of (ascent - descent) of its letters
  mid_to_base <- function(g) 0.5 * (grobAscent(g) - grobDescent(g))
  if (rot == 90) {                               # one line, centred on (x, y)
    g <- textGrob(s, gp = gp)
    add(textGrob(s, unit(x, "native") + mid_to_base(g), unit(y, "native"), rot = 90,
                 hjust = 0.5, vjust = 0, gp = gp), z)
    return(invisible(NULL))
  }
  lines <- strsplit(s, "\n", fixed = TRUE)[[1]]
  n <- length(lines)
  top <- if (va == "top") y else y + n * L / 2
  for (i in seq_len(n)) {
    g <- textGrob(lines[i], gp = gp)
    add(textGrob(lines[i], unit(x, "native"),
                 unit(top - (i - 0.5) * L, "native") - mid_to_base(g),
                 hjust = hj, vjust = 0, gp = gp), z)
  }
}

# -----------------------------------------------------------------------------
# Bands: INPUT, PROCESS, OUTPUT
# -----------------------------------------------------------------------------
X0 <- 0.42; X1 <- 6.95
IN <- c(5.10, 6.40); PR <- c(1.00, 4.98); OUT <- c(0.05, 0.90)
bands <- list(INPUT = IN, PROCESS = PR, OUTPUT = OUT)
for (lab in names(bands)) {
  y0 <- bands[[lab]][1]; y1 <- bands[[lab]][2]
  box(0.05, y0, 0.27, y1 - y0, fc = BAND, ec = BAND, r = 0.03)
  txt(0.185, (y0 + y1) / 2, lab, size = 8, weight = "bold", color = "white",
      ha = "center", rot = 90)
}

# -----------------------------------------------------------------------------
# INPUT
# -----------------------------------------------------------------------------
iy0 <- IN[1]; iy1 <- IN[2]

# field samples and sample size of each region
box(X0, iy0, 2.70, iy1 - iy0)
txt(X0 + 0.08, iy1 - 0.13, "Field samples", size = 7.5, weight = "bold")
txt(X0 + 0.08, iy1 - 0.29, "3,676 AGC quadrats in 1,767 pixels, one record per sample",
    size = 6.5, color = GREY)
bx0 <- X0 + 0.80; bmax <- 1.45
for (i in seq_along(regions)) {
  yy <- iy1 - 0.47 - (i - 1) * 0.135
  bl <- bmax * n_reg[i] / max(n_reg)
  txt(bx0 - 0.05, yy, regions[i], size = 6.5, ha = "right")
  rect(bx0, yy - 0.045, bl, 0.09, fc = "#9AA3AF")
  txt(bx0 + bl + 0.04, yy, format(n_reg[i], big.mark = ","), size = 6.5)
}
txt(X0 + 2.62, iy1 - 0.47, "Proxy stages also use\nthe other samples\nin Table 2", size = 6,
    color = GREY, ha = "right", va = "top", linespacing = 1.25)

# photo quadrats along one diagonal transect inside 10-m pixels
x2 <- X0 + 2.83
box(x2, iy0, 1.62, iy1 - iy0)
txt(x2 + 0.08, iy1 - 0.13, "Quadrats and pixels", size = 7.5, weight = "bold")
cs <- 0.20; gx <- x2 + 0.10; gy <- iy0 + 0.27
for (c in 0:2) for (r in 0:2) rect(gx + c * cs, gy + r * cs, cs, cs, fc = "white", ec = EDGE)
add(segmentsGrob(gx + 0.08 * cs, gy + 0.08 * cs, gx + 2.92 * cs, gy + 2.92 * cs,
                 default.units = "native",
                 gp = gpar(col = LINE, lwd = lwd(0.8), lineend = "square")), 3.5)
for (c in 0:2) for (f in c(0.2, 0.4, 0.6, 0.8)) {    # 4 photo quadrats in each pixel
  add(circleGrob(gx + (c + f) * cs, gy + (c + f) * cs, r = unit(0.015, "inches"),
                 default.units = "native", gp = gpar(fill = TXT, col = NA)), 4)
}
txt(gx + 1.5 * cs, gy - 0.10, "10-m pixels", size = 6, color = GREY, ha = "center")
txt(gx + 3 * cs + 0.08, gy + 1.5 * cs,
    paste("3–5 photo quadrats", "per pixel along a", "transect; 71% of",
          "quadrats share a", "pixel and identical", "GSE values", sep = "\n"),
    size = 6.1, linespacing = 1.25)

# predictors and target
x3 <- x2 + 1.75
box(x3, iy0, X1 - x3, iy1 - iy0)
txt(x3 + 0.08, iy1 - 0.13, "Predictors", size = 7.5, weight = "bold")
txt(x3 + 0.08, iy1 - 0.42, "64 GSE bands (annual, 10 m)\nWater depth\nDistance to land",
    size = 6.5, linespacing = 1.35)
txt(x3 + 0.08, iy0 + 0.40, "Target", size = 7.5, weight = "bold")
# "AGC point estimate (gC m-2)" with a raised "-2"
local({
  size <- 6.5; gp <- gpar(fontsize = size, col = TXT, fontfamily = FONT)
  gp_sup <- gpar(fontsize = size * 0.7, col = TXT, fontfamily = FONT)
  raise <- unit(0.36 * size, "points")           # height of the raised "-2"
  t1 <- "AGC point estimate (gC m"; t2 <- "−2"; t3 <- ")"
  a <- max(grobAscent(textGrob(t1, gp = gp)), raise + grobAscent(textGrob(t2, gp = gp_sup)))
  d <- grobDescent(textGrob(t1, gp = gp))
  yb <- unit(iy0 + 0.19, "native") - 0.5 * (a - d)   # baseline
  x1 <- unit(x3 + 0.08, "native")
  g1 <- textGrob(t1, x1, yb, hjust = 0, vjust = 0, gp = gp)
  g2 <- textGrob(t2, x1 + grobWidth(g1) + unit(0.3, "points"), yb + raise,
                 hjust = 0, vjust = 0, gp = gp_sup)
  g3 <- textGrob(t3, x1 + grobWidth(g1) + grobWidth(g2) + unit(0.6, "points"), yb,
                 hjust = 0, vjust = 0, gp = gp)
  add(g1, 5); add(g2, 5); add(g3, 5)
})

# -----------------------------------------------------------------------------
# PROCESS
# -----------------------------------------------------------------------------
py1 <- PR[2]
txt(X0, py1 - 0.13, "Split of field samples", size = 7.5, weight = "bold")
lx <- X0 + 1.55
for (k in 1:2) {
  rect(lx, py1 - 0.18, 0.16, 0.10, fc = c(TRAIN, TEST)[k])
  txt(lx + 0.21, py1 - 0.13, c("Training", "Test")[k], size = 6.5)
  lx <- lx + 0.80
}
txt(X0 + 4.55, py1 - 0.13, "ordered by distance to training data", size = 6.3, color = GREY,
    style = "italic", ha = "right")

barx0 <- X0 + 0.07; barx1 <- X0 + 3.80
bw <- (barx1 - barx0) / 6
for (i in seq_along(regions)) {
  txt(barx0 + (i - 0.5) * bw, py1 - 0.34, regions[i], size = 6.3, ha = "center", color = GREY)
}

lanes <- data.frame(
  tag   = c("A", "B", "C", "D", "E"),
  title = c("Pixel hold-out", "Spatial blocks", "New survey", "Unseen region (LOLO-CV)",
            "Single-region transfer"),
  dist  = c("median distance 6 m", "median distance 43–663 m",
            "median distance 57 m–23 km", "median distance 411–976 km",
            "median distance 411–2,933 km"),
  rep   = c("×10 repeats", "×10 repeats", "×8 surveys", "×6 folds",
            "×30 pairs"),
  note  = c("30% of pixels in each region held out; quadrats of one pixel kept together",
            "30% of 100, 250, 500 or 1000 m blocks in each region held out",
            "One survey year of one region held out (one shown); trained on all other data",
            "Each region held out once (one fold shown)",
            "Each region used once as source; each other region scored separately"),
  stringsAsFactors = FALSE)
# example pattern of held-out cells (0-based positions within each region)
pix_test <- list(c(2, 6, 8), c(1, 4, 7), c(0, 5, 8), c(3, 6, 9), c(1, 2, 7), c(4, 5, 8))
blk_test <- c(2, 0, 3, 1, 2, 0)

lane_w <- 4.55; lane_h <- 0.58; gap_l <- 0.05
lane_top <- py1 - 0.43
lane_y <- numeric(0)
for (k in seq_len(nrow(lanes))) {
  tag <- lanes$tag[k]
  top <- lane_top - (k - 1) * (lane_h + gap_l)
  box(X0, top - lane_h, lane_w, lane_h, fc = "white", ec = EDGE, lw = 0.6)
  txt(X0 + 0.07, top - 0.11, tag, size = 8, weight = "bold")
  txt(X0 + 0.23, top - 0.11, lanes$title[k], size = 7, weight = "bold")
  txt(X0 + 2.05, top - 0.11, lanes$dist[k], size = 6.5, color = GREY, style = "italic")
  by <- top - 0.35; bh <- 0.15
  for (i in 1:6) {
    x <- barx0 + (i - 1) * bw
    if (tag == "A") {                             # 10 pixels per region, 3 held out
      pw <- bw / 10
      for (j in 0:9) cell(x + j * pw, by, pw, bh, j %in% pix_test[[i]])
    } else if (tag == "B") {                      # 4 blocks per region, 1 held out
      qw <- bw / 4
      for (j in 0:3) cell(x + j * qw, by, qw, bh, j == blk_test[i])
    } else if (tag == "C") {                      # surveys of Karimunjawa, one held out
      if (i == 3) {
        sw <- bw / 3
        for (j in 0:2) cell(x + j * sw, by, sw, bh, j == 0)
        txt(x + sw / 2, by + bh / 2, "2018", size = 5, color = TXT, ha = "center")
      } else {
        cell(x, by, bw, bh, FALSE)
      }
    } else {
      is_test <- if (tag == "D") i == 5 else i != 4
      cell(x, by, bw, bh, is_test)
    }
  }
  txt(barx1 + 0.08, by + bh / 2, lanes$rep[k], size = 7, weight = "bold")
  txt(barx0, by - 0.10, lanes$note[k], size = 6, color = GREY)
  lane_y <- c(lane_y, by + bh / 2)
}
lanes_bot <- lane_top - 5 * lane_h - 4 * gap_l

# shared box: steps in every repeat, fold, survey and pair
sx0 <- X0 + lane_w + 0.15; sx1 <- X1
s_top <- py1 - 0.02; s_bot <- lanes_bot
box(sx0, s_bot, sx1 - sx0, s_top - s_bot, fc = SHARED, ec = EDGE)
cx <- (sx0 + sx1) / 2
txt(cx, s_top - 0.15, "In every repeat, fold,\nsurvey and pair", size = 7, weight = "bold",
    ha = "center", linespacing = 1.15)
steps <- list(
  list("Training samples only", NA, "-"),
  list("Refit the 5 proxy models", "presence, morphology, cover,\nAGB, carbon index", "-"),
  list("Cross-fitted proxy values",
       "extra run in A and B only:\nout-of-fold values for training samples", "--"),
  list("Fit the AGC models", NA, "-"),
  list("Predict the test samples", NA, "-"))
sy <- s_top - 0.40
for (j in seq_along(steps)) {
  s <- steps[[j]][[1]]; sub <- steps[[j]][[2]]; ls <- steps[[j]][[3]]
  hh <- if (!is.na(sub)) 0.34 else 0.19
  box(sx0 + 0.10, sy - hh, sx1 - sx0 - 0.20, hh, fc = "white", ec = EDGE, lw = 0.6, r = 0.03,
      z = 2, ls = ls)
  if (!is.na(sub)) {
    txt(cx, sy - 0.10, s, size = 6.5, weight = "bold", ha = "center")
    txt(cx, sy - 0.24, sub, size = 6, color = GREY, ha = "center", linespacing = 1.15)
  } else {
    txt(cx, sy - hh / 2, s, size = 6.5, weight = "bold", ha = "center")
  }
  if (j < length(steps)) arrow_line(cx, sy - hh, cx, sy - hh - 0.08, ms = 5)
  sy <- sy - hh - 0.08
}
txt(sx0 + 0.12, sy - 0.02, "Models compared on\nthe same test samples",
    size = 6.3, weight = "bold", va = "top", linespacing = 1.2)
txt(sx0 + 0.12, sy - 0.32,
    paste("• Proxy workflow", "• Proxy workflow, cross-fitted", "   (A and B)",
          "• GSE-only model", "• Mean-value baseline", "   (training mean)", sep = "\n"),
    size = 6.2, va = "top", linespacing = 1.3)
for (ym in lane_y) arrow_line(X0 + lane_w, ym, sx0, ym, ms = 6)   # lanes -> shared box

# sensitivity run
d_top <- lanes_bot - 0.08; d_bot <- PR[1] + 0.05
box(X0, d_bot, X1 - X0, d_top - d_bot, fc = "white", ec = EDGE, lw = 0.6, ls = "--")
txt(X0 + 0.07, (d_top + d_bot) / 2, "Sensitivity", size = 7, weight = "bold")
txt(X0 + 0.85, (d_top + d_bot) / 2,
    paste("Design A repeated with the proxy models trained on all samples, including the test samples,",
          "to show how much refitting inside each split matters"), size = 6.3, color = GREY)

# -----------------------------------------------------------------------------
# OUTPUT
# -----------------------------------------------------------------------------
outs <- list(
  c("Metrics", "RMSE, MAE, R², MEC and\nregional mean bias (%),\nscored against the field\npoint estimates"),
  c("A and B", "Mean and 2.5th–97.5th\npercentiles over 10 repeats,\noverall and per region;\npaired test of RMSE"),
  c("C, D and E", "Mean and range over\n8 surveys, 6 regions\nor 30 region pairs,\nequal weight"),
  c("Distance", "Median distance to the\nnearest training sample;\nmapped area by distance\nto the nearest AGC sample"),
  c("Sensitivity", "Apparent accuracy with\nproxy models trained\non all samples"))
oy0 <- OUT[1]; oy1 <- OUT[2]
og <- 0.125
ow <- (X1 - X0 - 4 * og) / 5
for (k in seq_along(outs)) {
  x <- X0 + (k - 1) * (ow + og)
  box(x, oy0, ow, oy1 - oy0, fc = if (k == 1) SHARED else FILL)
  txt(x + 0.08, oy1 - 0.12, outs[[k]][1], size = 7.5, weight = "bold")
  txt(x + 0.08, oy1 - 0.25, outs[[k]][2], size = 6.2, va = "top", linespacing = 1.3)
  arrow_line(x + ow / 2, PR[1], x + ow / 2, oy1)
}

# input -> process arrows
arrow_line(X0 + 1.35, iy0, X0 + 1.35, py1)
arrow_line(x2 + 0.81, iy0, x2 + 0.81, py1)
arrow_line(x3 + (X1 - x3) / 2, iy0, x3 + (X1 - x3) / 2, s_top)

# -----------------------------------------------------------------------------
# Draw and save
# -----------------------------------------------------------------------------
figure3 <- gTree(children = do.call(gList, grobs[order(zs)]),
                 vp = viewport(xscale = c(0, W), yscale = c(0, H)))
draw <- function() { grid.newpage(); grid.draw(figure3) }
stem <- file.path(out_dir, "Figure3_validation_design")

ragg::agg_png(paste0(stem, ".png"), width = W, height = H, units = "in", res = fig_dpi,
              background = "white"); draw(); invisible(dev.off())
ragg::agg_tiff(paste0(stem, ".tiff"), width = W, height = H, units = "in", res = fig_dpi,
               background = "white", compression = "lzw"); draw(); invisible(dev.off())
cairo_pdf(paste0(stem, ".pdf"), width = W, height = H, bg = "white"); draw(); invisible(dev.off())
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(paste0(stem, ".svg"), width = W, height = H, bg = "white"); draw(); invisible(dev.off())
}
cat("Figure 3 saved in", out_dir, "\n")
