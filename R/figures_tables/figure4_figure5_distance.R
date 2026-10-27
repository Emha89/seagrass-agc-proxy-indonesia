# =============================================================================
# figure4_figure5_distance.R
# Figure 4: MEC of the GSE-only model against the median distance from the
#           test samples to the nearest training sample, for all designs.
# Figure 5: cumulative share of the mapped seagrass area within a given distance
#           of the nearest AGC sample, for all regions and for each region.
#
# Inputs (result/)  agc_block_summary_ext.csv, agc_survey_transfer_ext.csv,
#                   agc_lolocv_pointTarget.csv, agc_cross_site_pointTarget.csv,
#                   agc_design_distances_ext.csv (R/validation/),
#                   seagrass_area_by_distance_class.csv (gee/11)
# Outputs (result/paper_tables/)  Figure4_mec_vs_distance.*, Figure5_area_vs_distance.*
# =============================================================================
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

in_dir  <- if (exists("in_dir"))  in_dir  else here::here("result")
out_dir <- if (exists("out_dir")) out_dir else here::here("result", "paper_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
f <- function(name) {
  p <- list.files(in_dir, pattern = paste0(name, "$"), full.names = TRUE)
  if (length(p) == 0) stop("File not found: ", name)
  read_csv(p[1], show_col_types = FALSE)
}

blk <- f("agc_block_summary_ext.csv")
sv  <- f("agc_survey_transfer_ext.csv")
lo  <- f("agc_lolocv_pointTarget.csv")
cs  <- f("agc_cross_site_pointTarget.csv")
dd  <- f("agc_design_distances_ext.csv")
ar  <- f("seagrass_area_by_distance_class.csv")

gse <- "GSE bands only"

# ---- (a) points ----------------------------------------------------------
blk_line <- blk %>%
  filter(Scope == "All regions", Model %in% c(gse, "Mean-value baseline")) %>%
  transmute(dist = median_nn_m, MEC = MEC_mean,
            Series = ifelse(Model == gse, "Pixel and blocks (GSE-only)", "Pixel and blocks (baseline)"))

survey <- sv %>% filter(Scope == "All other data", Model == gse) %>%
  mutate(Test = paste(Region, TestYear)) %>%
  inner_join(filter(dd, Design == "Survey transfer") %>% select(Test, median_nn_m), by = "Test") %>%
  transmute(dist = median_nn_m, MEC, Design = "New survey")

lolo <- lo %>% filter(Model == gse) %>%
  inner_join(filter(dd, Design == "Unseen regions (LOLO-CV)") %>% select(Test, median_nn_m),
             by = c("TestLocation" = "Test")) %>%
  transmute(dist = median_nn_m, MEC, Design = "Unseen region (LOLO-CV)")

trans <- cs %>% filter(Type == "Transfer", Model == gse) %>%
  inner_join(filter(dd, Design == "Single-region transfer") %>% select(Test, Source, median_nn_m),
             by = c("TestSite" = "Test", "ModelSource" = "Source")) %>%
  transmute(dist = median_nn_m, MEC, Design = "Single-region transfer")

pts <- bind_rows(survey, lolo, trans) %>%
  mutate(Design = factor(Design, levels = c("New survey", "Unseen region (LOLO-CV)", "Single-region transfer")))
stopifnot(nrow(survey) == 8, nrow(lolo) == 6, nrow(trans) == 30)

# Colours (validated categorical slots 1-3; shapes as second encoding)
cols   <- c("New survey" = "#2a78d6", "Unseen region (LOLO-CV)" = "#eb6834",
            "Single-region transfer" = "#1baf7a")
shapes <- c("New survey" = 16, "Unseen region (LOLO-CV)" = 17, "Single-region transfer" = 15)
ink <- "#2b2b29"; muted <- "#8a8984"; grid <- "#e6e5e1"

x_breaks <- c(1, 10, 100, 1e3, 1e4, 1e5, 1e6)
x_labels <- c("1 m", "10 m", "100 m", "1 km", "10 km", "100 km", "1000 km")
x_lim    <- c(3, 4e6)

theme_fig <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = grid, linewidth = 0.3),
        axis.text = element_text(colour = "#52514e"),
        axis.title = element_text(colour = ink),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(colour = ink),
        plot.tag = element_text(face = "bold", colour = ink),
        plot.background = element_rect(fill = "white", colour = NA))

pa <- ggplot() +
  annotate("rect", xmin = x_lim[1], xmax = 1000, ymin = -Inf, ymax = Inf, fill = "#f1f0ec") +
  annotate("text", x = 4, y = -1.1, label = "distance range of the pixel and block tests", colour = muted, size = 2.6, hjust = 0) +
  geom_hline(yintercept = 0, colour = muted, linewidth = 0.4, linetype = "dashed") +
  geom_line(data = filter(blk_line, Series == "Pixel and blocks (baseline)"),
            aes(dist, MEC), colour = muted, linewidth = 0.6) +
  geom_point(data = filter(blk_line, Series == "Pixel and blocks (baseline)"),
             aes(dist, MEC), colour = muted, size = 1.6) +
  geom_line(data = filter(blk_line, Series == "Pixel and blocks (GSE-only)"),
            aes(dist, MEC), colour = ink, linewidth = 0.8) +
  geom_point(data = filter(blk_line, Series == "Pixel and blocks (GSE-only)"),
             aes(dist, MEC), colour = ink, size = 2.2) +
  annotate("text", x = 50, y = 0.34, label = "GSE-only, pixel and blocks", colour = ink, size = 2.8, hjust = 0) +
  annotate("text", x = 8, y = -0.075, label = "Mean-value baseline", colour = muted, size = 2.8, hjust = 0) +
  geom_point(data = pts, aes(dist, MEC, colour = Design, shape = Design), size = 2, alpha = 0.9) +
  scale_colour_manual(values = cols) + scale_shape_manual(values = shapes) +
  scale_x_log10(breaks = x_breaks, labels = x_labels, limits = x_lim) +
  scale_y_continuous(breaks = seq(-1, 0.5, 0.25)) +
  labs(x = "Distance to the nearest training sample", y = "MEC") +
  theme_fig

# ---- (b) cumulative mapped area, all regions and each region -----------------
# Figure 5 shows all regions (thick line) and each region (thin lines).
# The values of each region are given in Table S11.
cum_of <- function(d) {
  d %>% group_by(class_id, upper_m) %>% summarise(area = sum(area_ha), .groups = "drop") %>%
    arrange(class_id) %>% mutate(cum_pct = 100 * cumsum(area) / sum(area))
}
cum <- cum_of(ar) %>% filter(upper_m <= 5000)
beyond <- 100 - max(cum$cum_pct)
cum_reg <- ar %>% group_by(region) %>% group_modify(~ cum_of(.x)) %>% ungroup() %>%
  filter(upper_m <= 5000)

reg_order <- c("Ayau", "Bintan", "Karimunjawa", "Komodo", "Manado", "Rote")
reg_cols  <- c(Ayau = "#2a78d6", Bintan = "#eb6834", Karimunjawa = "#1baf7a",
               Komodo = "#eda100", Manado = "#e87ba4", Rote = "#008300")
cum_reg$region <- factor(cum_reg$region, levels = reg_order)

x_lim_b <- c(10, 5e4)   # Figure 5 is a separate figure, so its axis is zoomed to the data
pb <- ggplot() +
  annotate("rect", xmin = x_lim_b[1], xmax = 1000, ymin = -Inf, ymax = Inf, fill = "#f1f0ec") +
  annotate("text", x = 12, y = 97, label = "distance range of the pixel and block tests",
           colour = muted, size = 2.6, hjust = 0) +
  geom_line(data = cum_reg, aes(upper_m, cum_pct, colour = region), linewidth = 0.45) +
  geom_point(data = cum_reg, aes(upper_m, cum_pct, colour = region), size = 0.9) +
  geom_line(data = cum, aes(upper_m, cum_pct), colour = ink, linewidth = 1.0) +
  geom_point(data = cum, aes(upper_m, cum_pct), colour = ink, size = 2) +
  # label below the Bintan line, so that it does not cross any region line
  annotate("text", x = 1000 * 1.12, y = 2.5,
           label = sprintf("All regions: %.0f%% within 1 km",
                           cum$cum_pct[cum$upper_m == 1000]),
           colour = ink, size = 2.6, hjust = 0) +
  annotate("text", x = 5000 * 1.3, y = max(cum$cum_pct),
           label = sprintf("%.0f%% beyond 5 km", beyond), colour = ink, size = 2.8, hjust = 0) +
  scale_colour_manual(values = reg_cols, name = NULL) +
  scale_x_log10(breaks = x_breaks, labels = x_labels, limits = x_lim_b) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Distance to the nearest AGC sample", y = "Cumulative mapped area (%)") +
  theme_fig + theme(legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.8, size = 1.5)))

# Figure 4 (MEC against distance) and Figure 5 (cumulative mapped area) are saved
# as two separate figures, without panel letters.
save_fig <- function(p, stem, w_mm, h_mm) {
  w <- w_mm / 25.4; h <- h_mm / 25.4
  ggsave(file.path(out_dir, paste0(stem, ".png")),  p, width = w, height = h, dpi = 300, bg = "white",
         device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(stem, ".tiff")), p, width = w, height = h, dpi = 300, bg = "white",
         device = ragg::agg_tiff, compression = "lzw")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")),  p, width = w, height = h, bg = "white", device = cairo_pdf)
  if (requireNamespace("svglite", quietly = TRUE))
    ggsave(file.path(out_dir, paste0(stem, ".svg")), p, width = w, height = h, bg = "white")
}
save_fig(pa, "Figure4_mec_vs_distance", 174, 95)
save_fig(pb, "Figure5_area_vs_distance", 140, 80)
cat("Figures 4 and 5 saved in", out_dir, "\n")
