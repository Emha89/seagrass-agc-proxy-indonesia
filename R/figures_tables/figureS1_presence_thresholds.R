# =============================================================================
# figureS1_presence_thresholds.R
# Figure S1: probability threshold that gave the highest F1-score of the
# presence model in each year (2017-2024). For each year, the labels of all
# field samples (2018-2023) were compared with the probabilities predicted from
# the GSE layer of that year. The values are printed by
# R/workflow/04_presence_thresholds.R (plot_pa_prob_distribution_optimal(),
# method = "f1") and typed in below.
#
# Output (result/paper_tables/)  FigureS1_presence_thresholds.png / .tiff / .pdf / .svg
# For a higher resolution, set fig_dpi before running, e.g.
#   fig_dpi <- 600; source(here::here("R", "figures_tables", "figureS1_presence_thresholds.R"))
# =============================================================================
library(ggplot2)

out_dir <- if (exists("out_dir")) out_dir else here::here("result", "paper_tables")
fig_dpi <- if (exists("fig_dpi")) fig_dpi else 300
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

thr <- data.frame(
  year = 2017:2024,
  threshold = c(0.46, 0.68, 0.59, 0.38, 0.21, 0.30, 0.52, 0.56)
)
lab <- thr[thr$threshold %in% range(thr$threshold), ]

ink <- "#2b2b29"; muted <- "#8a8984"; grid <- "#e6e5e1"
p <- ggplot(thr, aes(year, threshold)) +
  geom_line(colour = ink, linewidth = 0.6) +
  geom_point(colour = ink, size = 2.2) +
  geom_text(data = lab, aes(label = sprintf("%.2f", threshold)),
            vjust = ifelse(lab$threshold == max(thr$threshold), -1.0, 1.9),
            colour = ink, size = 2.9) +
  scale_x_continuous(breaks = 2017:2024) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Year of the GSE layer", y = "Threshold with the highest F1-score") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = grid, linewidth = 0.3),
        axis.text = element_text(colour = "#52514e"),
        axis.title = element_text(colour = ink),
        plot.background = element_rect(fill = "white", colour = NA))

w <- 120 / 25.4; h <- 70 / 25.4
ggsave(file.path(out_dir, "FigureS1_presence_thresholds.png"),  p, width = w, height = h, dpi = fig_dpi, bg = "white")
ggsave(file.path(out_dir, "FigureS1_presence_thresholds.tiff"), p, width = w, height = h, dpi = fig_dpi, bg = "white",
       compression = "lzw")
ggsave(file.path(out_dir, "FigureS1_presence_thresholds.pdf"),  p, width = w, height = h, bg = "white",
       device = cairo_pdf)
if (requireNamespace("svglite", quietly = TRUE))
  ggsave(file.path(out_dir, "FigureS1_presence_thresholds.svg"), p, width = w, height = h, bg = "white")
cat("Figure S1 saved in", out_dir, "\n")
