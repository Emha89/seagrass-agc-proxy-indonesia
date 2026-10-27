# =============================================================================
# figureS3_agc_importance.R
# Figure S3: predictor importance of the AGC model of the proxy workflow.
#   (a) Share of the total importance (increase in node purity) of the 20 most
#       important predictors in the final AGC model (R/workflow/18_agc_model.R).
#   (b) Share of the proxy layers in the importance in each of the ten repeats,
#       with in-sample and with cross-fitted proxy values, in the pixel hold-out
#       and 250 m block designs (R/validation/03, PART L).
#
# Inputs (result/)  var_importance_rf_AGC.csv, agc_crossfit_iterations_ext.csv
# Outputs (result/paper_tables/)  FigureS3_importance_AGC.png / .tiff / .pdf / .svg
# For a higher resolution, set fig_dpi before running, e.g.
#   fig_dpi <- 600; source(here::here("R", "figures_tables", "figureS3_agc_importance.R"))
# =============================================================================
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

in_dir  <- if (exists("in_dir"))  in_dir  else here::here("result")
out_dir <- if (exists("out_dir")) out_dir else here::here("result", "paper_tables")
fig_dpi <- if (exists("fig_dpi")) fig_dpi else 300
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
f <- function(name) {
  p <- list.files(in_dir, pattern = paste0(name, "$"), full.names = TRUE)
  if (length(p) == 0) stop("File not found: ", name)
  read_csv(p[1], show_col_types = FALSE)
}

proxy_vars <- c(tAGB_pred = "AGB", carbon_index = "Carbon index", tSPC_pred = "Percent cover",
                PA_prob = "Presence probability",
                P_mixed_short_plus_mono_short = "Short-leaf probability",
                P_mixed_long = "Mixed-leaf probability", P_mono_Ea = "Long-leaf probability")

ink <- "#2b2b29"; muted <- "#8a8984"; grid <- "#e6e5e1"
c_in <- "#eb6834"; c_cf <- "#2a78d6"; c_gse <- "#b9b8b3"

theme_fig <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = grid, linewidth = 0.3),
        axis.text = element_text(colour = "#52514e"),
        axis.title = element_text(colour = ink),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(colour = ink),
        plot.tag = element_text(face = "bold", colour = ink),
        plot.background = element_rect(fill = "white", colour = NA))

# ---- (a) final model, top 20 predictors --------------------------------------
imp <- f("var_importance_rf_AGC.csv") %>%
  mutate(share = 100 * importance / sum(importance),
         is_proxy = variable %in% names(proxy_vars),
         label = case_when(is_proxy ~ proxy_vars[variable],
                           variable == "depth" ~ "Depth",
                           variable == "distToLand" ~ "Distance to land",
                           TRUE ~ sub("^GSE_", "GSE ", variable)))
proxy_total <- sum(imp$share[imp$is_proxy])
top <- imp %>% arrange(desc(share)) %>% slice(1:20) %>%
  mutate(label = factor(label, levels = rev(label)),
         group = ifelse(is_proxy, "Proxy layer (in-sample)", "GSE band or environmental layer"))

pa <- ggplot(top, aes(share, label, fill = group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Proxy layer (in-sample)" = c_in,
                               "GSE band or environmental layer" = c_gse)) +
  annotate("text", x = max(top$share), y = 3, hjust = 1, size = 2.8, colour = ink,
           label = sprintf("Proxy layers (orange): %.0f%% of the\ntotal importance of 73 predictors", proxy_total)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Share of the total importance (%)", y = NULL, tag = "(a)") +
  theme_fig + theme(panel.grid.major.y = element_blank(), legend.position = "none")

# ---- (b) share of proxy layers per repeat -------------------------------------
cf <- f("agc_crossfit_iterations_ext.csv") %>%
  filter(Scope == "All regions", !is.na(share_proxy_imp)) %>%
  mutate(Design = recode(Design, "Pixel-grouped" = "Pixel hold-out", "Block 250 m" = "250 m blocks"),
         Design = factor(Design, levels = c("Pixel hold-out", "250 m blocks")),
         Proxies = ifelse(Model == "Proxy workflow", "In-sample proxy values", "Cross-fitted proxy values"),
         Proxies = factor(Proxies, levels = c("In-sample proxy values", "Cross-fitted proxy values")),
         share = 100 * share_proxy_imp)
cf_mean <- cf %>% group_by(Design, Proxies) %>% summarise(share = mean(share), .groups = "drop")

pb <- ggplot(cf, aes(Design, share, colour = Proxies, shape = Proxies, group = Proxies)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.7, seed = 1),
             size = 1.6, alpha = 0.75) +
  geom_errorbar(data = cf_mean, aes(ymin = share, ymax = share), width = 0.28, linewidth = 0.8,
                position = position_dodge(width = 0.7), show.legend = FALSE) +
  geom_text(data = cf_mean, aes(y = share + 3.2, label = sprintf("%.0f%%", share)),
            position = position_dodge(width = 0.7), size = 2.8, colour = ink, show.legend = FALSE) +
  scale_colour_manual(values = c("In-sample proxy values" = c_in, "Cross-fitted proxy values" = c_cf)) +
  scale_shape_manual(values = c("In-sample proxy values" = 16, "Cross-fitted proxy values" = 17)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  labs(x = NULL, y = "Proxy share of the importance (%)", tag = "(b)") +
  theme_fig + theme(panel.grid.major.x = element_blank()) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

fig <- pa + pb + plot_layout(widths = c(1.35, 1))

w <- 174 / 25.4; h <- 95 / 25.4
ggsave(file.path(out_dir, "FigureS3_importance_AGC.png"),  fig, width = w, height = h, dpi = fig_dpi, bg = "white")
ggsave(file.path(out_dir, "FigureS3_importance_AGC.tiff"), fig, width = w, height = h, dpi = fig_dpi, bg = "white",
       compression = "lzw")
ggsave(file.path(out_dir, "FigureS3_importance_AGC.pdf"),  fig, width = w, height = h, bg = "white",
       device = cairo_pdf)
if (requireNamespace("svglite", quietly = TRUE))
  ggsave(file.path(out_dir, "FigureS3_importance_AGC.svg"), fig, width = w, height = h, bg = "white")
cat(sprintf("Proxy share in final model: %.1f%%\n", proxy_total))
print(cf_mean)
cat("Figure S3 saved in", out_dir, "\n")
