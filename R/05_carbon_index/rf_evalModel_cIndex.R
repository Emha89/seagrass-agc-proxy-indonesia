# =============================================================================
# rf_evalModel_cIndex.R
# Evaluate Carbon Index RF predictions (2017-2024) with tables and plots.
#
# Pipeline position: stage 5 of 6, final step. Follows
# rf_applyModel_cIndex.R.
#
# Self-contained by design: defines its own metric helpers (rmse_vec,
# mae_vec, r2_vec, bias_vec) rather than reusing rf_func_reg.R, and does
# not source any shared file. r2_vec() computes R2 manually
# (1 - SS_res/SS_tot), the same approach used in
# eval_regression_performance() in rf_func_reg.R -- see the note there
# about this differing from caret::R2() used elsewhere in the codebase.
# tidyr::pivot_longer() is used via namespace only (tidyr is not
# library()-loaded above, but must be installed).
#
# Inputs:
#   1) predicted_indexC_<year>.csv  (from rf_applyModel_cIndex.R)
#      - must contain: gee_id, carbon_index_pred (and ideally lat/lon for maps)
#   2) dataCARBON_<year>.csv (truth)
#      - must contain: gee_id, carbon_index (0-1)
#
# Outputs (saved to result/):
#   - cindex_performance_by_year.csv
#   - cindex_eval_points_all_years.csv
#   - 10 diagnostic figures saved as PNG (metrics by year, scatter,
#     residuals, calibration, spatial maps if coordinates exist, etc.)
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(here)   # install.packages("here") if you don't have it yet

# -----------------------------------------------------------------------------
# SETUP
# -----------------------------------------------------------------------------
result_dir <- here("result")
data_years <- as.character(2017:2024)

out_perf_csv   <- file.path(result_dir, "cindex_performance_by_year.csv")
out_points_csv <- file.path(result_dir, "cindex_eval_points_all_years.csv")

# -----------------------------------------------------------------------------
# METRICS (robust, no external deps)
# -----------------------------------------------------------------------------
rmse_vec <- function(pred, obs) sqrt(mean((pred - obs)^2, na.rm = TRUE))
mae_vec  <- function(pred, obs) mean(abs(pred - obs), na.rm = TRUE)
r2_vec   <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  pred <- pred[ok]; obs <- obs[ok]
  if (length(obs) < 2) return(NA_real_)
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}
bias_vec <- function(pred, obs) mean(pred - obs, na.rm = TRUE)

safe_quant <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

# -----------------------------------------------------------------------------
# BLOCK 1: EVALUATE PER-YEAR METRICS
# -----------------------------------------------------------------------------
cat("Evaluating Carbon Index predictions per year...\n")

eval_points_all <- list()

perf_by_year <- purrr::map_dfr(data_years, function(yr) {
  
  pred_path  <- file.path(result_dir, paste0("predicted_indexC_", yr, ".csv"))
  truth_path <- file.path(result_dir, paste0("dataCARBON_", yr, ".csv"))
  
  if (!file.exists(pred_path)) {
    message(sprintf("Year %s skipped: missing %s", yr, basename(pred_path)))
    return(NULL)
  }
  if (!file.exists(truth_path)) {
    message(sprintf("Year %s skipped: missing %s", yr, basename(truth_path)))
    return(NULL)
  }
  
  df_pred <- read_csv(pred_path, show_col_types = FALSE, guess_max = 100000)
  df_truth <- read_csv(truth_path, show_col_types = FALSE)
  
  # Standardize minimum required columns
  if (!all(c("gee_id") %in% names(df_pred))) {
    warning(sprintf("%s: gee_id missing in prediction file.", yr))
    return(NULL)
  }
  if (!all(c("gee_id") %in% names(df_truth))) {
    warning(sprintf("%s: gee_id missing in truth file.", yr))
    return(NULL)
  }
  
  # Identify pred/truth columns
  pred_col <- if ("carbon_index_pred" %in% names(df_pred)) "carbon_index_pred" else NULL
  truth_col <- if ("carbon_index" %in% names(df_truth)) "carbon_index" else NULL
  
  if (is.null(pred_col)) {
    warning(sprintf("%s: carbon_index_pred missing in %s", yr, basename(pred_path)))
    return(NULL)
  }
  if (is.null(truth_col)) {
    warning(sprintf("%s: carbon_index missing in %s", yr, basename(truth_path)))
    return(NULL)
  }
  
  # Join
  df_eval <- inner_join(
    df_truth %>% select(gee_id, carbon_index),
    df_pred  %>% select(gee_id, carbon_index_pred, any_of(c("lat","lon","xcoord","ycoord","loc","location","area","year"))),
    by = "gee_id"
  ) %>%
    mutate(
      carbon_index = suppressWarnings(as.numeric(carbon_index)),
      carbon_index_pred = suppressWarnings(as.numeric(carbon_index_pred))
    ) %>%
    filter(
      is.finite(carbon_index), is.finite(carbon_index_pred),
      carbon_index > 0, carbon_index <= 1
    )
  
  if (nrow(df_eval) == 0) {
    message(sprintf("Year %s: 0 matched/valid rows after join+filter.", yr))
    return(NULL)
  }
  
  # Residuals
  df_eval <- df_eval %>%
    mutate(
      year_eval = yr,
      resid = carbon_index_pred - carbon_index,
      abs_err = abs(resid)
    )
  
  eval_points_all[[yr]] <<- df_eval
  
  # Metrics
  rmse <- rmse_vec(df_eval$carbon_index_pred, df_eval$carbon_index)
  mae  <- mae_vec(df_eval$carbon_index_pred, df_eval$carbon_index)
  r2   <- r2_vec(df_eval$carbon_index_pred, df_eval$carbon_index)
  bias <- bias_vec(df_eval$carbon_index_pred, df_eval$carbon_index)
  
  q <- safe_quant(df_eval$abs_err, c(0.5, 0.75, 0.95))
  
  tibble::tibble(
    year = yr,
    n = nrow(df_eval),
    rmse = round(rmse, 4),
    mae  = round(mae, 4),
    r2   = round(r2, 4),
    bias = round(bias, 4),
    abs_err_p50 = round(q[1], 4),
    abs_err_p75 = round(q[2], 4),
    abs_err_p95 = round(q[3], 4)
  )
})

if (nrow(perf_by_year) == 0) {
  stop("No evaluation results produced. Check existence/columns of predicted_indexC_YYYY and dataCARBON_YYYY.")
}

write_csv(perf_by_year, out_perf_csv)
cat(sprintf("Saved per-year performance: %s\n", out_perf_csv))
print(perf_by_year)

points_df <- bind_rows(eval_points_all)
write_csv(points_df, out_points_csv)
cat(sprintf("Saved joined evaluation points: %s (%d rows)\n", out_points_csv, nrow(points_df)))

# -----------------------------------------------------------------------------
# BLOCK 2: VISUALIZATIONS (saved to PNG)
# -----------------------------------------------------------------------------
cat("\nGenerating plots...\n")

# 2A) Metrics by year
metrics_long <- perf_by_year %>%
  select(year, rmse, mae, r2, bias) %>%
  tidyr::pivot_longer(cols = c(rmse, mae, r2, bias), names_to = "metric", values_to = "value") %>%
  mutate(year = factor(year, levels = data_years))

p_metrics <- ggplot(metrics_long, aes(x = year, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.35, size = 3) +
  labs(title = "Carbon Index Model Performance per Year", x = "Year", y = "Metric value") +
  theme_minimal(base_size = 13)

ggsave(file.path(result_dir, "fig_cindex_metrics_by_year.png"), p_metrics, width = 11, height = 5, dpi = 300)
print(p_metrics)

# 2B) Scatter (all years)
points_df <- points_df %>%
  mutate(year_eval = factor(year_eval, levels = data_years))

p_scatter <- ggplot(points_df, aes(x = carbon_index, y = carbon_index_pred, color = year_eval)) +
  geom_point(alpha = 0.45, size = 1.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Carbon Index: Observed vs Predicted (All Years)",
    x = "Observed carbon_index",
    y = "Predicted carbon_index_pred",
    color = "Year"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(result_dir, "fig_cindex_scatter_all_years.png"), p_scatter, width = 7.5, height = 6, dpi = 300)
print(p_scatter)

# 2C) Residual distribution (all years)
p_resid <- ggplot(points_df, aes(x = resid)) +
  geom_histogram(bins = 40, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Residuals (Pred - Obs) Distribution (All Years)", x = "Residual", y = "Count") +
  theme_minimal(base_size = 13)

ggsave(file.path(result_dir, "fig_cindex_residuals_all_years.png"), p_resid, width = 7.5, height = 5, dpi = 300)
print(p_resid)

# -----------------------------------------------------------------------------
# EXTRA APPENDIX PLOTS (Continuous Regression Diagnostics)
# -----------------------------------------------------------------------------

# 2D-1) Residual vs Observed (heteroscedasticity)
p_resid_vs_obs <- ggplot(points_df, aes(x = carbon_index, y = resid)) +
  geom_point(alpha = 0.35, size = 1.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residuals vs Observed Carbon Index",
    x = "Observed carbon_index",
    y = "Residual (Pred - Obs)"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(result_dir, "fig_cindex_resid_vs_obs.png"),
       p_resid_vs_obs, width = 7.5, height = 5.5, dpi = 300)
print(p_resid_vs_obs)


# 2D-2) Absolute Error Distribution per Year (Boxplot)
p_abs_err_year <- ggplot(points_df, aes(x = year_eval, y = abs_err)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.7) +
  labs(
    title = "Absolute Prediction Error per Year",
    x = "Year",
    y = "|Residual|"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(result_dir, "fig_cindex_absError_by_year.png"),
       p_abs_err_year, width = 9, height = 5.5, dpi = 300)
print(p_abs_err_year)


# 2D-3) Observed vs Predicted Density Overlay
p_density <- ggplot(points_df) +
  geom_density(aes(x = carbon_index, color = "Observed"), linewidth = 1.1) +
  geom_density(aes(x = carbon_index_pred, color = "Predicted"), linewidth = 1.1) +
  labs(
    title = "Distribution of Observed vs Predicted Carbon Index",
    x = "Carbon Index Value",
    y = "Density",
    color = ""
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(result_dir, "fig_cindex_density_obs_pred.png"),
       p_density, width = 7.5, height = 5.5, dpi = 300)
print(p_density)


# 2D-4) Scatter Faceted by Year (Cleaner than color encoding)
p_scatter_facet <- ggplot(points_df, aes(x = carbon_index, y = carbon_index_pred)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  facet_wrap(~ year_eval) +
  labs(
    title = "Observed vs Predicted Carbon Index (Per Year)",
    x = "Observed",
    y = "Predicted"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(result_dir, "fig_cindex_scatter_facet_by_year.png"),
       p_scatter_facet, width = 11, height = 7, dpi = 300)
print(p_scatter_facet)


# 2D-5) Calibration Plot (Binned Mean Prediction vs Observation)
df_cal <- points_df %>%
  mutate(bin = cut(carbon_index, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    obs_mean = mean(carbon_index, na.rm = TRUE),
    pred_mean = mean(carbon_index_pred, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_cal <- ggplot(df_cal, aes(x = obs_mean, y = pred_mean)) +
  geom_point(size = 2.5) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Calibration Plot (Binned Means)",
    x = "Mean Observed Carbon Index",
    y = "Mean Predicted Carbon Index"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(result_dir, "fig_cindex_calibration_plot.png"),
       p_cal, width = 7.5, height = 6, dpi = 300)
print(p_cal)

cat("Extra appendix diagnostic plots saved.\n")


# 2D) Spatial maps (if coordinates available)
# Prefer lat/lon; fallback xcoord/ycoord
coord_mode <- NULL
if (all(c("lon","lat") %in% names(points_df))) coord_mode <- c("lon","lat")
if (is.null(coord_mode) && all(c("xcoord","ycoord") %in% names(points_df))) coord_mode <- c("xcoord","ycoord")

if (!is.null(coord_mode)) {
  xcol <- coord_mode[1]
  ycol <- coord_mode[2]
  
  df_map <- points_df %>%
    mutate(
      x = suppressWarnings(as.numeric(.data[[xcol]])),
      y = suppressWarnings(as.numeric(.data[[ycol]]))
    ) %>%
    filter(is.finite(x), is.finite(y))
  
  if (nrow(df_map) > 0) {
    p_sp_pred <- ggplot(df_map, aes(x = x, y = y, color = carbon_index_pred)) +
      geom_point(size = 1.8, alpha = 0.85) +
      coord_fixed() +
      labs(
        title = "Spatial Distribution of Predicted Carbon Index (All Years)",
        x = xcol, y = ycol, color = "Pred"
      ) +
      theme_minimal(base_size = 13)
    
    ggsave(file.path(result_dir, "fig_cindex_spatial_pred_all_years.png"), p_sp_pred, width = 7.5, height = 6, dpi = 300)
    print(p_sp_pred)
    
    p_sp_res <- ggplot(df_map, aes(x = x, y = y, color = resid)) +
      geom_point(size = 1.8, alpha = 0.85) +
      coord_fixed() +
      labs(
        title = "Spatial Distribution of Residuals (Pred - Obs) (All Years)",
        x = xcol, y = ycol, color = "Residual"
      ) +
      theme_minimal(base_size = 13)
    
    ggsave(file.path(result_dir, "fig_cindex_spatial_resid_all_years.png"), p_sp_res, width = 7.5, height = 6, dpi = 300)
    print(p_sp_res)
  } else {
    warning("Spatial plotting skipped: no valid coordinates after coercion.")
  }
} else {
  warning("Spatial plotting skipped: no (lon,lat) or (xcoord,ycoord) columns found in joined evaluation points.")
}

cat("\nCarbon Index evaluation completed.\n")