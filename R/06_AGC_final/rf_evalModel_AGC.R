# =============================================================================
# rf_evalModel_AGC.R
# Evaluate AGC predictions by year and location, using ground-truth
# AGC_pred from dataAGC_<year>.csv, and generate diagnostic visualizations.
# This is the final evaluation step of the full Study 2 pipeline.
#
# Pipeline position: stage 6 of 6, final step. Follows rf_applyModel_AGC3.R.
#
# FIX APPLIED (Block 3): the location-plot guard checked for a column
# named "location", but Block 2 above it (in this same script) groups by
# "loc", and plot_regression_performance_by_location_agc() in
# rf_func_vis.R also expects "loc" specifically. As written, the check
# would always fail and silently skip that plot even with valid location
# data present. Changed to check for "loc".
#
# Note: years covered are 2017-2023 only (no 2024), consistent with the
# AGC data-prep range (which itself follows the AGB training range).
#
# Output: agc_performance_by_year.csv, agc_performance_by_location.csv,
# written to result/.
#
# Data availability: this script's inputs (predicted_tAGC_<year>.csv and
# dataAGC_<year>.csv) are not included in this repository.
# =============================================================================

library(dplyr)
library(readr)
library(Metrics)
library(ggplot2)
library(caret)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_vis.R"))  # AGC visualization functions

rmse <- Metrics::rmse
mae  <- Metrics::mae

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
data_dir <- here("result")
years <- 2017:2023
output_csv_year <- file.path(data_dir, "agc_performance_by_year.csv")
output_csv_loc  <- file.path(data_dir, "agc_performance_by_location.csv")

# -----------------------------------------------------------------------------
# BLOCK 1: Evaluate Per-Year Metrics
# -----------------------------------------------------------------------------
results_yearly <- lapply(years, function(yr) {
  pred_path <- file.path(data_dir, sprintf("predicted_tAGC_%d.csv", yr))
  ref_path  <- file.path(data_dir, sprintf("dataAGC_%d.csv", yr))
  
  if (!file.exists(pred_path) || !file.exists(ref_path)) {
    message(sprintf("Missing file for year %d. Skipping...", yr))
    return(NULL)
  }
  
  pred_df <- read_csv(pred_path, show_col_types = FALSE) %>%
    filter(year == yr)  # only keep rows matching this year
  
  ref_df <- read_csv(ref_path, show_col_types = FALSE)
  
  df <- left_join(
    pred_df,
    ref_df %>% select(gee_id, AGC_pred),
    by = "gee_id"
  ) %>%
    filter(!is.na(AGC_pred), !is.na(tAGC_pred), AGC_pred > 0)
  
  if (nrow(df) == 0) {
    message(sprintf("No valid matched data for year %d. Skipping...", yr))
    return(NULL)
  }
  
  cat(sprintf("Year %d: %d matched samples used for evaluation.\n", yr, nrow(df)))
  
  tibble(
    year = as.character(yr),
    n = nrow(df),
    RMSE = rmse(df$AGC_pred, df$tAGC_pred),
    MAE  = mae(df$AGC_pred, df$tAGC_pred),
    R2   = R2(df$tAGC_pred, df$AGC_pred)
  )
}) %>% bind_rows()

# Save and display results
if (nrow(results_yearly) > 0) {
  write_csv(results_yearly, output_csv_year)
  cat(sprintf("Saved yearly evaluation to: %s\n", output_csv_year))
  print(results_yearly)
} else {
  warning("No yearly evaluation results produced.")
}

# -----------------------------------------------------------------------------
# BLOCK 2: Evaluate Per-Location Metrics (Across All Years)
# -----------------------------------------------------------------------------
df_all <- lapply(years, function(yr) {
  pred_path <- file.path(data_dir, sprintf("predicted_tAGC_%d.csv", yr))
  ref_path  <- file.path(data_dir, sprintf("dataAGC_%d.csv", yr))
  
  if (!file.exists(pred_path) || !file.exists(ref_path)) return(NULL)
  
  pred_df <- read_csv(pred_path, show_col_types = FALSE) %>%
    filter(year == yr) %>%
    mutate(year = as.integer(yr))
  
  ref_df <- read_csv(ref_path, show_col_types = FALSE)
  
  df_joined <- left_join(
    pred_df,
    ref_df %>% select(gee_id, AGC_pred),
    by = "gee_id"
  ) %>%
    filter(!is.na(loc), !is.na(AGC_pred), !is.na(tAGC_pred), AGC_pred > 0)
  
  if (nrow(df_joined) == 0) return(NULL)
  df_joined
}) %>% bind_rows()

if (nrow(df_all) > 0) {
  loc_results <- df_all %>%
    group_by(loc) %>%
    summarise(
      n = n(),
      RMSE = rmse(AGC_pred, tAGC_pred),
      MAE  = mae(AGC_pred, tAGC_pred),
      R2   = R2(tAGC_pred, AGC_pred),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
  
  write_csv(loc_results, output_csv_loc)
  cat(sprintf("Saved location-level evaluation to: %s\n", output_csv_loc))
  print(loc_results)
} else {
  warning("No valid data for per-location evaluation.")
}

# -----------------------------------------------------------------------------
# BLOCK 3: Visualizations (AGC)
# -----------------------------------------------------------------------------

# Per-year visualization
if (!"year" %in% colnames(df_all)) {
  warning("Column 'year' is missing in df_all. Skipping year-based plots.")
} else {
  if (exists("plot_regression_performance_by_year_agc")) {
    print(plot_regression_performance_by_year_agc(df_all))
  }
}

# Per-location visualization
if (!"loc" %in% colnames(df_all)) {
  warning("Column 'loc' is missing in df_all. Skipping location-based plots.")
} else {
  if (exists("plot_regression_performance_by_location_agc")) {
    print(plot_regression_performance_by_location_agc(df_all))
  }
}

# Error distribution & prediction visualization
if (exists("plot_regression_boxplot_actual_vs_pred_bin_agc")) {
  print(plot_regression_boxplot_actual_vs_pred_bin_agc(df_all))
}

if (exists("plot_regression_error_by_actual_bin_agc")) {
  print(plot_regression_error_by_actual_bin_agc(df_all))
}

if (exists("plot_regression_spatial_map_agc")) {
  print(plot_regression_spatial_map_agc(df_all))
}

if (exists("plot_agc_rmse_with_range_by_location")) {
  print(plot_agc_rmse_with_range_by_location(df_all))
}

cat("\nAGC Evaluation & Visualization Completed.\n")