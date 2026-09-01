# =============================================================================
# rf_evalModel_AGB.R
# Evaluate AGB predictions by year and location, using ground-truth AGB_pred
# from dataAGB_<year>.csv, and generate diagnostic visualizations.
#
# Pipeline position: stage 4 of 6, final step. Follows rf_applyModel_AGB.R.
#
# Output: agb_performance_by_year.csv, agb_performance_by_location.csv,
# written to result/.
#
# Data availability: this script's inputs (predicted_tAGB_<year>.csv and
# dataAGB_<year>.csv) are not included in this repository.
# =============================================================================

library(dplyr)
library(readr)
library(Metrics)
library(ggplot2)
library(caret)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_vis.R"))  # AGB visualization functions

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
data_dir <- here("result")
years <- 2017:2024
output_csv_year <- file.path(data_dir, "agb_performance_by_year.csv")
output_csv_loc  <- file.path(data_dir, "agb_performance_by_location.csv")

# -----------------------------------------------------------------------------
# BLOCK 1: Evaluate Per-Year Metrics
# -----------------------------------------------------------------------------
results_yearly <- lapply(years, function(yr) {
  pred_path <- file.path(data_dir, sprintf("predicted_tAGB_%d.csv", yr))
  ref_path  <- file.path(data_dir, sprintf("dataAGB_%d.csv", yr))
  
  # Check that both files exist
  if (!file.exists(pred_path) || !file.exists(ref_path)) {
    message(sprintf("Missing file for year %d. Skipping...", yr))
    return(NULL)
  }
  
  # --- Load data ---
  pred_df <- read_csv(pred_path, show_col_types = FALSE) %>%
    filter(year == yr)  # only keep rows matching this year
  
  ref_df  <- read_csv(ref_path, show_col_types = FALSE)
  
  # --- Join & filter valid samples ---
  df <- left_join(
    pred_df,
    ref_df %>% select(gee_id, AGB_pred),
    by = "gee_id"
  ) %>%
    filter(!is.na(AGB_pred), !is.na(tAGB_pred), AGB_pred > 0)
  
  if (nrow(df) == 0) {
    message(sprintf("No valid matched data for year %d. Skipping...", yr))
    return(NULL)
  }
  
  cat(sprintf("Year %d: %d matched samples used for evaluation.\n", yr, nrow(df)))
  
  # --- Metrics ---
  tibble(
    year = as.character(yr),
    n = nrow(df),
    RMSE = Metrics::rmse(df$AGB_pred, df$tAGB_pred),
    MAE  = Metrics::mae(df$AGB_pred, df$tAGB_pred),
    R2   = caret::R2(df$tAGB_pred, df$AGB_pred)
  )
}) %>% bind_rows()

# --- Save and display results ---
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
  pred_path <- file.path(data_dir, sprintf("predicted_tAGB_%d.csv", yr))
  ref_path  <- file.path(data_dir, sprintf("dataAGB_%d.csv", yr))
  
  if (!file.exists(pred_path) || !file.exists(ref_path)) return(NULL)
  
  pred_df <- read_csv(pred_path, show_col_types = FALSE) %>%
    filter(year == yr) %>%
    mutate(year = as.integer(yr))  # make sure the year column is kept
  
  ref_df <- read_csv(ref_path, show_col_types = FALSE)
  
  df_joined <- left_join(
    pred_df,
    ref_df %>% select(gee_id, AGB_pred),
    by = "gee_id"
  ) %>%
    filter(!is.na(loc), !is.na(AGB_pred), !is.na(tAGB_pred), AGB_pred > 0)
  
  if (nrow(df_joined) == 0) return(NULL)
  df_joined
}) %>% bind_rows()

if (nrow(df_all) > 0) {
  loc_results <- df_all %>%
    group_by(loc) %>%
    summarise(
      n = n(),
      RMSE = Metrics::rmse(AGB_pred, tAGB_pred),
      MAE  = Metrics::mae(AGB_pred, tAGB_pred),
      R2   = caret::R2(tAGB_pred, AGB_pred),
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
# BLOCK 3: Visualizations (AGB)
# -----------------------------------------------------------------------------

# Make sure df_all has valid 'year' and 'loc' columns
if (!"year" %in% colnames(df_all)) {
  warning("Column 'year' is missing in df_all. Skipping year-based plots.")
} else {
  if (exists("plot_regression_performance_by_year_agb")) {
    print(plot_regression_performance_by_year_agb(df_all))
  }
}

if (!"loc" %in% colnames(df_all)) {
  warning("Column 'loc' is missing in df_all. Skipping location-based plots.")
} else {
  if (exists("plot_regression_performance_by_location_agb")) {
    print(plot_regression_performance_by_location_agb(df_all))
  }
}

if (exists("plot_regression_boxplot_actual_vs_pred_bin_agb")) {
  print(plot_regression_boxplot_actual_vs_pred_bin_agb(df_all))
}

if (exists("plot_regression_error_by_actual_bin_agb")) {
  print(plot_regression_error_by_actual_bin_agb(df_all))
}

if (exists("plot_regression_spatial_map_agb")) {
  print(plot_regression_spatial_map_agb(df_all))
}

if (exists("plot_agb_rmse_with_range_by_location")) {
  print(plot_agb_rmse_with_range_by_location(df_all))
}

cat("\nAGB Evaluation & Visualization Completed.\n")