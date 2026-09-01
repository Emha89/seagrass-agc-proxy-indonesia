# =============================================================================
# rf_evalModel_SPC.R
# Evaluate SPC (percent cover) regression predictions by year and location,
# and generate diagnostic visualizations.
#
# Pipeline position: stage 3 of 6, final step. Follows rf_applyModel_SPC.R.
#
# Output: pct_performance_by_year.csv, pct_performance_by_location.csv,
# written to result/.
#
# Data availability: this script's input (predicted_tSPC_<year>.csv,
# produced by rf_applyModel_SPC.R) is not included in this repository.
# =============================================================================

library(dplyr)
library(readr)
library(Metrics)
library(ggplot2)
library(yardstick)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_vis.R"))  # visualizations

data_dir <- here("result")
data_years <- as.character(2017:2024)
out_csv_year <- file.path(data_dir, "pct_performance_by_year.csv")
out_csv_loc <- file.path(data_dir, "pct_performance_by_location.csv")

# -----------------------------------------------------------------------------
# Block 1: Evaluate Per-Year Metrics (2017-2024)
# -----------------------------------------------------------------------------
results_yearly <- lapply(data_years, function(yr) {
  path <- file.path(data_dir, paste0("predicted_tSPC_", yr, ".csv"))
  if (!file.exists(path)) {
    message(sprintf("File not found for year %s. Skipped.", yr))
    return(NULL)
  }
  
  df <- read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(tSPC), !is.na(tSPC_pred), year == yr, tSPC > 0, tSPC <= 100)
  
  if (nrow(df) == 0) {
    message(sprintf("No valid GT data for year %s. Skipped.", yr))
    return(NULL)
  }
  
  tibble(
    year = yr,
    n = nrow(df),
    RMSE = rmse_vec(truth = df$tSPC, estimate = df$tSPC_pred),
    MAE  = mae_vec(truth = df$tSPC, estimate = df$tSPC_pred),
    R2 = rsq_vec(truth = df$tSPC, estimate = df$tSPC_pred)
  )
}) %>% bind_rows()

write_csv(results_yearly, out_csv_year)
cat(sprintf("Saved yearly performance to: %s\n", out_csv_year))
print(results_yearly)

# -----------------------------------------------------------------------------
# Block 2: Evaluate Per-Location Metrics
# -----------------------------------------------------------------------------
data_all <- lapply(data_years, function(yr) {
  path <- file.path(data_dir, paste0("predicted_tSPC_", yr, ".csv"))
  if (!file.exists(path)) return(NULL)
  
  read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(tSPC), !is.na(tSPC_pred), year == yr)
}) %>% bind_rows()

if (nrow(data_all) > 0) {
  location_metrics <- data_all %>%
    group_by(loc) %>%
    summarise(
      n = n(),
      RMSE = rmse_vec(truth = tSPC, estimate = tSPC_pred),
      MAE  = mae_vec(truth = tSPC, estimate = tSPC_pred),
      R2   = rsq_vec(truth = tSPC, estimate = tSPC_pred),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
  
  write_csv(location_metrics, out_csv_loc)
  cat("Saved location-level evaluation to: pct_performance_by_location.csv\n")
  print(location_metrics)
} else {
  warning("No data available for per-location evaluation.")
}

# -----------------------------------------------------------------------------
# Block 3: Visualizations
# -----------------------------------------------------------------------------
print(plot_regression_performance_by_year(data_dir))
print(plot_regression_performance_by_location(data_dir))
print(plot_regression_boxplot_actual_vs_pred_bin(data_dir))
print(plot_regression_error_by_actual_bin(data_dir))
print(plot_regression_spatial_map(data_dir))

cat("\nAll evaluations and visualizations completed.\n")