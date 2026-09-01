# =============================================================================
# rf_applyModel_AGB.R
# STEP 3: Apply the final AGB model to each yearly dataset (2017-2024),
# producing predicted_tAGB_<year>.csv across the full study extent.
#
# Pipeline position: stage 4 of 6, after rf_model_AGB2.R. Reads
# predicted_tSPC_<year>.csv (from 03_SPC).
#
# Data availability: this script's inputs (predicted_tSPC_<year>.csv, and
# the model/predictor-structure files from rf_model_AGB2.R) are not
# included in this repository.
# =============================================================================

library(dplyr)
library(readr)
library(randomForest)
library(Metrics)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_reg.R"))  # includes apply_model_regression()

result_dir <- here("result")

# -----------------------------------------------------------------------------
# Load Final Model and Predictor Structure
# -----------------------------------------------------------------------------
model <- readRDS(file.path(result_dir, "final_rf_model_AGB.rds"))
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_AGB.rds"))
agb_predictors <- names(predictor_structure$classes)
data_years <- 2017:2024

# -----------------------------------------------------------------------------
# Loop Over Each Year
# -----------------------------------------------------------------------------
for (yr in data_years) {
  cat(sprintf("Processing year: %s\n", yr))
  
  input_path <- file.path(result_dir, sprintf("predicted_tSPC_%d.csv", yr))
  if (!file.exists(input_path)) {
    warning(sprintf("File not found: %s. Skipping year %d.\n", input_path, yr))
    next
  }
  
  df_raw <- read_csv(input_path, show_col_types = FALSE, guess_max = 100000)
  
  # Check that predictors are available
  missing_vars <- setdiff(agb_predictors, names(df_raw))
  if (length(missing_vars) > 0) {
    warning(sprintf("Missing predictors in year %d: %s", yr, paste(missing_vars, collapse = ", ")))
    next
  }
  
  # No year conversion needed
  n_before <- nrow(df_raw)
  df_filtered <- df_raw %>%
    filter(if_all(all_of(agb_predictors), ~ !is.na(.x)))
  
  n_after <- nrow(df_filtered)
  cat(sprintf("Filtered out %d rows with NA predictors (%.1f%% loss)\n",
              n_before - n_after, 100 * (n_before - n_after) / n_before))
  
  if (nrow(df_filtered) == 0) {
    warning(sprintf("No valid data left after filtering for year %d. Skipping...\n", yr))
    next
  }
  
  # Predict AGB
  df_filtered$tAGB_pred <- apply_model_regression(
    df = df_filtered,
    predictors = agb_predictors,
    model = model,
    predictor_structure = predictor_structure
  )
  
  # Save results
  out_path <- file.path(result_dir, sprintf("predicted_tAGB_%d.csv", yr))
  write_csv(df_filtered, out_path)
  cat(sprintf("Saved predicted AGB for %d -> %s\n", yr, out_path))
}