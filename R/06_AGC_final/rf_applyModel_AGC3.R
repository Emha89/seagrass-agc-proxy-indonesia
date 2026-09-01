# =============================================================================
# rf_applyModel_AGC3.R
# Apply the final AGC model to each yearly dataset (2017-2024), producing
# predicted_tAGC_<year>.csv -- the final pixel-level AGC estimates across
# the full study extent.
#
# Pipeline position: stage 6 of 6 (final), after rf_model_AGC3.R. Reads
# predicted_tAGB_<year>.csv (from 04_AGB) and predicted_indexC_<year>.csv
# (from 05_carbon_index).
#
# NOTE: carbon_index = carbon_index_pred (see below) is a deliberate
# rename, not a bug. The AGC model was trained on a column literally
# named "carbon_index" (the field-truth value, merged into dataPA_<year>
# during its pass-2 run). At prediction time across the full grid there
# is no field-truth value -- only carbon_index_pred, the carbon-index
# stage's own model output -- so it's renamed here to match the exact
# predictor name the AGC model expects.
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr)
library(readr)
library(randomForest)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_reg.R"))  # includes apply_model_regression()

result_dir <- here("result")

# -----------------------------------------------------------------------------
# Load Model + Predictor Structure
# -----------------------------------------------------------------------------
model <- readRDS(file.path(result_dir, "final_rf_model_AGC.rds"))
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_AGC.rds"))
agc_predictors <- names(predictor_structure$classes)
data_years <- 2017:2024

# -----------------------------------------------------------------------------
# Loop Over Each Year
# -----------------------------------------------------------------------------
for (yr in data_years) {
  
  cat(sprintf("\nProcessing year: %d\n", yr))
  
  input_path <- file.path(result_dir, sprintf("predicted_tAGB_%d.csv", yr))
  index_path <- file.path(result_dir, sprintf("predicted_indexC_%d.csv", yr))
  
  if (!file.exists(input_path)) {
    warning(sprintf("File not found: %s. Skipping.", input_path))
    next
  }
  
  if (!file.exists(index_path)) {
    warning(sprintf("Index file not found: %s. Skipping.", index_path))
    next
  }
  
  df_raw   <- read_csv(input_path, show_col_types = FALSE)
  df_index <- read_csv(index_path, show_col_types = FALSE)
  
  # -----------------------------------------------------------------------------
  # Merge carbon_index_pred + Morph3 probabilities
  # -----------------------------------------------------------------------------
  df_index_small <- df_index %>%
    select(
      gee_id,
      carbon_index_pred,
      P_mixed_short_plus_mono_short,
      P_mixed_long,
      P_mono_Ea
    )
  
  df_raw <- df_raw %>%
    left_join(df_index_small, by = "gee_id") %>%
    mutate(
      carbon_index = carbon_index_pred
    )
  
  cat("Joined carbon_index + Morph3 probabilities\n")
  
  # -----------------------------------------------------------------------------
  # Check predictors availability
  # -----------------------------------------------------------------------------
  missing_vars <- setdiff(agc_predictors, names(df_raw))
  if (length(missing_vars) > 0) {
    warning(sprintf("Missing predictors in year %d: %s",
                    yr, paste(missing_vars, collapse = ", ")))
    next
  }
  
  # -----------------------------------------------------------------------------
  # Filter NA predictors
  # -----------------------------------------------------------------------------
  df_filtered <- df_raw %>%
    filter(if_all(all_of(agc_predictors), ~ !is.na(.x)))
  
  if (nrow(df_filtered) == 0) {
    warning(sprintf("No valid rows left after NA filtering for year %d.", yr))
    next
  }
  
  # -----------------------------------------------------------------------------
  # Predict AGC
  # -----------------------------------------------------------------------------
  df_filtered$tAGC_pred <- apply_model_regression(
    df = df_filtered,
    predictors = agc_predictors,
    model = model,
    predictor_structure = predictor_structure
  )
  
  # -----------------------------------------------------------------------------
  # Save output
  # -----------------------------------------------------------------------------
  out_path <- file.path(result_dir, sprintf("predicted_tAGC_%d.csv", yr))
  write_csv(df_filtered, out_path)
  
  cat(sprintf("Saved predicted AGC -> %s (%d rows)\n",
              out_path, nrow(df_filtered)))
}

cat("\nAGC prediction completed for all years.\n")