# =============================================================================
# rf_applyModel_CINDEX.R
# Apply the final RF regression model (Carbon Index) to annual data,
# producing predicted_indexC_<year>.csv across the full study extent.
#
# Pipeline position: stage 5 of 6, after rf_model_combined_CINDEX.R.
# Reads predicted_PA_<year>.csv (from 01_occurrence_PA).
#
# Data availability: this script's inputs (predicted_PA_<year>.csv, and
# the model/predictor-structure files from rf_model_combined_CINDEX.R)
# are not included in this repository.
# =============================================================================

library(dplyr)
library(readr)
library(randomForest)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_reg.R"))  # includes apply_model_regression()

cat("Loading model and predictor structure...\n")

result_dir <- here("result")

model <- readRDS(file.path(result_dir, "final_rf_model_CINDEX.rds"))
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_CINDEX.rds"))
cindex_predictors <- names(predictor_structure$classes)
data_years <- as.character(2017:2024)

for (yr in data_years) {
  cat(sprintf("\nRunning prediction for year: %s\n", yr))
  
  in_path <- file.path(result_dir, paste0("predicted_PA_", yr, ".csv"))
  out_path <- file.path(result_dir, paste0("predicted_indexC_", yr, ".csv"))
  
  if (!file.exists(in_path)) {
    warning(sprintf("File not found: %s. Skipping year %s.\n", in_path, yr))
    next
  }
  
  df <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)
  
  # Ensure all required predictors are available
  missing <- setdiff(cindex_predictors, colnames(df))
  if (length(missing) > 0) {
    warning(sprintf("Missing predictors in %s: %s", yr, paste(missing, collapse = ", ")))
    next
  }
  
  valid_idx <- df %>%
    select(all_of(cindex_predictors)) %>%
    complete.cases()
  
  df_valid <- df[valid_idx, ]
  
  # Predict Carbon Index
  df$carbon_index_pred <- NA_real_
  df$carbon_index_pred[valid_idx] <- apply_model_regression(
    df = df_valid,
    predictors = cindex_predictors,
    model = model,
    predictor_structure = predictor_structure
  )
  
  # Clamp predictions to 0-1
  n_below_0 <- sum(df$carbon_index_pred < 0, na.rm = TRUE)
  n_above_1 <- sum(df$carbon_index_pred > 1, na.rm = TRUE)
  if (n_below_0 + n_above_1 > 0) {
    cat(sprintf("Clamped %d predictions below 0 and %d above 1 to [0-1].\n",
                n_below_0, n_above_1))
  }
  df$carbon_index_pred <- pmin(pmax(df$carbon_index_pred, 0), 1)
  
  # Save output
  write_csv(df, out_path)
  cat(sprintf("Saved: %s (%d rows)\n", out_path, nrow(df)))
}