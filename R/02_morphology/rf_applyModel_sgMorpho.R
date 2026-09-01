# =============================================================================
# rf_applyModel_MORPH3.R
# Apply the final RF morphology classification model (3-class) to annual
# data, producing predicted morphology probabilities across the full
# study extent (not just training points).
#
# Pipeline position: stage 2 of 6, after rf_model_combined_MORPH3.R.
# Reads predicted_PA_<year>.csv (from 01_occurrence_PA) and appends P_*
# probability columns without touching any existing column -- so the
# lon/lat fix applied in rf_applyModel_PA.R carries through unchanged.
#
# Output: predicted_MORPH3_probs_<year>.csv, one per year, written to
# result/. This is the exact file dataPrep_func_PA.R and
# rf_applyModel_PA.R are waiting on to backfill their own outputs -- see
# the pipeline ordering note in dataPrep_func_PA.R. Also used as an
# additional predictor input for the AGC regression stage.
#
# Data availability: this script's inputs (predicted_PA_<year>.csv, and
# the model/predictor-structure files from rf_model_combined_MORPH3.R)
# are not included in this repository.
# =============================================================================

library(dplyr)
library(readr)
library(randomForest)
library(here)   # install.packages("here") if you don't have it yet

cat("Loading morphology model and predictor structure...\n")

result_dir <- here("result")

# -----------------------------------------------------------------------------
# Load trained morphology RF model
# -----------------------------------------------------------------------------
model <- readRDS(file.path(result_dir, "final_rf_model_MORPH3.rds"))

# Predictor structure (same approach as Carbon Index)
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_MORPH3.rds"))

# Predictor names required by the model
morph_predictors <- names(predictor_structure$classes)

# Years to apply the model
data_years <- as.character(2017:2024)

# -----------------------------------------------------------------------------
# Apply morphology model year-by-year
# -----------------------------------------------------------------------------
for (yr in data_years) {
  
  cat(sprintf("\nRunning morphology probability prediction for year: %s\n", yr))
  
  # Input file (same source as Carbon Index apply)
  in_path <- file.path(result_dir, paste0("predicted_PA_", yr, ".csv"))
  
  # Output file (new: morphology probabilities)
  out_path <- file.path(result_dir, paste0("predicted_MORPH3_probs_", yr, ".csv"))
  
  # Skip missing years
  if (!file.exists(in_path)) {
    warning(sprintf("File not found: %s. Skipping year %s.\n", in_path, yr))
    next
  }
  
  # Load yearly dataset
  df <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)
  
  # ---------------------------------------------------------------------------
  # Ensure all required predictors exist
  # ---------------------------------------------------------------------------
  missing <- setdiff(morph_predictors, colnames(df))
  if (length(missing) > 0) {
    warning(sprintf("Missing predictors in %s: %s",
                    yr, paste(missing, collapse = ", ")))
    next
  }
  
  # Complete-case filtering
  valid_idx <- df %>%
    select(all_of(morph_predictors)) %>%
    complete.cases()
  
  df_valid <- df[valid_idx, ]
  
  # ---------------------------------------------------------------------------
  # Predict morphology probabilities (multi-class RF)
  # ---------------------------------------------------------------------------
  cat("Predicting morphology probabilities...\n")
  
  prob_matrix <- predict(
    model,
    newdata = df_valid[, morph_predictors],
    type = "prob"
  )
  
  prob_df <- as.data.frame(prob_matrix)
  
  # Rename columns to paper-ready probability predictors
  colnames(prob_df) <- paste0("P_", colnames(prob_df))
  
  # Initialize probability columns in full dataframe
  for (cc in colnames(prob_df)) {
    df[[cc]] <- NA_real_
  }
  
  # Fill only valid rows
  df[valid_idx, colnames(prob_df)] <- prob_df
  
  # ---------------------------------------------------------------------------
  # Save output file
  # ---------------------------------------------------------------------------
  write_csv(df, out_path)
  
  cat(sprintf("Saved: %s (%d rows)\n", out_path, nrow(df)))
}

cat("\nMORPH3 probability prediction completed for all years.\n")