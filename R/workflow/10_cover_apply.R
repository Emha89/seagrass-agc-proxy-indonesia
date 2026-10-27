# =============================================================================
# 10_cover_apply.R
# Percent cover (tSPC_pred) at the field points of every year.
#
# Inputs   result/predicted_PA_<year>.csv (second run of step 03),
#          result/final_rf_model_PCT.rds, result/rf_predictor_structure_PCT.rds
# Outputs  result/predicted_tSPC_<year>.csv
# =============================================================================

library(dplyr)
library(readr)
library(randomForest)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "functions", "rf_regression_functions.R"))  # includes apply_model_regression()

cat("Loading model and predictor structure...\n")

result_dir <- here("result")

model <- readRDS(file.path(result_dir, "final_rf_model_PCT.rds"))
predictor_structure <- readRDS(file.path(result_dir, "rf_predictor_structure_PCT.rds"))
pct_predictors <- names(predictor_structure$classes)

data_years <- as.character(2017:2024)

for (yr in data_years) {
  cat(sprintf("\nRunning prediction for year: %s\n", yr))
  
  in_path <- file.path(result_dir, paste0("predicted_PA_", yr, ".csv"))
  out_path <- file.path(result_dir, paste0("predicted_tSPC_", yr, ".csv"))
  
  if (!file.exists(in_path)) {
    warning(sprintf("File not found: %s. Skipping year %s.\n", in_path, yr))
    next
  }
  
  df <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)
  
  # Ensure predictors exist and valid
  missing <- setdiff(pct_predictors, colnames(df))
  if (length(missing) > 0) {
    warning(sprintf("Missing predictors in %s: %s", yr, paste(missing, collapse = ", ")))
    next
  }
  
  valid_idx <- df %>%
    select(all_of(pct_predictors)) %>%
    complete.cases()
  
  df_valid <- df[valid_idx, ]
  
  # Predict SPC
  df$tSPC_pred <- NA_real_
  df$tSPC_pred[valid_idx] <- apply_model_regression(
    df = df_valid,
    predictors = pct_predictors,
    model = model,
    predictor_structure = predictor_structure
  )
  
  # Clamp predictions to 0-100
  n_below_0 <- sum(df$tSPC_pred < 0, na.rm = TRUE)
  n_above_100 <- sum(df$tSPC_pred > 100, na.rm = TRUE)
  if (n_below_0 + n_above_100 > 0) {
    cat(sprintf("Clamped %d predictions below 0 and %d above 100 to [0-100].\n",
                n_below_0, n_above_100))
  }
  df$tSPC_pred <- pmin(pmax(df$tSPC_pred, 0), 100)
  
  # Save output
  write_csv(df, out_path)
  cat(sprintf("Saved: %s (%d rows)\n", out_path, nrow(df)))
}