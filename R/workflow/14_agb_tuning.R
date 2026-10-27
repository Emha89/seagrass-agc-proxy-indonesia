# =============================================================================
# 14_agb_tuning.R
# Grid search of the Random Forest settings of the AGB model (Table S3 of the
# paper): five random 70/30 splits for every combination of settings.
#
# For the grid search, the target of each sample is one random value within
# its 95% interval. The final model (15_agb_model.R) is trained on the point
# estimate (AGB_pred) with the settings selected here.
#
# Inputs   result/dataAGB_<year>.csv (built by 13_agb_data.R)
# Outputs  result/grid_rf_tuning_AGB.csv, result/grid_rf_best_AGB.csv
# =============================================================================

# -----------------------------------------------------------------------------
# STEP 0: Libraries and functions
# -----------------------------------------------------------------------------
library(dplyr)
library(readr)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "functions", "rf_regression_functions.R"))   # grid search
source(here("R", "workflow", "13_agb_data.R"))                # prepare_agb_training_data()

result_dir <- here("result")

# -----------------------------------------------------------------------------
# STEP 1: Training rows with a valid point estimate and 95% interval
# -----------------------------------------------------------------------------
training_df <- prepare_agb_training_data()

training_df <- training_df %>%
  filter(!is.na(AGB_pred), AGB_pred > 0,
         !is.na(AGB_low), !is.na(AGB_up),
         AGB_low >= 0, AGB_up > AGB_low)

# -----------------------------------------------------------------------------
# STEP 2: Target for the grid search: one random value within the 95% interval
# -----------------------------------------------------------------------------
set.seed(42)
training_df_final <- training_df %>%
  filter(!is.na(AGB_pred), AGB_pred > 0,
         !is.na(AGB_low), !is.na(AGB_up),
         AGB_low >= 0, AGB_up > AGB_low) %>%
  rowwise() %>%
  mutate(
    AGB_simulated = runif(1, AGB_low, AGB_up)
  ) %>%
  ungroup()

cat(sprintf("AGB simulation complete: %d rows.\n", nrow(training_df_final)))

# -----------------------------------------------------------------------------
# STEP 3: Predictors and response
# -----------------------------------------------------------------------------
env_vars <- c("depth", "distToLand")
gse_vars <- paste0("GSE_A", sprintf("%02d", 0:63))

agb_predictors <- c(env_vars, gse_vars)
response_var <- "AGB_simulated"

training_df_final <- training_df_final %>%
  filter(!if_any(all_of(c(response_var, agb_predictors)), is.na))

cat(sprintf("After filtering NA in predictor and response: %d rows remaining\n", nrow(training_df_final)))

# -----------------------------------------------------------------------------
# STEP 4: Grid search
# -----------------------------------------------------------------------------
grid_rds <- file.path(result_dir, "grid_rf_tuning_AGB.rds")

tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning results...\n")
  readRDS(grid_rds)
} else {
  cat("Running tuning...\n")
  res <- evaluate_rf_grid_ranger_regression(
    data = training_df_final,
    covars = agb_predictors,
    response = response_var,
    ntree_values = c(300, 500, 700, 900),
    mtry_values = 4:10,
    nodesize_values = c(3, 5, 7, 9),
    sample_fraction_values = c(0.6, 0.7, 0.8, 0.9),
    n_iter = 5,
    train_frac = 0.7,
    seed = 42
  )
  saveRDS(res, grid_rds)
  res
}

write_csv(tune_result$grid_results, file.path(result_dir, "grid_rf_tuning_AGB.csv"))
write_csv(tune_result$best_params,  file.path(result_dir, "grid_rf_best_AGB.csv"))

cat("Selected AGB settings (grid_rf_best_AGB.csv):\n")
print(tune_result$best_params)
