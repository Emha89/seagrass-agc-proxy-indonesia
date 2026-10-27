# =============================================================================
# 17_agc_tuning.R
# Grid search of the Random Forest settings of the AGC model (Table S3 of the
# paper): five random 70/30 splits for every combination of settings.
#
# For the grid search, the target of each sample is one random value within
# its 95% interval. The final model (18_agc_model.R) is trained on the point
# estimate (AGC_pred) with the selected settings (900 trees, mtry 5, node
# size 9, sampling fraction 0.6).
#
# Inputs   result/dataAGC_<year>.csv (built by 16_agc_data.R)
# Outputs  result/grid_rf_tuning_AGC.csv, result/grid_rf_best_AGC.csv
# =============================================================================

# -----------------------------------------------------------------------------
# STEP 0: Libraries and functions
# -----------------------------------------------------------------------------
library(dplyr)
library(readr)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "functions", "rf_regression_functions.R"))   # grid search
source(here("R", "workflow", "16_agc_data.R"))                # prepare_agc_training_data()

result_dir <- here("result")

# -----------------------------------------------------------------------------
# STEP 1: Training rows with a valid point estimate and 95% interval
# -----------------------------------------------------------------------------
training_df <- prepare_agc_training_data()

training_df <- training_df %>%
  filter(!is.na(AGC_pred), AGC_pred > 0,
         !is.na(AGC_low), !is.na(AGC_up),
         AGC_low >= 0, AGC_up > AGC_low)

# -----------------------------------------------------------------------------
# STEP 2: Target for the grid search: one random value within the 95% interval
# -----------------------------------------------------------------------------
set.seed(42)
training_df_final <- training_df %>%
  filter(!is.na(AGC_pred), AGC_pred > 0,
         !is.na(AGC_low), !is.na(AGC_up),
         AGC_low >= 0, AGC_up > AGC_low) %>%
  rowwise() %>%
  mutate(AGC_simulated = runif(1, AGC_low, AGC_up)) %>%
  ungroup()

cat(sprintf("AGC simulation complete: %d rows.\n", nrow(training_df_final)))

# -----------------------------------------------------------------------------
# STEP 3: Predictors and response
# -----------------------------------------------------------------------------
env_vars   <- c("depth", "distToLand")
gse_vars   <- paste0("GSE_A", sprintf("%02d", 0:63))
proxy_vars <- c("PA_prob", "tSPC_pred", "tAGB_pred", "carbon_index")
morph_pred <- c("P_mixed_short_plus_mono_short", "P_mixed_long", "P_mono_Ea")

agc_predictors <- c(gse_vars, env_vars, proxy_vars, morph_pred)
response_var <- "AGC_simulated"

missing_vars <- setdiff(c(response_var, agc_predictors), names(training_df_final))
if (length(missing_vars) > 0) {
  stop(sprintf("Missing required variables in training data: %s",
               paste(missing_vars, collapse = ", ")))
}

training_df_final <- training_df_final %>%
  filter(!if_any(all_of(c(response_var, agc_predictors)), is.na))

cat(sprintf("After filtering NA in predictor and response: %d rows remaining\n", nrow(training_df_final)))

# -----------------------------------------------------------------------------
# STEP 4: Grid search
#   The mtry range is narrower than for the proxy models (4:10), because the
#   AGC model has more predictors that are strongly correlated.
# -----------------------------------------------------------------------------
grid_rds <- file.path(result_dir, "grid_rf_tuning_AGC.rds")

tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning results...\n")
  readRDS(grid_rds)
} else {
  cat("Running tuning...\n")
  res <- evaluate_rf_grid_ranger_regression(
    data = training_df_final,
    covars = agc_predictors,
    response = response_var,
    ntree_values = c(300, 500, 700, 900),
    mtry_values = 2:5,
    nodesize_values = c(3, 5, 7, 9),
    sample_fraction_values = c(0.6, 0.7, 0.8, 0.9),
    n_iter = 5,
    train_frac = 0.7,
    seed = 42
  )
  saveRDS(res, grid_rds)
  res
}

write_csv(tune_result$grid_results, file.path(result_dir, "grid_rf_tuning_AGC.csv"))
write_csv(tune_result$best_params,  file.path(result_dir, "grid_rf_best_AGC.csv"))

cat("Selected AGC settings (grid_rf_best_AGC.csv):\n")
print(tune_result$best_params)
