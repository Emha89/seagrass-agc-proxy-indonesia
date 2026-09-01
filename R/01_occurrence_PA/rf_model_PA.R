# =============================================================================
# rf_model_PA.R
# Random Forest model training and confidence-interval evaluation for
# seagrass presence-absence (PA / occurrence probability) classification.
#
# Pipeline position: stage 1 of 6. Follows dataPrep_func_PA.R (which
# produces dataPA_<year>.csv). Its outputs are consumed by the later
# apply-model and evaluation scripts in this folder.
#
# Key outputs written to result/:
#   - final_rf_model_PA.rds          the trained model
#   - rf_predictor_structure_PA.rds  predictor classes/factor levels, needed
#                                    by apply_model_regression()-style
#                                    functions when scoring new data later
#   - training_data_for_GEE_PA.csv   training rows + predictions, formatted
#                                    for transfer to the GEE deployment step
#   - grid_rf_tuning_PA.csv / .rds, threshold_eval_PA.csv/.rds,
#     var_importance_rf_PA.csv, rf_pa_mc_results.csv,
#     rf_pa_mc_summary_ci.csv       supporting tuning/evaluation outputs
#
# Data availability: this script's input (dataPA_<year>.csv, produced by
# dataPrep_func_PA.R) is not included in this repository.
# =============================================================================

# -----------------------------------------------------------------------------
# STEP 0: Load Libraries and Utility Functions
# -----------------------------------------------------------------------------
library(dplyr)
library(readr)
library(randomForest)
library(caret)
library(ggplot2)  # used directly in STEP 10 below
library(here)     # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_class.R"))  # prepare_training_data(), evaluate_rf_grid_ranger_classification(), etc.
source(here("R", "00_shared_functions", "rf_func_vis.R"))    # plot_pa_threshold_accuracy(), plot_rf_tuning_heatmap()

# Result directory (used for both intermediate stage outputs and this
# script's own output)
result_dir <- here("result")
dir.create(result_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# STEP 1: Define Global Predictor Variables and Configuration
# -----------------------------------------------------------------------------
s2_vars <- c("B2_p60", "B2_stddev", "B4_p20", "B4_p60", "B4_p80",
             "B4_stddev", "B8_p20", "B8_p40", "B8_p80", "B8_stddev",
             "B3_p60_corr", "B3_p60_ent", "B3_stddev_neigh",
             "rb_median", "rg_median", "ndvi", "ndwi1")

env_vars <- c("depth", "distToLand")
wave_vars <- c("elevation", "mean_wave_period", "sig_wave_height")
temporal_var <- "year" # not used: some years have very few data points, not enough to model on
gse_vars <- paste0("GSE_A", sprintf("%02d", 0:63))

# Final list of predictors (s2_vars and wave_vars are defined above but
# excluded from the model -- kept here for reference/future use)
pa_predictors <- c(env_vars, gse_vars)
response_var <- "PA"

# -----------------------------------------------------------------------------
# STEP 2: Load Training Data (with optional balancing)
# -----------------------------------------------------------------------------
# Toggle: TRUE = balance classes if PA=0 >> PA=1, FALSE = use the data as-is
use_balanced_data <- TRUE

training_df <- prepare_training_data(input_dir = result_dir, balance = use_balanced_data)

# -----------------------------------------------------------------------------
# DATA CLEANING: Remove rows with NA values in predictors or response
# -----------------------------------------------------------------------------
na_before <- nrow(training_df)

training_df <- training_df %>%
  dplyr::filter(if_all(all_of(c(response_var, pa_predictors)), ~ !is.na(.)))

na_after <- nrow(training_df)

cat(sprintf("Removed %d rows with NA in predictors/response. Final rows: %d\n\n",
            na_before - na_after, na_after))

# Additional info
cat("PA distribution (after loading):\n")
print(table(training_df$PA))

cat("Years in training data:\n")
print(table(training_df$year))

cat(sprintf("Final training data ready: %d rows\n\n", nrow(training_df)))

# -----------------------------------------------------------------------------
# STEP 3: Run or Load Hyperparameter Tuning
# -----------------------------------------------------------------------------
grid_rds <- file.path(result_dir, "grid_rf_tuning_PA.rds")

if (file.exists(grid_rds)) {
  cat("Loading cached tuning results...\n")
  tune_result <- readRDS(grid_rds)
} else {
  cat("Running Random Forest tuning...\n")
  tune_result <- evaluate_rf_grid_ranger_classification(
    data = training_df,
    covars = pa_predictors,
    response = response_var,
    ntree_values = c(300, 500, 700),
    mtry_values = 4:10,
    nodesize_values = c(3, 5, 7),
    sample_fraction_values = c(0.6, 0.7, 0.8),
    n_iter = 5,
    train_frac = 0.7,
    seed = 42
  )
  saveRDS(tune_result, grid_rds)
}

# Save tuning results
write_csv(tune_result$grid_results, file.path(result_dir, "grid_rf_tuning_PA.csv"))
write_csv(tune_result$best_params,  file.path(result_dir, "grid_rf_best_PA.csv"))

# Extract best parameters
best_params <- tune_result$best_params
model_params <- list(
  ntree = best_params$ntree[1],
  mtry = best_params$mtry[1],
  nodesize = best_params$nodesize[1],
  sample_fraction = best_params$sample_fraction[1]
)

print(model_params)

# -----------------------------------------------------------------------------
# STEP 4: Evaluate Classification Threshold
# -----------------------------------------------------------------------------
threshold_result <- evaluate_thresholds_classification(
  data = training_df,
  covars = pa_predictors,
  response = response_var,
  thresholds = seq(0.1, 0.9, 0.1),
  n_iter = 10,
  ntree = model_params$ntree,
  mtry = model_params$mtry
)

saveRDS(threshold_result, file.path(result_dir, "threshold_eval_PA.rds"))
write_csv(threshold_result, file.path(result_dir, "threshold_eval_PA.csv"))

# The automatically computed best-accuracy threshold is reported below for
# comparison, but the operating threshold actually used is set manually to
# 0.6 (a deliberate choice, not derived from auto_thresh).
auto_thresh <- threshold_result$threshold[which.max(threshold_result$accuracy)]
model_params$threshold <- 0.6

cat(sprintf("Auto best threshold: %.2f | Final threshold used: %.2f\n",
            auto_thresh, model_params$threshold))

# -----------------------------------------------------------------------------
# STEP 5: Train Final RF Model
# -----------------------------------------------------------------------------
rf_formula <- build_rf_formula(pa_predictors, response_var)

final_rf <- randomForest(
  formula = rf_formula,
  data = training_df[, c(response_var, pa_predictors)],
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df)),
  replace = FALSE,
  importance = TRUE
)

saveRDS(final_rf, file.path(result_dir, "final_rf_model_PA.rds"))

# -----------------------------------------------------------------------------
# STEP 6: Predict on Training Data and Export
# -----------------------------------------------------------------------------
training_df$PA_prob <- predict(final_rf, newdata = training_df, type = "prob")[, "1"]
training_df$PA_pred <- ifelse(training_df$PA_prob > model_params$threshold, 1, 0)

extra_cols <- c("lat", "lon", "location", "gee_id", "composition_id", "year")
cols_to_save <- c("PA", pa_predictors, "PA_prob", "PA_pred", intersect(extra_cols, colnames(training_df)))

write_csv(training_df[, cols_to_save], file.path(result_dir, "training_data_for_GEE_PA.csv"))
cat("Saved predictions to 'training_data_for_GEE_PA.csv'\n")

# -----------------------------------------------------------------------------
# STEP 7: Final Model Performance (Confusion Matrix)
# -----------------------------------------------------------------------------
rf_pred <- training_df$PA_pred
actuals <- as.numeric(as.character(training_df$PA))

conf_mat <- confusionMatrix(
  factor(rf_pred, levels = c(0, 1)),
  factor(actuals, levels = c(0, 1)),
  positive = "1"
)

cat("\n==============================\n")
cat("Final PA Model Performance (Training Data)\n")
cat("==============================\n")
print(conf_mat)

# -----------------------------------------------------------------------------
# STEP 8: Variable Importance and Predictor Structure
# -----------------------------------------------------------------------------
var_imp_df <- as.data.frame(importance(final_rf))
var_imp_df$variable <- rownames(var_imp_df)

var_imp_df <- var_imp_df %>%
  select(variable, importance = MeanDecreaseGini) %>%
  arrange(desc(importance))

write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_PA.csv"))

# Predictor classes and factor levels recorded here so that new data can be
# coerced to match exactly at prediction time (see apply_model_regression()
# in rf_func_reg.R for why this matters).
predictor_structure <- list(
  classes = sapply(training_df[, pa_predictors], class),
  levels = lapply(training_df[, pa_predictors], function(x) if (is.factor(x)) levels(x) else NULL)
)
saveRDS(predictor_structure, file.path(result_dir, "rf_predictor_structure_PA.rds"))

cat("Saved variable importance and predictor metadata\n")

# -----------------------------------------------------------------------------
# STEP 9: Monte Carlo Evaluation for CI Estimation
# -----------------------------------------------------------------------------
cat("\nRunning Monte Carlo Evaluation (n_iter = 100)...\n")
mc_result <- run_montecarlo_rf_classification(
  data = training_df,
  covars = pa_predictors,
  response = response_var,
  n_iter = 100,
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  threshold = model_params$threshold,
  train_frac = 0.7,
  seed = 42
)

write_csv(mc_result, file.path(result_dir, "rf_pa_mc_results.csv"))

# CI Summary
summary_ci <- summarize_ci_classification(mc_result)
write_csv(summary_ci, file.path(result_dir, "rf_pa_mc_summary_ci.csv"))

cat("\nMonte Carlo CI summary saved to rf_pa_mc_summary_ci.csv\n")
print(summary_ci)

# -----------------------------------------------------------------------------
# STEP 10: Visualization
# -----------------------------------------------------------------------------
if (file.exists(file.path(result_dir, "threshold_eval_PA.csv"))) {
  print(plot_pa_threshold_accuracy(result_dir, threshold = model_params$threshold))
}

# Heatmap of tuning grid results
grid_file <- file.path(result_dir, "grid_rf_tuning_PA.csv")
if (file.exists(grid_file)) {
  grid_df <- read_csv(grid_file, show_col_types = FALSE)
  cat("\nPlotting tuning heatmap (accuracy)...\n")
  print(plot_rf_tuning_heatmap(grid_df, fill_var = "accuracy_mean",
                               title = "Grid Tuning Heatmap (Accuracy)"))
}

# Plot variable importance
varimp_file <- file.path(result_dir, "var_importance_rf_PA.csv")
if (file.exists(varimp_file)) {
  var_imp_df <- read_csv(varimp_file, show_col_types = FALSE) %>%
    arrange(desc(importance)) %>%
    slice_head(n = 20)  # top 20 only
  
  cat("\nPlotting top 20 variable importance...\n")
  p_varimp <- ggplot(var_imp_df, aes(x = reorder(variable, importance), y = importance)) +
    geom_col(fill = "grey") +
    coord_flip() +
    labs(title = "Top 20 Variable Importance (MeanDecreaseGini)",
         x = "Predictor", y = "Importance") +
    theme_minimal(base_size = 13)
  
  print(p_varimp)
}