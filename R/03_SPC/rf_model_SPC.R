# =============================================================================
# rf_model_combined_SPC.R
# Random Forest regression for seagrass percent cover (SPC / tSPC).
# Combined script: data prep, tuning, training, prediction (on training
# data), performance summary, and Monte Carlo CI estimation all in one
# file (same "combined" pattern as rf_model_combined_MORPH3.R).
#
# CONFIRMED VERSION: SPC predicts from morphology class probabilities
# (P_mixed_short_plus_mono_short, P_mixed_long, P_mono_Ea), not from
# carbon_index or PA_prob -- both are explicitly excluded. An earlier
# draft of this script used carbon_index as a predictor instead; that
# version has been replaced by this one.
#
# Pipeline position: stage 3 of 6. Trains on dataPA_<year>.csv (from
# 01_occurrence_PA), which needs to already include the morphology P_*
# columns -- i.e. dataPrep_func_PA.R must have been re-run after
# 02_morphology/rf_applyModel_MORPH3.R has produced
# predicted_MORPH3_probs_<year>.csv. See the pipeline ordering note in
# dataPrep_func_PA.R for this same multi-pass requirement.
#
# Key output: tSPC_pred, used as a predictor in later stages (AGB, carbon
# index regression, AGC).
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

# -----------------------------------------------------------------------------
# STEP 0: Load Libraries and Utility Functions
# -----------------------------------------------------------------------------
library(dplyr)
library(readr)
library(randomForest)
library(caret)
library(ranger)
library(purrr)
library(ggplot2)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_reg.R"))   # formula builder, tuning, Monte Carlo, CI summarizer
source(here("R", "00_shared_functions", "rf_func_vis.R"))   # visualization functions

result_dir <- here("result")

# -----------------------------------------------------------------------------
# STEP 1: Prepare Training Data from Annual dataPA Files (GT + Predictors Ready)
# -----------------------------------------------------------------------------
prepare_training_data <- function(
    input_dir = result_dir,
    years = c("2017", "2018", "2019", "2020", "2021", "2022", "2023"),
    response = "tSPC"
) {
  
  df_all <- purrr::map_dfr(years, function(yr) {
    
    # -------------------------------------------------------------------------
    # Load annual prepared dataset: dataPA_YYYY.csv
    # This file already contains Morph3 probability predictors (P_*)
    # -------------------------------------------------------------------------
    path <- file.path(input_dir, paste0("dataPA_", yr, ".csv"))
    
    if (!file.exists(path)) {
      message(sprintf("File not found: %s", path))
      return(NULL)
    }
    
    df_tmp <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
    
    # -------------------------------------------------------------------------
    # Filter valid GT samples for SPC regression
    # -------------------------------------------------------------------------
    df_tmp <- df_tmp %>%
      dplyr::filter(!is.na(.data[[response]])) %>%   # response must exist
      dplyr::filter(.data[[response]] >= 0,
                    .data[[response]] <= 100) %>%
      dplyr::mutate(year = as.factor(yr))
    
    return(df_tmp)
  })
  
  # ---------------------------------------------------------------------------
  # Report sample count per year
  # ---------------------------------------------------------------------------
  year_counts <- df_all %>%
    dplyr::count(year, name = "n_samples") %>%
    dplyr::arrange(year)
  
  cat("GT samples retained per year for SPC model:\n")
  print(year_counts)
  
  message(sprintf("Total training samples retained: %d rows.", nrow(df_all)))
  
  return(df_all)
}

training_df <- prepare_training_data(
  input_dir = result_dir,
  years = c("2017", "2018", "2019", "2020", "2021", "2022", "2023"),
  response = "tSPC"
)

# -----------------------------------------------------------------------------
# STEP 2: Define Predictors and Response Variable
# -----------------------------------------------------------------------------
s2_vars <- c("B2_p60", "B2_stddev", "B4_p20", "B4_p60", "B4_p80",
             "B4_stddev", "B8_p20", "B8_p40", "B8_p80", "B8_stddev",
             "B3_p60_corr", "B3_p60_ent", "B3_stddev_neigh",
             "rb_median", "rg_median", "ndvi", "ndwi1")

env_vars <- c("depth", "distToLand")

wave_vars <- c("elevation", "mean_wave_period", "sig_wave_height")
temporal <- "year"
gse_vars <- paste0("GSE_A", sprintf("%02d", 0:63))
pred <- "PA_prob"  # defined, deliberately excluded below (see pct_predictors)

morph_pred <- c("P_mixed_short_plus_mono_short", "P_mixed_long", "P_mono_Ea")

# wave_vars and pred (PA_prob) are deliberately excluded, to avoid predictor
# redundancy across proxy stages (same principle applied in morphology)
pct_predictors <- c(env_vars, gse_vars, morph_pred)
response_var <- "tSPC"

# -----------------------------------------------------------------------------
# STEP 3: Grid Tuning using ranger
# -----------------------------------------------------------------------------
grid_rds <- file.path(result_dir, "grid_rf_tuning_PCT.rds")

tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning results...\n")
  readRDS(grid_rds)
} else {
  cat("Running tuning...\n")
  res <- evaluate_rf_grid_ranger_regression(
    data = training_df,
    covars = pct_predictors,
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

write_csv(tune_result$grid_results, file.path(result_dir, "grid_rf_tuning_PCT.csv"))
write_csv(tune_result$best_params,  file.path(result_dir, "grid_rf_best_PCT.csv"))

best_params <- tune_result$best_params
model_params <- list(
  ntree = best_params$ntree[1],
  mtry = best_params$mtry[1],
  nodesize = best_params$nodesize[1],
  sample_fraction = best_params$sample_fraction[1]
)

# -----------------------------------------------------------------------------
# STEP 4: Train Final Model using randomForest
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Remove rows with NA in response or predictors (required by randomForest)
# -----------------------------------------------------------------------------
missing_cols <- c(response_var, pct_predictors)

na_before <- sum(!complete.cases(training_df[, missing_cols]))
cat("Rows with NA before training:", na_before, "\n")

training_df <- training_df %>%
  filter(complete.cases(across(all_of(missing_cols))))

cat("Rows remaining after NA removal:", nrow(training_df), "\n")

# -----------------------------------------------------------------------------

rf_formula <- build_rf_formula(pct_predictors, response_var)

final_rf <- randomForest(
  formula = rf_formula,
  data = training_df[, c(response_var, pct_predictors)],
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df)),
  replace = FALSE,
  importance = TRUE
)

saveRDS(final_rf, file.path(result_dir, "final_rf_model_PCT.rds"))

# -----------------------------------------------------------------------------
# STEP 5: Predict on Training Set & Export for GEE
# -----------------------------------------------------------------------------
training_df$tSPC_pred <- predict(final_rf, newdata = training_df)

extra_cols <- c("lat", "lon", "location", "gee_id", "composition_id", "year")
cols_to_save <- c("tSPC", pct_predictors, "tSPC_pred", intersect(extra_cols, colnames(training_df)))

write_csv(training_df[, cols_to_save], file.path(result_dir, "training_data_for_GEE_PCT.csv"))

# -----------------------------------------------------------------------------
# STEP 6: Final Model Performance
# -----------------------------------------------------------------------------
training_df_eval <- training_df %>%
  filter(!is.na(tSPC_pred), tSPC_pred >= 0, tSPC_pred <= 100)

cat("\n==============================\n")
cat("Final SPC Model Performance (Training Data)\n")
cat("==============================\n")
cat(sprintf("RMSE: %.3f | MAE: %.3f | R2: %.3f\n",
            RMSE(training_df_eval$tSPC_pred, training_df_eval$tSPC),
            MAE(training_df_eval$tSPC_pred, training_df_eval$tSPC),
            R2(training_df_eval$tSPC_pred, training_df_eval$tSPC)))

# -----------------------------------------------------------------------------
# STEP 7: Variable Importance
# -----------------------------------------------------------------------------
var_imp_df <- as.data.frame(final_rf$importance)
var_imp_df$variable <- rownames(var_imp_df)

var_imp_df <- var_imp_df %>%
  select(variable, importance = IncNodePurity) %>%
  arrange(desc(importance))

write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_PCT.csv"))

# Save predictor metadata (for apply model step)
predictor_structure <- list(
  classes = sapply(training_df[, pct_predictors], class),
  levels = lapply(training_df[, pct_predictors], function(x) if (is.factor(x)) levels(x) else NULL)
)
saveRDS(predictor_structure, file.path(result_dir, "rf_predictor_structure_PCT.rds"))

# -----------------------------------------------------------------------------
# STEP 8: Visualization
# -----------------------------------------------------------------------------
if (file.exists(file.path(result_dir, "var_importance_rf_PCT.csv"))) {
  var_imp_top <- read_csv(file.path(result_dir, "var_importance_rf_PCT.csv"), show_col_types = FALSE) %>%
    arrange(desc(importance)) %>% slice_head(n = 20)
  p_varimp <- ggplot(var_imp_top, aes(x = reorder(variable, importance), y = importance)) +
    geom_col(fill = "grey") +
    coord_flip() +
    labs(title = "Top 20 Variable Importance (SPC)", x = "Variable", y = "IncNodePurity") +
    theme_minimal(base_size = 13)
  print(p_varimp)
}

# -----------------------------------------------------------------------------
# STEP 9: Monte Carlo Evaluation
# -----------------------------------------------------------------------------
cat("\nRunning Monte Carlo Evaluation (n_iter = 100)...\n")
training_df_mc <- training_df %>%
  filter(!is.na(tSPC), !is.na(tSPC_pred), tSPC > 0, tSPC <= 100)

mc_result <- run_montecarlo_rf_regression(
  data = training_df_mc,
  covars = pct_predictors,
  response = response_var,
  n_iter = 100,
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  nodesize = model_params$nodesize,
  sample_fraction = model_params$sample_fraction,
  train_frac = 0.7,
  seed = 42
)

write_csv(mc_result, file.path(result_dir, "rf_pct_mc_results.csv"))

# -----------------------------------------------------------------------------
# STEP 10: CI Summary
# -----------------------------------------------------------------------------
summary_ci <- summarize_ci_regression(mc_result)
write_csv(summary_ci, file.path(result_dir, "rf_pct_mc_summary_ci.csv"))

cat("\nMonte Carlo CI summary saved to rf_pct_mc_summary_ci.csv\n")
print(summary_ci)