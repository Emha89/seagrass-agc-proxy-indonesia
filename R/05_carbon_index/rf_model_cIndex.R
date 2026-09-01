# =============================================================================
# rf_model_combined_CINDEX.R
# Random Forest regression for carbon index (0-1), predicting it from
# satellite/environmental predictors so it can be estimated across the
# full study extent, not just at points with field species-composition
# data. Combined script: data prep, tuning, training, prediction (on
# training data), performance summary, and Monte Carlo CI estimation all
# in one file (same "combined" pattern as morphology and SPC).
#
# Pipeline position: stage 5 of 6. Reads predicted_PA_<year>.csv (from
# 01_occurrence_PA) and dataCARBON_<year>.csv (from
# carbon_scaling_join_SP_finalDell.R in this folder, for the carbon_index
# response).
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

source(here("R", "00_shared_functions", "rf_func_reg.R"))
source(here("R", "00_shared_functions", "rf_func_vis.R"))

result_dir <- here("result")

# -----------------------------------------------------------------------------
# STEP 1: Prepare Training Data (Filtered by GT Year)
# -----------------------------------------------------------------------------
prepare_training_data_cindex <- function(
    input_dir = result_dir,
    years = c("2017","2018","2019","2020","2021","2022","2023"),
    response = "carbon_index"
) {
  
  df_all <- purrr::map_dfr(years, function(yr) {
    
    path <- file.path(input_dir, paste0("predicted_PA_", yr, ".csv"))
    carbon_path <- file.path(input_dir, paste0("dataCARBON_", yr, ".csv"))
    
    # -------------------------------------
    # CHECK predicted_PA file
    # -------------------------------------
    if (!file.exists(path)) {
      message(sprintf("File not found: %s", path))
      return(NULL)
    }
    
    df_tmp <- readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
    
    df_tmp <- df_tmp %>%
      dplyr::filter(year == yr)
    
    # -------------------------------------
    # JOIN carbon_index
    # -------------------------------------
    if (!file.exists(carbon_path)) {
      message(sprintf("Carbon file not found: %s", carbon_path))
      return(NULL)   # Important: do not continue if the carbon file is missing
    }
    
    carbon_df <- readr::read_csv(carbon_path, show_col_types = FALSE) %>%
      dplyr::select(gee_id, carbon_index)
    
    df_tmp <- df_tmp %>% left_join(carbon_df, by = "gee_id")
    
    # -------------------------------------
    # FILTER carbon_index range 0-1
    # -------------------------------------
    n_before <- nrow(df_tmp)
    
    df_tmp <- df_tmp %>%
      dplyr::filter(
        !is.na(carbon_index),
        carbon_index > 0,
        carbon_index <= 1
      ) %>%
      dplyr::mutate(year = as.factor(yr))
    
    n_after <- nrow(df_tmp)
    
    cat(sprintf(
      "Removed %d rows with NA or invalid carbon_index for year %s\n",
      n_before - n_after, yr
    ))
    
    return(df_tmp)
  })
  
  # -------------------------------------
  # Summary stats
  # -------------------------------------
  if (nrow(df_all) == 0) {
    stop("ERROR: No training samples available after filtering. Please check data availability.")
  }
  
  year_counts <- df_all %>%
    dplyr::count(year, name = "n_samples") %>%
    dplyr::arrange(year)
  
  cat("Number of GT samples per year (after carbon_index filter):\n")
  print(year_counts)
  
  message(sprintf("Total training samples retained: %d rows.", nrow(df_all)))
  
  return(df_all)
}

training_df <- prepare_training_data_cindex()

# -----------------------------------------------------------------------------
# STEP 2: Define Predictors and Response (carbon_index)
# -----------------------------------------------------------------------------
s2_vars <- c("B2_p60", "B2_stddev", "B4_p20", "B4_p60", "B4_p80",
             "B4_stddev", "B8_p20", "B8_p40", "B8_p80", "B8_stddev",
             "B3_p60_corr", "B3_p60_ent", "B3_stddev_neigh",
             "rb_median", "rg_median", "ndvi", "ndwi1")

env_vars <- c("depth", "distToLand")   

wave_vars <- c("elevation", "mean_wave_period", "sig_wave_height")
temporal <- "year"
gse_vars <- paste0("GSE_A", sprintf("%02d", 0:63))
pred <- "PA_prob"  # defined, deliberately excluded below (see cindex_predictors)

# wave_vars and pred (PA_prob) are deliberately excluded, to avoid
# predictor redundancy across proxy stages (same principle applied in
# morphology, SPC, and AGB)
cindex_predictors <- c(env_vars, gse_vars)
response_var <- "carbon_index"

# -----------------------------------------------------------------------------
# STEP 3: Grid Tuning using ranger
# -----------------------------------------------------------------------------
grid_rds <- file.path(result_dir, "grid_rf_tuning_CINDEX.rds")

tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning results...\n")
  readRDS(grid_rds)
} else {
  cat("Running tuning...\n")
  res <- evaluate_rf_grid_ranger_regression(
    data = training_df,
    covars = cindex_predictors,
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

write_csv(tune_result$grid_results, file.path(result_dir, "grid_rf_tuning_CINDEX.csv"))
write_csv(tune_result$best_params,  file.path(result_dir, "grid_rf_best_CINDEX.csv"))

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
rf_formula <- build_rf_formula(cindex_predictors, response_var)

final_rf <- randomForest(
  formula = rf_formula,
  data = training_df[, c(response_var, cindex_predictors)],
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df)),
  replace = FALSE,
  importance = TRUE
)

saveRDS(final_rf, file.path(result_dir, "final_rf_model_CINDEX.rds"))

# -----------------------------------------------------------------------------
# STEP 5: Predict on Training Set & Export
# -----------------------------------------------------------------------------
training_df$carbon_index_pred <- predict(final_rf, newdata = training_df)

extra_cols <- c("lat", "lon", "location", "gee_id", "composition_id", "year")
cols_to_save <- c("carbon_index", cindex_predictors,
                  "carbon_index_pred",
                  intersect(extra_cols, colnames(training_df)))

write_csv(training_df[, cols_to_save],
          file.path(result_dir, "training_data_for_GEE_CINDEX.csv"))

# -----------------------------------------------------------------------------
# STEP 6: Final Model Performance
# -----------------------------------------------------------------------------
training_df_eval <- training_df %>%
  filter(!is.na(carbon_index_pred),
         carbon_index_pred >= 0,
         carbon_index_pred <= 1)

cat("\n==============================\n")
cat("Final Carbon Index Model Performance (Training Data)\n")
cat("==============================\n")
cat(sprintf("RMSE: %.4f | MAE: %.4f | R2: %.4f\n",
            RMSE(training_df_eval$carbon_index_pred, training_df_eval$carbon_index),
            MAE(training_df_eval$carbon_index_pred, training_df_eval$carbon_index),
            R2(training_df_eval$carbon_index_pred, training_df_eval$carbon_index)))

# -----------------------------------------------------------------------------
# STEP 7: Variable Importance
# -----------------------------------------------------------------------------
var_imp_df <- as.data.frame(final_rf$importance)
var_imp_df$variable <- rownames(var_imp_df)

var_imp_df <- var_imp_df %>%
  select(variable, importance = IncNodePurity) %>%
  arrange(desc(importance))

write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_CINDEX.csv"))

# Predictor Structure (needed for apply model)
predictor_structure <- list(
  classes = sapply(training_df[, cindex_predictors], class),
  levels = lapply(training_df[, cindex_predictors],
                  function(x) if (is.factor(x)) levels(x) else NULL)
)
saveRDS(predictor_structure,
        file.path(result_dir, "rf_predictor_structure_CINDEX.rds"))

# -----------------------------------------------------------------------------
# STEP 8: Visualization
# -----------------------------------------------------------------------------
if (file.exists(file.path(result_dir, "var_importance_rf_CINDEX.csv"))) {
  
  var_imp_top <- read_csv(file.path(result_dir, "var_importance_rf_CINDEX.csv"),
                          show_col_types = FALSE) %>%
    arrange(desc(importance)) %>%
    slice_head(n = 20)
  
  p_varimp <- ggplot(var_imp_top,
                     aes(x = reorder(variable, importance), y = importance)) +
    geom_col(fill = "grey") +
    coord_flip() +
    labs(title = "Top 20 Variable Importance (Carbon Index)",
         x = "Variable", y = "IncNodePurity") +
    theme_minimal(base_size = 13)
  
  print(p_varimp)
}

# -----------------------------------------------------------------------------
# STEP 9: Monte Carlo Evaluation
# -----------------------------------------------------------------------------
cat("\nRunning Monte Carlo Evaluation (n_iter = 100)...\n")

training_df_mc <- training_df %>%
  filter(!is.na(carbon_index),
         !is.na(carbon_index_pred),
         carbon_index > 0,
         carbon_index <= 1)

mc_result <- run_montecarlo_rf_regression(
  data = training_df_mc,
  covars = cindex_predictors,
  response = response_var,
  n_iter = 100,
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  nodesize = model_params$nodesize,
  sample_fraction = model_params$sample_fraction,
  train_frac = 0.7,
  seed = 42
)

write_csv(mc_result, file.path(result_dir, "rf_cindex_mc_results.csv"))

# -----------------------------------------------------------------------------
# STEP 10: CI Summary
# -----------------------------------------------------------------------------
summary_ci <- summarize_ci_regression(mc_result)
write_csv(summary_ci, file.path(result_dir, "rf_cindex_mc_summary_ci.csv"))

cat("\nMonte Carlo CI summary saved to rf_cindex_mc_summary_ci.csv\n")
print(summary_ci)