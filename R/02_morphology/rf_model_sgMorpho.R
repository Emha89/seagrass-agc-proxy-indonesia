# =============================================================================
# rf_model_combined_MORPH3.R
# Random Forest multi-class classification for leaf morphology (morph3).
# Combined script: data prep, tuning, training, prediction (on training
# data), evaluation, and Monte Carlo CI estimation all in one file (unlike
# the PA stage, which is split across 4 separate scripts).
#
# Pipeline position: stage 2 of 6. Reads predicted_PA_<year>.csv (from
# 01_occurrence_PA) and dataCARBON_<year>.csv (from 05_carbon_index, for
# the morph3 label) -- another instance of a later-stage output being
# needed as input here; see the pipeline ordering note in
# dataPrep_func_PA.R for the same kind of dependency.
#
# Key output: P_* class probability columns
# (training_data_for_GEE_MORPH3_probs.csv), used as predictors in the AGC
# regression stage, and training_MORPH3_points_ALL.csv for GEE upload.
#
# OPEN QUESTIONS (flagged during review, not changed here):
#   - filter(year == yr) is applied to predicted_PA_<year>.csv; this
#     assumes that file has its own `year` column. rf_applyModel_PA.R
#     doesn't appear to add one explicitly, so this only works if the raw
#     GSE_training_<year>_CSV.csv already includes it.
#   - `pred <- "PA_prob"` is defined but not included in morph_predictors
#     (only env_vars + gse_vars are used) -- confirm whether PA_prob is
#     deliberately excluded as a predictor here.
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
library(tibble)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_class.R"))  # classification utilities

set.seed(42)

result_dir <- here("result")

# -----------------------------------------------------------------------------
# STEP 1: Prepare Training Data (Filtered by GT Year)
# -----------------------------------------------------------------------------
prepare_training_data_morph3 <- function(
    input_dir = result_dir,
    years = c("2017","2018","2019","2020","2021","2022","2023"),
    response = "morph3"
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
      filter(year == yr)
    
    # -------------------------------------
    # JOIN morph3 label from Carbon dataset
    # -------------------------------------
    if (!file.exists(carbon_path)) {
      message(sprintf("Carbon file not found: %s", carbon_path))
      return(NULL)
    }
    
    carbon_df <- readr::read_csv(carbon_path, show_col_types = FALSE) %>%
      select(gee_id, morph3)
    
    df_tmp <- df_tmp %>% left_join(carbon_df, by = "gee_id")
    
    # -------------------------------------
    # FILTER valid morph3 classes
    # -------------------------------------
    n_before <- nrow(df_tmp)
    
    df_tmp <- df_tmp %>%
      filter(!is.na(morph3)) %>%
      mutate(
        year = as.factor(yr),
        morph3 = as.factor(morph3)
      )
    
    n_after <- nrow(df_tmp)
    
    cat(sprintf(
      "Removed %d rows with NA morph3 for year %s\n",
      n_before - n_after, yr
    ))
    
    return(df_tmp)
  })
  
  # -------------------------------------
  # Summary stats
  # -------------------------------------
  if (nrow(df_all) == 0) {
    stop("ERROR: No morphology samples available after filtering.")
  }
  
  year_counts <- df_all %>%
    count(year, name = "n_samples") %>%
    arrange(year)
  
  cat("\nTraining samples per year (morph3):\n")
  print(year_counts)
  
  cat("\nMorphology class distribution (overall):\n")
  print(df_all %>% count(morph3) %>% mutate(prop = n / sum(n)))
  
  message(sprintf("Total training samples retained: %d rows.", nrow(df_all)))
  
  return(df_all)
}

training_df <- prepare_training_data_morph3()

# =============================================================================
# STEP 1B: Export MORPH3 Training Points (Combined File for GEE Upload)
# Supports both coordinate styles:
#   - lat/lon
#   - xcoord/ycoord
# =============================================================================

cat("\nPreparing MORPH3 training points for GEE upload...\n")

# ---------------------------------------------------------------------------
# Detect coordinate column names automatically
# ---------------------------------------------------------------------------
coord_cols <- NULL

if (all(c("lon", "lat") %in% names(training_df))) {
  coord_cols <- c("lon", "lat")
  cat("Coordinates detected: lon/lat\n")
} else if (all(c("xcoord", "ycoord") %in% names(training_df))) {
  coord_cols <- c("xcoord", "ycoord")
  cat("Coordinates detected: xcoord/ycoord\n")
} else {
  stop("No coordinate columns found (expected lon/lat or xcoord/ycoord).")
}

# ---------------------------------------------------------------------------
# Export only required fields (label + metadata)
# ---------------------------------------------------------------------------
gee_training_morph3 <- training_df %>%
  
  select(
    gee_id,
    year,
    loc,
    morph3,
    any_of(coord_cols),
    any_of(c("area", "location"))
  ) %>%
  
  mutate(
    year   = as.integer(as.character(year)),
    morph3 = as.character(morph3),
    loc    = as.character(loc)
  ) %>%
  
  # Remove rows without label or coordinates
  filter(
    !is.na(morph3),
    !is.na(.data[[coord_cols[1]]]),
    !is.na(.data[[coord_cols[2]]])
  ) %>%
  
  distinct(gee_id, .keep_all = TRUE)

cat(sprintf("Final MORPH3 training dataset ready: %d rows\n",
            nrow(gee_training_morph3)))

# ---------------------------------------------------------------------------
# Save ONE combined CSV file (multi-year)
# ---------------------------------------------------------------------------
out_path <- file.path(result_dir, "training_MORPH3_points_ALL.csv")
write_csv(gee_training_morph3, out_path)

cat(sprintf("\nSaved combined MORPH3 training file:\n   %s\n", out_path))

cat("\nNext step: Upload this CSV to GEE as a FeatureCollection.\n")
# =============================================================================

# -----------------------------------------------------------------------------
# STEP 2: Define Predictors and Response (morph3)
# -----------------------------------------------------------------------------
env_vars <- c("depth", "distToLand")
gse_vars <- paste0("GSE_A", sprintf("%02d", 0:63))
pred     <- "PA_prob"  # defined but not currently used in morph_predictors -- see note at top of file

morph_predictors <- c(env_vars, gse_vars)
response_var <- "morph3"

# -----------------------------------------------------------------------------
# STEP 3: Grid Tuning using ranger (Classification mode)
# -----------------------------------------------------------------------------
grid_rds <- file.path(result_dir, "grid_rf_tuning_MORPH3.rds")

tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning results...\n")
  readRDS(grid_rds)
} else {
  cat("Running morphology classification tuning...\n")
  
  res <- evaluate_rf_grid_ranger_multiclass(
    data = training_df,
    covars = morph_predictors,
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

write_csv(tune_result$grid_results,
          file.path(result_dir, "grid_rf_tuning_MORPH3.csv"))

write_csv(tune_result$best_params,
          file.path(result_dir, "grid_rf_best_MORPH3.csv"))

best_params <- tune_result$best_params

model_params <- list(
  ntree = best_params$ntree[1],
  mtry = best_params$mtry[1],
  nodesize = best_params$nodesize[1],
  sample_fraction = best_params$sample_fraction[1]
)

# -----------------------------------------------------------------------------
# STEP 4: Train Final Model using randomForest (Multi-class)
# -----------------------------------------------------------------------------
rf_formula <- build_rf_formula(morph_predictors, response_var)

final_rf <- randomForest(
  formula = rf_formula,
  data = training_df[, c(response_var, morph_predictors)],
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df)),
  replace = FALSE,
  importance = TRUE
)

saveRDS(final_rf,
        file.path(result_dir, "final_rf_model_MORPH3.rds"))

# -----------------------------------------------------------------------------
# STEP 5: Predict Probabilities & Export (Key Paper Output)
# -----------------------------------------------------------------------------
cat("\nPredicting morphology class probabilities...\n")

prob_matrix <- predict(final_rf,
                       newdata = training_df,
                       type = "prob")

prob_df <- as.data.frame(prob_matrix)
colnames(prob_df) <- paste0("P_", colnames(prob_df))

training_df_probs <- bind_cols(training_df, prob_df)

# Save for AGC regression predictor input
write_csv(training_df_probs,
          file.path(result_dir, "training_data_for_GEE_MORPH3_probs.csv"))

# -----------------------------------------------------------------------------
# STEP 6: Final Model Performance (Classification Accuracy)
# -----------------------------------------------------------------------------
pred_class <- predict(final_rf, newdata = training_df)

cm <- caret::confusionMatrix(pred_class, training_df$morph3)

cat("\n==============================\n")
cat("Final Morphology Model Performance (Training Data)\n")
cat("==============================\n")

cat(sprintf("Accuracy: %.4f\n", cm$overall["Accuracy"]))
cat(sprintf("Kappa   : %.4f\n", cm$overall["Kappa"]))

cat("\nClass-wise accuracy:\n")
print(round(prop.table(cm$table, margin = 2), 2))

# Save confusion matrix
write_csv(as.data.frame(cm$table),
          file.path(result_dir, "confusion_matrix_MORPH3.csv"))

# -----------------------------------------------------------------------------
# STEP 7: Variable Importance
# -----------------------------------------------------------------------------
var_imp_df <- as.data.frame(final_rf$importance)
var_imp_df$variable <- rownames(var_imp_df)

var_imp_df <- var_imp_df %>%
  arrange(desc(MeanDecreaseGini))

write_csv(var_imp_df,
          file.path(result_dir, "var_importance_rf_MORPH3.csv"))

# -----------------------------------------------------------------------------
# STEP 7b: Save Predictor Structure (needed for apply model)
# -----------------------------------------------------------------------------
predictor_structure <- list(
  classes = sapply(training_df[, morph_predictors], class),
  
  # Store factor levels (important for consistent prediction later)
  levels = lapply(training_df[, morph_predictors],
                  function(x) if (is.factor(x)) levels(x) else NULL)
)

saveRDS(
  predictor_structure,
  file.path(result_dir, "rf_predictor_structure_MORPH3.rds")
)

cat("Predictor structure saved: rf_predictor_structure_MORPH3.rds\n")

# -----------------------------------------------------------------------------
# STEP 8: Visualization (Top 20 predictors)
# -----------------------------------------------------------------------------
var_imp_top <- var_imp_df %>% slice_head(n = 20)

p_varimp <- ggplot(var_imp_top,
                   aes(x = reorder(variable, MeanDecreaseGini),
                       y = MeanDecreaseGini)) +
  geom_col(fill = "grey") +
  coord_flip() +
  labs(title = "Top 20 Variable Importance (Morphology RF)",
       x = "Predictor",
       y = "MeanDecreaseGini") +
  theme_minimal(base_size = 13)

print(p_varimp)

# -----------------------------------------------------------------------------
# STEP 9: Monte Carlo Evaluation (Multi-class)
# -----------------------------------------------------------------------------
cat("\nRunning Monte Carlo Evaluation (n_iter = 100)...\n")

mc_result <- run_montecarlo_rf_multiclass(
  data = training_df,
  covars = morph_predictors,
  response = response_var,
  n_iter = 100,
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  seed = 42
)

write_csv(mc_result,
          file.path(result_dir, "rf_morph3_mc_results.csv"))

# -----------------------------------------------------------------------------
# STEP 10: CI Summary (Accuracy + Kappa)
# -----------------------------------------------------------------------------
summary_ci <- mc_result %>%
  summarise(
    accuracy_mean = mean(accuracy),
    accuracy_sd   = sd(accuracy),
    kappa_mean    = mean(kappa),
    kappa_sd      = sd(kappa),
    acc_lower     = quantile(accuracy, 0.025),
    acc_upper     = quantile(accuracy, 0.975)
  )

write_csv(summary_ci,
          file.path(result_dir, "rf_morph3_mc_summary_ci.csv"))

cat("\nMonte Carlo CI summary saved.\n")
print(summary_ci)

cat("\nMORPH3 RF classification pipeline complete -- probabilities ready for AGC.\n")