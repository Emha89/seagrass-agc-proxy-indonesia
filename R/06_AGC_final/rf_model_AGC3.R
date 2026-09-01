# =============================================================================
# rf_model_AGC3.R
# Random Forest regression for above-ground carbon (AGC) -- the final
# synthesis stage, with CI-based bootstrap error propagation, grid tuning,
# and Monte Carlo evaluation.
#
# NAMING NOTE: this script's internal header originally read
# "rf_model_AGC.R"; named rf_model_AGC3.R here to match the AGC3 naming
# used across this stage's original file listing. This is a separate
# naming choice from which data-prep CONTENT is used (see below) -- the
# model script's own filename isn't confirmed by any direct reference the
# way the data-prep filename is.
#
# CONFIRMED: this script sources dataPrep_funcGSE_AGC2.R -- that source()
# reference is the direct evidence for which data-prep version is
# actually used (the "check-only" version, not the alternative draft
# that explicitly re-merges carbon_index/morph3). See the header of
# dataPrep_funcGSE_AGC2.R for what that implies about dataPA_<year>.csv
# needing to already be the pass-2 output of dataPrep_func_PA.R.
#
# NOTE (not changed, same pattern as rf_model_AGB2.R): training_df_final
# is built once from the bootstrap set (sim_id == 1) in STEP 1A, then
# immediately rebuilt from scratch in STEP 1C using a different sampling
# approach. The STEP 1A version is therefore never actually used.
#
# NOTE: unlike morphology/SPC/AGB/carbon-index (which each exclude the
# immediately preceding stage's probability/prediction output, to avoid
# predictor redundancy), this final AGC stage deliberately includes ALL
# upstream proxy outputs as predictors (PA_prob, tSPC_pred, tAGB_pred,
# carbon_index, plus the morphology P_* columns) -- expected here, since
# this stage is meant to synthesize the full proxy chain.
#
# Pipeline position: stage 6 of 6 (final). Sources
# dataPrep_funcGSE_AGC3.R for prepare_training_data_agc_noFilter() (this
# script uses the "noFilter" variant, not "withFilter").
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
library(ggpubr)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_reg.R"))            # tuning, formula builder, Monte Carlo, etc.
source(here("R", "00_shared_functions", "rf_func_vis.R"))            # visualization functions
source(here("R", "06_AGC_final", "dataPrep_funcGSE_AGC2.R"))          # for prepare_training_data_agc_noFilter()

result_dir <- here("result")

# -----------------------------------------------------------------------------
# STEP 1: Prepare Training Data (without filtering)
# -----------------------------------------------------------------------------
training_df <- prepare_training_data_agc_noFilter()

# Make sure the data is valid within the CI bounds
training_df <- training_df %>%
  filter(!is.na(AGC_pred), AGC_pred > 0,
         !is.na(AGC_low), !is.na(AGC_up),
         AGC_low >= 0, AGC_up > AGC_low)

# -----------------------------------------------------------------------------
# STEP 1A: Bootstrap Simulation (Empirical Error Propagation)
# -----------------------------------------------------------------------------
n_boot <- 10
set.seed(42)

training_df_boot <- training_df %>%
  slice(rep(1:n(), each = n_boot)) %>%
  mutate(
    sim_id = rep(1:n_boot, times = nrow(training_df)),
    AGC_simulated = runif(n(), AGC_low, AGC_up)
  ) %>%
  filter(!is.na(AGC_simulated))

cat(sprintf("Bootstrap simulated data created: %d rows (n = %d x %d)\n",
            nrow(training_df_boot), nrow(training_df), n_boot))

# Superseded by STEP 1C below -- see the note at the top of this file.
training_df_final <- training_df_boot %>% filter(sim_id == 1)

cat(sprintf("Data for final model: %d rows (sim_id = 1)\n", nrow(training_df_final)))

# -----------------------------------------------------------------------------
# STEP 1B: Normality Check
# -----------------------------------------------------------------------------
p1 <- ggplot(training_df, aes(x = AGC_pred)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "skyblue", color = "white", alpha = 0.8) +
  geom_density(color = "red", linewidth = 1) +
  labs(title = "Distribution of AGC_pred", x = "AGC_pred", y = "Density") +
  theme_minimal()

p2 <- ggplot(training_df, aes(x = AGC_CIwidt)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "lightgreen", color = "white", alpha = 0.8) +
  geom_density(color = "darkgreen", linewidth = 1) +
  labs(title = "Distribution of AGC_CIwidt", x = "CI Width", y = "Density") +
  theme_minimal()

set.seed(42)
shapiro_sample <- training_df %>%
  filter(!is.na(AGC_pred), !is.na(AGC_CIwidt)) %>%
  sample_n(min(nrow(.), 5000))

shapiro_pred <- shapiro.test(shapiro_sample$AGC_pred)
shapiro_ciwd <- shapiro.test(shapiro_sample$AGC_CIwidt)

cat("\nShapiro-Wilk normality test results:\n")
cat(sprintf("AGC_pred  => W = %.4f | p-value = %.4f\n", shapiro_pred$statistic, shapiro_pred$p.value))
cat(sprintf("AGC_CIwidt => W = %.4f | p-value = %.4f\n", shapiro_ciwd$statistic, shapiro_ciwd$p.value))

ggpubr::ggarrange(p1, p2, ncol = 2)

# -----------------------------------------------------------------------------
# STEP 1C: AGC_simulated Simulation (single sampling)
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
# STEP 2: Define Predictors and Response
# -----------------------------------------------------------------------------
env_vars <- c("depth", "distToLand")
gse_vars <- paste0("GSE_A", sprintf("%02d", 0:63))
pred <- c("PA_prob", "tSPC_pred", "tAGB_pred", "carbon_index")
morph_pred <- c("P_mixed_short_plus_mono_short", "P_mixed_long",  "P_mono_Ea")

agc_predictors <- c(gse_vars, env_vars, pred, morph_pred) # alternative ordering considered: c(env_vars, gse_vars, pred, morph_pred)
response_var <- "AGC_simulated"

#------------------------------------------------------------------------------
required_vars <- c(response_var, agc_predictors)

missing_vars <- setdiff(required_vars, names(training_df_final))

if (length(missing_vars) > 0) {
  stop(sprintf(
    "ERROR: Missing required variables in training data: %s",
    paste(missing_vars, collapse = ", ")
  ))
} else {
  cat("All predictor + response variables are available for AGC model.\n")
}
#------------------------------------------------------------------------------

na_count <- sum(is.na(training_df_final[, c(response_var, agc_predictors)]))
cat("Number of NA in training data: ", na_count, "\n")

training_df_final <- training_df_final %>%
  filter(!if_any(all_of(c(response_var, agc_predictors)), is.na))

cat(sprintf("After filtering NA in predictor and response: %d rows remaining\n", nrow(training_df_final)))

# -----------------------------------------------------------------------------
# STEP 3: Grid Tuning
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
    mtry_values = 2:5, # narrower than earlier stages' 4:10, given the larger/more correlated predictor set here
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

best_params <- tune_result$best_params
model_params <- list(
  ntree = best_params$ntree[1],
  mtry = best_params$mtry[1],
  nodesize = best_params$nodesize[1],
  sample_fraction = best_params$sample_fraction[1]
)

# -----------------------------------------------------------------------------
# STEP 4: Train Final RF Model
# -----------------------------------------------------------------------------
rf_formula <- build_rf_formula(agc_predictors, response_var)

final_rf <- randomForest(
  formula = rf_formula,
  data = training_df_final[, c(response_var, agc_predictors)],
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df_final)),
  replace = FALSE,
  importance = TRUE
)

saveRDS(final_rf, file.path(result_dir, "final_rf_model_AGC.rds"))

# -----------------------------------------------------------------------------
# STEP 5: Predict on Training Data & Export for GEE
# -----------------------------------------------------------------------------
training_df_final$tAGC_pred <- predict(final_rf, newdata = training_df_final)

extra_cols <- c("X", "Y", "location", 'loc', "gee_id", "composition_id", "year")
cols_to_save <- c("AGC_pred", "AGC_CIwidt", "AGC_simulated", agc_predictors,
                  "tAGC_pred", intersect(extra_cols, colnames(training_df_final)))

write_csv(training_df_final[, cols_to_save],
          file.path(result_dir, "training_data_for_GEE_AGC.csv"))

# -----------------------------------------------------------------------------
# STEP 6: Save Predictor Structure
# -----------------------------------------------------------------------------
predictor_structure <- list(
  classes = sapply(training_df_final[agc_predictors], class),
  levels = lapply(training_df_final[agc_predictors], function(x) if (is.factor(x)) levels(x) else NULL)
)

saveRDS(predictor_structure, file.path(result_dir, "rf_predictor_structure_AGC.rds"))
cat("Saved rf_predictor_structure_AGC.rds\n")

# -----------------------------------------------------------------------------
# STEP 7: Model Performance
# -----------------------------------------------------------------------------
cat("\n==============================\n")
cat("Final AGC Model Performance (Training Data)\n")
cat("==============================\n")
cat(sprintf("RMSE: %.3f | MAE: %.3f | R2: %.3f\n",
            RMSE(training_df_final$tAGC_pred, training_df_final$AGC_pred),
            MAE(training_df_final$tAGC_pred, training_df_final$AGC_pred),
            R2(training_df_final$tAGC_pred, training_df_final$AGC_pred)))

# -----------------------------------------------------------------------------
# STEP 8: Variable Importance
# -----------------------------------------------------------------------------
if (!is.null(final_rf$importance)) {
  var_imp_df <- as.data.frame(final_rf$importance) %>%
    tibble::rownames_to_column("variable") %>%
    arrange(desc(IncNodePurity)) %>%
    rename(importance = IncNodePurity)
  
  write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_AGC.csv"))
  cat("Variable importance saved: var_importance_rf_AGC.csv\n")
} else {
  warning("Variable importance not found in model object.")
}

# -----------------------------------------------------------------------------
# STEP 9: Visualization
# -----------------------------------------------------------------------------
if (file.exists(file.path(result_dir, "var_importance_rf_AGC.csv"))) {
  var_imp_top <- read_csv(file.path(result_dir, "var_importance_rf_AGC.csv"), show_col_types = FALSE) %>%
    arrange(desc(importance)) %>% slice_head(n = 20)
  p_varimp <- ggplot(var_imp_top, aes(x = reorder(variable, importance), y = importance)) +
    geom_col(fill = "grey") +
    coord_flip() +
    labs(title = "Top 20 Variable Importance (AGC)", x = "Variable", y = "IncNodePurity") +
    theme_minimal(base_size = 13)
  print(p_varimp)
}

# -----------------------------------------------------------------------------
# STEP 10: Monte Carlo Evaluation
# -----------------------------------------------------------------------------
training_df_boot <- training_df_boot %>%
  filter(!if_any(all_of(c(response_var, agc_predictors)), is.na))

cat(sprintf("After filtering NA in bootstrap data: %d rows remaining\n", nrow(training_df_boot)))
cat("\nRunning Monte Carlo Evaluation (n_iter = 100, using bootstrap data)...\n")

mc_result <- run_montecarlo_rf_regression(
  data = training_df_boot,
  covars = agc_predictors,
  response = response_var,
  n_iter = 100,
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  nodesize = model_params$nodesize,
  sample_fraction = model_params$sample_fraction,
  train_frac = 0.7,
  seed = 42
)

write_csv(mc_result, file.path(result_dir, "rf_agc_mc_results.csv"))

summary_ci <- summarize_ci_regression(mc_result)
write_csv(summary_ci, file.path(result_dir, "rf_agc_mc_summary_ci.csv"))

cat("\nMonte Carlo CI summary saved to rf_agc_mc_summary_ci.csv\n")
print(summary_ci)