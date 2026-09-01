# =============================================================================
# rf_model_AGB2.R
# Random Forest regression for above-ground biomass (AGB), with CI-based
# bootstrap error propagation, grid tuning, and Monte Carlo evaluation.
#
# NAMING NOTE: this script's internal header originally read
# "rf_model_AGB.R", but the file it sources (originally referenced as
# "dataPrep_funcGSE_AGB2.R") only matches the on-disk name established
# from the original file listing, not the header text used when that
# script was pasted here ("dataPrep_func_AGB.R"). Following the same
# precedent, this file is named rf_model_AGB2.R (matching the original
# file listing) rather than rf_model_AGB.R -- flag if this guess is wrong.
#
# FIX APPLIED (STEP 4): sampsize was computed from nrow(training_df), but
# the model is fit on training_df_final, which has additional NA filtering
# (STEP 2) not applied to training_df. With replace = FALSE, a sampsize
# larger than the actual fitted data would error. Changed to
# nrow(training_df_final).
#
# NOTE (not changed): training_df_final is built once in STEP 1
# (sim_id == 1 from the bootstrap set), then immediately rebuilt from
# scratch in STEP 1C using a different sampling approach. The STEP 1
# version is therefore never actually used -- STEP 1C's version is what
# flows into tuning and training.
#
# Pipeline position: stage 4 of 6. Sources dataPrep_funcGSE_AGB2.R for
# prepare_training_data_agb_noFilter() (this script uses the "noFilter"
# variant, not "withFilter").
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
source(here("R", "04_AGB", "dataPrep_funcGSE_AGB2.R"))                # for prepare_training_data_agb_noFilter()

result_dir <- here("result")

# -----------------------------------------------------------------------------
# STEP 1: Prepare Training Data with CI-Based Bootstrapping (Error Propagation)
# -----------------------------------------------------------------------------
training_df <- prepare_training_data_agb_noFilter()

# Filter to valid data and set up CI bounds from the AGB_low and AGB_up columns
training_df <- training_df %>%
  filter(!is.na(AGB_pred), AGB_pred > 0,
         !is.na(AGB_low), !is.na(AGB_up),
         AGB_low >= 0, AGB_up > AGB_low)

# -----------------------------------------------------------------------------
# Generate bootstrap simulations (empirical error propagation)
# -----------------------------------------------------------------------------
n_boot <- 10  # number of bootstrap simulations per location
set.seed(42)

# Expand each row n_boot times, then randomly sample AGB within the CI range
training_df_boot <- training_df %>%
  slice(rep(1:n(), each = n_boot)) %>%
  mutate(
    sim_id = rep(1:n_boot, times = nrow(training_df)),
    AGB_simulated = runif(n(), AGB_low, AGB_up)
  ) %>%
  filter(!is.na(AGB_simulated))

cat(sprintf("Bootstrap simulated data created: %d rows (n = %d x %d)\n",
            nrow(training_df_boot), nrow(training_df), n_boot))

# Select one simulation (sim_id = 1) for the final model (STEP 4) --
# NOTE: this is superseded by STEP 1C below, see the note at the top of
# this file.
training_df_final <- training_df_boot %>% filter(sim_id == 1)

cat(sprintf("Data for final model: %d rows (sim_id = 1)\n", nrow(training_df_final)))

# -----------------------------------------------------------------------------
# STEP 1B: Normality Check (AGB_pred & AGB_CIwidt)
# -----------------------------------------------------------------------------

# Plot the AGB_pred and AGB_CIwidt distributions BEFORE simulation
p1 <- ggplot(training_df, aes(x = AGB_pred)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "skyblue", color = "white", alpha = 0.8) +
  geom_density(color = "red", linewidth = 1) +
  labs(title = "Distribution of AGB_pred", x = "AGB_pred", y = "Density") +
  theme_minimal()

p2 <- ggplot(training_df, aes(x = AGB_CIwidt)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "lightgreen", color = "white", alpha = 0.8) +
  geom_density(color = "darkgreen", linewidth = 1) +
  labs(title = "Distribution of AGB_CIwidt", x = "CI Width", y = "Density") +
  theme_minimal()

# Shapiro-Wilk test
set.seed(42)
shapiro_sample <- training_df %>%
  filter(!is.na(AGB_pred), !is.na(AGB_CIwidt)) %>%
  sample_n(min(nrow(.), 5000))

shapiro_pred <- shapiro.test(shapiro_sample$AGB_pred)
shapiro_ciwd <- shapiro.test(shapiro_sample$AGB_CIwidt)

cat("\nShapiro-Wilk normality test results:\n")
cat(sprintf("AGB_pred  => W = %.4f | p-value = %.4f\n", shapiro_pred$statistic, shapiro_pred$p.value))
cat(sprintf("AGB_CIwidt => W = %.4f | p-value = %.4f\n", shapiro_ciwd$statistic, shapiro_ciwd$p.value))

# Display the plot
ggpubr::ggarrange(p1, p2, ncol = 2)

# -----------------------------------------------------------------------------
# STEP 1C: AGB_simulated Simulation via Empirical Bootstrap
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
pred <- c("PA_prob", "tSPC_pred") # , "carbon_index"

# wave_vars and pred are deliberately excluded, to avoid predictor
# redundancy across proxy stages (same principle applied in morphology and SPC)
agb_predictors <- c(env_vars, gse_vars)
response_var <- "AGB_simulated"

# Check NA after response & predictor are defined
na_count <- sum(is.na(training_df_final[, c(response_var, agb_predictors)]))
cat("Number of NA in training data: ", na_count, "\n")

colSums(is.na(training_df_final[, c(response_var, agb_predictors)]))

# Filter NA rows
training_df_final <- training_df_final %>%
  filter(!if_any(all_of(c(response_var, agb_predictors)), is.na))

cat(sprintf("After filtering NA in predictor and response: %d rows remaining\n", nrow(training_df_final)))

# -----------------------------------------------------------------------------
# STEP 3: Grid Tuning
# -----------------------------------------------------------------------------
grid_rds <- file.path(result_dir, "grid_rf_tuning_AGB.rds")

tune_result <- if (file.exists(grid_rds)) {
  cat("Loading cached tuning results...\n")
  readRDS(grid_rds)
} else {
  cat("Running tuning...\n")
  res <- evaluate_rf_grid_ranger_regression(
    data = training_df_final,     # use the final data with AGB_simulated
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
rf_formula <- build_rf_formula(agb_predictors, response_var)

final_rf <- randomForest(
  formula = rf_formula,
  data = training_df_final[, c(response_var, agb_predictors)],
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  nodesize = model_params$nodesize,
  sampsize = floor(model_params$sample_fraction * nrow(training_df_final)),
  replace = FALSE,
  importance = TRUE
)

saveRDS(final_rf, file.path(result_dir, "final_rf_model_AGB.rds"))

# -----------------------------------------------------------------------------
# STEP 5: Predict on Training Data & Export for GEE
# -----------------------------------------------------------------------------
training_df_final$tAGB_pred <- predict(final_rf, newdata = training_df_final)

extra_cols <- c("X", "Y", "location", "gee_id", "composition_id", "year")
cols_to_save <- c("AGB_pred", "AGB_CIwidt", "AGB_simulated", agb_predictors,
                  "tAGB_pred", intersect(extra_cols, colnames(training_df)))

write_csv(training_df_final[, cols_to_save],
          file.path(result_dir, "training_data_for_GEE_AGB.csv"))

# -----------------------------------------------------------------------------
# STEP 5b: Save Predictor Structure for Apply Model
# -----------------------------------------------------------------------------
predictor_structure <- list(
  classes = sapply(training_df[agb_predictors], class),
  levels = lapply(training_df[agb_predictors], function(x) if (is.factor(x)) levels(x) else NULL)
)

saveRDS(predictor_structure, file.path(result_dir, "rf_predictor_structure_AGB.rds"))
cat("Saved rf_predictor_structure_AGB.rds\n")

# -----------------------------------------------------------------------------
# STEP 6: Final Model Performance
# -----------------------------------------------------------------------------
cat("\n==============================\n")
cat("Final AGB Model Performance (Training Data)\n")
cat("==============================\n")
cat(sprintf("RMSE: %.3f | MAE: %.3f | R2: %.3f\n",
            RMSE(training_df_final$tAGB_pred, training_df_final$AGB_pred),
            MAE(training_df_final$tAGB_pred, training_df_final$AGB_pred),
            R2(training_df_final$tAGB_pred, training_df_final$AGB_pred)))

# -----------------------------------------------------------------------------
# STEP 7: Variable Importance
# -----------------------------------------------------------------------------
if (!is.null(final_rf$importance)) {
  var_imp_df <- as.data.frame(final_rf$importance) %>%
    tibble::rownames_to_column("variable") %>%
    arrange(desc(IncNodePurity)) %>%
    rename(importance = IncNodePurity)
  
  write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_AGB.csv"))
  
  cat("Variable importance saved: var_importance_rf_AGB.csv\n")
} else {
  warning("Variable importance not found in model object.")
}

# -----------------------------------------------------------------------------
# STEP 8: Visualization
# -----------------------------------------------------------------------------
if (file.exists(file.path(result_dir, "var_importance_rf_AGB.csv"))) {
  var_imp_top <- read_csv(file.path(result_dir, "var_importance_rf_AGB.csv"), show_col_types = FALSE) %>%
    arrange(desc(importance)) %>% slice_head(n = 20)
  p_varimp <- ggplot(var_imp_top, aes(x = reorder(variable, importance), y = importance)) +
    geom_col(fill = "grey") +
    coord_flip() +
    labs(title = "Top 20 Variable Importance (AGB)", x = "Variable", y = "IncNodePurity") +
    theme_minimal(base_size = 13)
  print(p_varimp)
}

# -----------------------------------------------------------------------------
# STEP 8b: Filter NA on Bootstrapped Data Before Monte Carlo
# -----------------------------------------------------------------------------
training_df_boot <- training_df_boot %>%
  filter(!if_any(all_of(c(response_var, agb_predictors)), is.na))

cat(sprintf("After filtering NA in bootstrap data: %d rows remaining\n", nrow(training_df_boot)))

# -----------------------------------------------------------------------------
# STEP 9: Monte Carlo Evaluation
# -----------------------------------------------------------------------------
cat("\nRunning Monte Carlo Evaluation (n_iter = 100, using bootstrap data)...\n")

mc_result <- run_montecarlo_rf_regression(
  data = training_df_boot,
  covars = agb_predictors,
  response = response_var,
  n_iter = 100,
  ntree = model_params$ntree,
  mtry = model_params$mtry,
  nodesize = model_params$nodesize,
  sample_fraction = model_params$sample_fraction,
  train_frac = 0.7,
  seed = 42
)

write_csv(mc_result, file.path(result_dir, "rf_agb_mc_results.csv"))

# -----------------------------------------------------------------------------
# STEP 10: CI Summary
# -----------------------------------------------------------------------------
summary_ci <- summarize_ci_regression(mc_result)
write_csv(summary_ci, file.path(result_dir, "rf_agb_mc_summary_ci.csv"))

cat("\nMonte Carlo CI summary saved to rf_agb_mc_summary_ci.csv\n")
print(summary_ci)