# =============================================================================
# 15_agb_model.R
# Final AGB model, trained on the field point estimate (AGB_pred) with the
# settings selected by 14_agb_tuning.R. The script also predicts the AGB
# (tAGB_pred) at the field points of every year.
#
# The fit on the training samples is not an accuracy estimate. Accuracy comes
# from the validation scripts in R/validation/.
#
# Inputs   result/dataAGB_<year>.csv (built by 13_agb_data.R),
#          result/grid_rf_best_AGB.csv, result/predicted_tSPC_<year>.csv
# Outputs  result/final_rf_model_AGB.rds, training_data_for_GEE_AGB.csv,
#          var_importance_rf_AGB.csv, final_model_AGB_summary.csv,
#          predicted_tAGB_<year>.csv
# =============================================================================


# -----------------------------------------------------------------------------
# STEP 0: Libraries and paths
# -----------------------------------------------------------------------------
library(dplyr)
library(readr)
library(purrr)
library(randomForest)
library(ggplot2)

# Paths: open seagrass-agc-proxy-indonesia.Rproj in RStudio (or start R in the
# repository root), so that here() returns the repository folder.
library(here)
result_dir <- here("result")   # input tables and outputs (not distributed, see README)


# -----------------------------------------------------------------------------
# STEP 1: Training rows with a valid point estimate and 95% interval
#   13_agb_data.R builds dataAGB_YYYY (dataPA + predicted SPC) and provides
#   prepare_agb_training_data().
# -----------------------------------------------------------------------------
source(here("R", "workflow", "13_agb_data.R"))

env_vars       <- c("depth", "distToLand")
gse_vars       <- paste0("GSE_A", sprintf("%02d", 0:63))
agb_predictors <- c(env_vars, gse_vars)
response_var   <- "AGB_pred"   # field point estimate

training_df <- prepare_agb_training_data() %>%
  filter(!is.na(AGB_pred), AGB_pred > 0,
         !is.na(AGB_low), !is.na(AGB_up),
         AGB_low >= 0, AGB_up > AGB_low)

training_df_final <- training_df %>%
  filter(!if_any(all_of(c(response_var, agb_predictors)), is.na))

cat(sprintf("AGB training rows: %d (%d removed for missing predictors)\n",
            nrow(training_df_final), nrow(training_df) - nrow(training_df_final)))
if ("loc" %in% names(training_df_final)) print(count(training_df_final, loc))


# -----------------------------------------------------------------------------
# STEP 2: Settings selected by the grid search (14_agb_tuning.R)
# -----------------------------------------------------------------------------
best <- read_csv(file.path(result_dir, "grid_rf_best_AGB.csv"), show_col_types = FALSE)
model_params <- list(ntree           = best$ntree[1],
                     mtry            = best$mtry[1],
                     nodesize        = best$nodesize[1],
                     sample_fraction = best$sample_fraction[1])

cat("AGB hyperparameters (grid_rf_best_AGB.csv):\n")
print(as_tibble(model_params))


# -----------------------------------------------------------------------------
# STEP 3: Final model
# -----------------------------------------------------------------------------
set.seed(42)
final_rf <- randomForest(
  formula    = as.formula(paste(response_var, "~", paste(agb_predictors, collapse = " + "))),
  data       = training_df_final[, c(response_var, agb_predictors)],
  ntree      = model_params$ntree,
  mtry       = model_params$mtry,
  nodesize   = model_params$nodesize,
  sampsize   = floor(model_params$sample_fraction * nrow(training_df_final)),
  replace    = FALSE,
  importance = TRUE
)

saveRDS(final_rf, file.path(result_dir, "final_rf_model_AGB.rds"))
cat("Saved final_rf_model_AGB.rds\n")


# -----------------------------------------------------------------------------
# STEP 4: Training table (target AGB_pred, as in gee/08_agb_model.js)
# -----------------------------------------------------------------------------
training_df_final$tAGB_pred <- predict(final_rf, newdata = training_df_final)

extra_cols   <- c("X", "Y", "location", "loc", "gee_id", "composition_id", "year")
cols_to_save <- c("AGB_pred", "AGB_low", "AGB_up", "AGB_CIwidt", agb_predictors,
                  "tAGB_pred", extra_cols)

write_csv(training_df_final %>% select(any_of(cols_to_save)),
          file.path(result_dir, "training_data_for_GEE_AGB.csv"))
cat("Saved training_data_for_GEE_AGB.csv (target AGB_pred)\n")


# -----------------------------------------------------------------------------
# STEP 5: Training fit and variable importance
#   The training fit is NOT an accuracy estimate. The AGB accuracy comes from
#   R/validation/02_agc_pixel_lolo_transfer.R (PART E).
# -----------------------------------------------------------------------------
obs  <- training_df_final$AGB_pred
fit  <- training_df_final$tAGB_pred
fit_rmse <- sqrt(mean((obs - fit)^2))
fit_r2   <- cor(obs, fit)^2

cat(sprintf("\nAGB training fit (in-sample, not accuracy): RMSE %.2f | R² %.3f\n",
            fit_rmse, fit_r2))

var_imp_df <- as.data.frame(final_rf$importance) %>%
  tibble::rownames_to_column("variable") %>%
  arrange(desc(IncNodePurity)) %>%
  rename(importance = IncNodePurity)
write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_AGB.csv"))

write_csv(tibble(model = "AGB", target = response_var, n_rows = nrow(training_df_final),
                 ntree = model_params$ntree, mtry = model_params$mtry,
                 nodesize = model_params$nodesize,
                 sample_fraction = model_params$sample_fraction,
                 train_fit_RMSE = fit_rmse, train_fit_R2 = fit_r2),
          file.path(result_dir, "final_model_AGB_summary.csv"))

p_varimp <- var_imp_df %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "grey") +
  coord_flip() +
  labs(title = "Top 20 variable importance (AGB, point-estimate target)",
       x = "Variable", y = "IncNodePurity") +
  theme_minimal(base_size = 13)
print(p_varimp)


# -----------------------------------------------------------------------------
# STEP 6: tAGB_pred at the field points of every year
#   Reads predicted_tSPC_YYYY.csv, keeps rows with all predictors, adds
#   tAGB_pred and writes predicted_tAGB_YYYY.csv (input of 16_agc_data.R).
# -----------------------------------------------------------------------------
for (yr in 2017:2024) {
  in_path <- file.path(result_dir, sprintf("predicted_tSPC_%d.csv", yr))
  if (!file.exists(in_path)) {
    cat(sprintf("   %s not found, year skipped\n", basename(in_path)))
    next
  }

  df_yr <- read_csv(in_path, show_col_types = FALSE, guess_max = 100000)

  missing_yr <- setdiff(agb_predictors, names(df_yr))
  if (length(missing_yr) > 0) {
    warning(sprintf("Year %d skipped, missing predictors: %s", yr, paste(missing_yr, collapse = ", ")))
    next
  }

  n_before <- nrow(df_yr)
  df_yr <- df_yr %>%
    mutate(across(all_of(agb_predictors), ~ suppressWarnings(as.numeric(.x)))) %>%
    filter(if_all(all_of(agb_predictors), ~ !is.na(.x)))

  if (nrow(df_yr) == 0) {
    warning(sprintf("Year %d skipped, no rows with all predictors", yr))
    next
  }

  df_yr$tAGB_pred <- predict(final_rf, newdata = as.data.frame(df_yr[, agb_predictors]))

  write_csv(df_yr, file.path(result_dir, sprintf("predicted_tAGB_%d.csv", yr)))
  cat(sprintf("predicted_tAGB_%d.csv: %d rows (%d removed for missing predictors)\n",
              yr, nrow(df_yr), n_before - nrow(df_yr)))
}

cat("\nAGB model done. Next: 16_agc_data.R\n")
