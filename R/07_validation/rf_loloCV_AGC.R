# =============================================================================
# rf_lolocv_AGC_full.R
# Leave-one-location-out cross-validation (LOLO-CV) for the AGC model,
# with lightweight per-fold hyperparameter tuning, and visualization.
#
# Pipeline position: 07_validation. Reads training_data_for_GEE_AGC.csv
# (from 06_AGC_final/rf_model_AGC3.R, STEP 5).
#
# OPEN ITEM: Visualization 2's y-axis is hardcoded to
# scale_y_continuous(limits = c(10, 20), ...) -- if actual RMSE values
# fall outside that range, points would be silently clipped from the
# plot. Worth checking against the real results, or switching to a
# data-driven range (e.g. removing the explicit limits).
#
# Output: agc_lolocv_rmse.csv, written to result/.
#
# Data availability: this script's input is not included in this
# repository.
# =============================================================================

library(dplyr)
library(readr)
library(randomForest)
library(Metrics)
library(tidyr)
library(ggplot2)
library(caret)
library(here)   # install.packages("here") if you don't have it yet

rmse <- Metrics::rmse
mae  <- Metrics::mae

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
result_dir <- here("result")
input_path <- file.path(result_dir, "training_data_for_GEE_AGC.csv")
output_csv <- file.path(result_dir, "agc_lolocv_rmse.csv")
min_samples <- 10

# Response and predictor variables (matches rf_model_AGC3.R exactly)
gse_vars <- paste0("GSE_A", sprintf("%02d", 0:63))
env_vars <- c("depth", "distToLand")
pred_vars <- c("PA_prob", "tSPC_pred", "tAGB_pred", "carbon_index")
morph_pred <- c("P_mixed_short_plus_mono_short", "P_mixed_long",  "P_mono_Ea")

predictors <- c(gse_vars, env_vars, pred_vars, morph_pred)

response_var <- "AGC_simulated"

# -----------------------------------------------------------------------------
# LOAD & PREPARE DATA
# -----------------------------------------------------------------------------
df <- read_csv(input_path, show_col_types = FALSE) %>%
  filter(!is.na(loc), !is.na(!!sym(response_var)), !!sym(response_var) > 0)

valid_locs <- df %>%
  count(loc) %>%
  filter(n >= min_samples) %>%
  pull(loc) %>%
  sort()

cat("Valid locations for LOLO-CV:\n")
print(valid_locs)

# -----------------------------------------------------------------------------
# LIGHTWEIGHT TUNING GRID
# -----------------------------------------------------------------------------
mtry_grid <- c(floor(sqrt(length(predictors))), floor(length(predictors) / 3), floor(length(predictors) / 2))
nodesize_grid <- c(3, 5, 7)

# -----------------------------------------------------------------------------
# LOOP: LEAVE-ONE-LOCATION-OUT CROSS VALIDATION
# -----------------------------------------------------------------------------
results <- list()

for (test_loc in valid_locs) {
  cat(sprintf("\nTest on: %s (excluded from training)\n", test_loc))
  
  train_df <- df %>% filter(loc != test_loc)
  test_df  <- df %>% filter(loc == test_loc)
  
  train_df <- train_df %>%
    filter(if_all(all_of(predictors), ~ !is.na(.x)))
  test_df <- test_df %>%
    filter(if_all(all_of(predictors), ~ !is.na(.x)))
  
  if (nrow(train_df) == 0 || nrow(test_df) == 0) {
    warning(sprintf("Skipping %s due to insufficient data.\n", test_loc))
    next
  }
  
  # -----------------------------------------------------------------------------
  # TUNING: find the best mtry and nodesize combination on the training set
  # -----------------------------------------------------------------------------
  best_model <- NULL
  best_rmse <- Inf
  best_params <- list()
  
  for (mtry_val in mtry_grid) {
    for (nodesize_val in nodesize_grid) {
      model_try <- randomForest(
        x = train_df[, predictors],
        y = train_df[[response_var]],
        ntree = 500,
        mtry = mtry_val,
        nodesize = nodesize_val
      )
      pred_try <- predict(model_try, newdata = train_df[, predictors])
      rmse_try <- rmse(train_df[[response_var]], pred_try)
      
      if (rmse_try < best_rmse) {
        best_model <- model_try
        best_rmse <- rmse_try
        best_params <- list(mtry = mtry_val, nodesize = nodesize_val)
      }
    }
  }
  
  rf_model <- best_model
  
  # -----------------------------------------------------------------------------
  # EVALUATE MODEL AT TEST LOCATION
  # -----------------------------------------------------------------------------
  pred <- predict(rf_model, newdata = test_df[, predictors])
  
  results[[length(results) + 1]] <- tibble(
    TestLocation = test_loc,
    n_train = nrow(train_df),
    n_test = nrow(test_df),
    RMSE = rmse(test_df[[response_var]], pred),
    MAE  = mae(test_df[[response_var]], pred),
    R2   = R2(pred, test_df[[response_var]]),
    mtry = best_params$mtry,
    nodesize = best_params$nodesize
  )
  
  cat(sprintf("%s | Train: %d | Test: %d | RMSE: %.2f | MAE: %.2f | R2: %.2f | mtry: %d | nodesize: %d\n",
              test_loc, nrow(train_df), nrow(test_df),
              results[[length(results)]]$RMSE,
              results[[length(results)]]$MAE,
              results[[length(results)]]$R2,
              best_params$mtry, best_params$nodesize))
}

results_df <- bind_rows(results)

# -----------------------------------------------------------------------------
# SAVE RESULT TABLE
# -----------------------------------------------------------------------------
write_csv(results_df, output_csv)
cat(sprintf("\nLOLO-CV results saved to: %s\n", output_csv))
print(results_df)

# -----------------------------------------------------------------------------
# VISUALIZATION 1: Barplot RMSE per Location
# -----------------------------------------------------------------------------
p1 <- results_df %>%
  arrange(desc(RMSE)) %>%
  mutate(TestLocation = factor(TestLocation, levels = TestLocation)) %>%
  ggplot(aes(x = TestLocation, y = RMSE)) +
  geom_col(fill = "grey") +
  geom_text(aes(label = sprintf("%.2f", RMSE)), vjust = -0.4, size = 3.8) +
  labs(
    title = "LOLO-CV RMSE per Test Location (AGC)",
    subtitle = "Model trained on all other locations",
    x = "Test Location",
    y = "RMSE (AGC_simulated)"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

print(p1)

# -----------------------------------------------------------------------------
# VISUALIZATION 2: RMSE vs R2 Scatter Plot
# -----------------------------------------------------------------------------
p2 <- results_df %>%
  ggplot(aes(x = R2, y = RMSE, label = TestLocation)) +
  geom_point(size = 4, color = "grey") +
  geom_text(vjust = -0.7, hjust = 0.5, size = 3.5) +
  labs(
    title = "LOLO-CV AGC Model: RMSE vs R2",
    x = expression(R^2 ~ "(Prediction Accuracy)"),
    y = "RMSE (Prediction Error)"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(10, 20), breaks = seq(0, 20, by = 2)) +
  theme_minimal(base_size = 14)

print(p2)