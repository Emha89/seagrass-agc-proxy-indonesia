# =============================================================================
# rf_crossSite_AGC.R
# Cross-site AGC model evaluation, with lightweight per-site hyperparameter
# tuning: trains a model on each location (and on all locations pooled),
# then tests every trained model against every location, to see how well
# a model generalizes across sites versus staying local.
#
# Pipeline position: 07_validation. Reads training_data_for_GEE_AGC.csv
# (from 06_AGC_final/rf_model_AGC3.R, STEP 5).
#
# Output: agc_cross_site_rmse.csv, agc_delta_rmse_local_vs_allsite.csv,
# agc_delta_rmse_barplot.png/.pdf, written to result/.
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
library(forcats)
library(here)   # install.packages("here") if you don't have it yet

# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
result_dir <- here("result")
input_path <- file.path(result_dir, "training_data_for_GEE_AGC.csv")
min_samples <- 10

# Predictors and response (matches rf_model_AGC3.R exactly)
gse_vars <- paste0("GSE_A", sprintf("%02d", 0:63))
env_vars <- c("depth", "distToLand")
pred_vars <- c("PA_prob", "tSPC_pred", "tAGB_pred", "carbon_index")
morph_pred <- c("P_mixed_short_plus_mono_short", "P_mixed_long",  "P_mono_Ea")

predictors <- c(gse_vars, env_vars, pred_vars, morph_pred)
response_var <- "AGC_simulated"

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------
df <- read_csv(input_path, show_col_types = FALSE) %>%
  filter(!is.na(loc), !is.na(!!sym(response_var)), !!sym(response_var) > 0)

# -----------------------------------------------------------------------------
# REPORT SAMPLE COUNTS BY LOCATION
# -----------------------------------------------------------------------------
cat("\nNumber of samples per location:\n")
df %>%
  count(loc, name = "n_sample") %>%
  arrange(desc(n_sample)) %>%
  print()

# -----------------------------------------------------------------------------
# SPLIT BY LOCATION
# -----------------------------------------------------------------------------
location_groups <- df %>%
  group_by(loc) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= min_samples)

loc_list <- sort(unique(location_groups$loc))
loc_list <- c(loc_list, "AllSite")

cat("\nLocations with sufficient samples (>=", min_samples, "):\n")
print(loc_list)

# -----------------------------------------------------------------------------
# LIGHTWEIGHT TUNING GRID
# -----------------------------------------------------------------------------
mtry_grid <- c(floor(sqrt(length(predictors))), floor(length(predictors)/3), floor(length(predictors)/2))
nodesize_grid <- c(3, 5, 7)

# -----------------------------------------------------------------------------
# BUILD AND EVALUATE MODELS
# -----------------------------------------------------------------------------
results <- list()

for (train_loc in loc_list) {
  cat(sprintf("\nTraining model from: %s\n", train_loc))
  
  train_df <- if (train_loc == "AllSite") df else df %>% filter(loc == train_loc)
  
  train_df <- train_df %>%
    filter(if_all(all_of(predictors), ~ !is.na(.x)))
  
  if (nrow(train_df) < min_samples) {
    cat(sprintf("Skip %s: sample size too small after filtering.\n", train_loc))
    next
  }
  
  # -----------------------------------------------------------------------------
  # TUNING: find the best combination based on RMSE on the training data
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
  
  cat(sprintf("Best params for %s -> mtry: %d, nodesize: %d, RMSE: %.2f\n",
              train_loc, best_params$mtry, best_params$nodesize, best_rmse))
  
  rf_model <- best_model
  
  # -----------------------------------------------------------------------------
  # PREDICT ACROSS ALL TEST SITES
  # -----------------------------------------------------------------------------
  for (test_loc in unique(df$loc)) {
    test_df <- df %>% filter(loc == test_loc)
    
    if (nrow(test_df) == 0) next
    
    test_df <- test_df %>%
      filter(if_all(all_of(predictors), ~ !is.na(.x)))
    
    if (nrow(test_df) == 0) next
    
    pred <- predict(rf_model, newdata = test_df[, predictors])
    rmse_val <- rmse(test_df[[response_var]], pred)
    
    results[[length(results) + 1]] <- tibble(
      ModelSource = train_loc,
      TestSite = test_loc,
      RMSE = rmse_val,
      mtry = best_params$mtry,
      nodesize = best_params$nodesize,
      train_n = nrow(train_df),
      test_n = nrow(test_df)
    )
    
    cat(sprintf("%s -> %s | RMSE: %.2f\n", train_loc, test_loc, rmse_val))
  }
}

# -----------------------------------------------------------------------------
# COMBINE RESULTS
# -----------------------------------------------------------------------------
results_df <- bind_rows(results)

write_csv(results_df, file.path(result_dir, "agc_cross_site_rmse.csv"))

# -----------------------------------------------------------------------------
# PLOT HEATMAP
# -----------------------------------------------------------------------------
heatmap_data <- results_df %>%
  mutate(
    ModelSource = fct_reorder(ModelSource, RMSE, .fun = median),
    TestSite = fct_reorder(TestSite, RMSE, .fun = median)
  )

ggplot(heatmap_data, aes(x = TestSite, y = ModelSource, fill = RMSE)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(RMSE, 1)), size = 4.2) +
  scale_fill_gradient(low = "white", high = "red", name = "RMSE") +
  labs(
    title = "Cross-Site AGC Model RMSE (with lightweight tuning)",
    x = "Test Data",
    y = "Model Source"
  ) +
  theme_minimal(base_size = 14)

cat("\nCross-site evaluation complete. Heatmap displayed and results saved to agc_cross_site_rmse.csv\n")


# =============================================================================
# OPTION 1 -- DELTA-RMSE BAR PLOT (AllSite vs Local)
# delta_RMSE = RMSE(AllSite->site) - RMSE(Local->same site)
# =============================================================================

# -----------------------------------------------------------------------------
# Extract local (diagonal) and AllSite rows
# -----------------------------------------------------------------------------
local_diag <- results_df %>%
  filter(ModelSource == TestSite) %>%
  select(TestSite, RMSE_local = RMSE)

allsite_row <- results_df %>%
  filter(ModelSource == "AllSite") %>%
  select(TestSite, RMSE_allsite = RMSE)

# -----------------------------------------------------------------------------
# Merge and compute delta_RMSE
# -----------------------------------------------------------------------------
delta_df <- local_diag %>%
  inner_join(allsite_row, by = "TestSite") %>%
  mutate(
    delta_RMSE = RMSE_allsite - RMSE_local,
    performance = case_when(
      delta_RMSE > 1  ~ "Local clearly better",
      delta_RMSE < -1 ~ "AllSite clearly better",
      TRUE            ~ "Comparable"
    )
  ) %>%
  arrange(delta_RMSE)

# Save comparison table
write_csv(delta_df, file.path(result_dir, "agc_delta_rmse_local_vs_allsite.csv"))

cat("\ndelta_RMSE table saved: agc_delta_rmse_local_vs_allsite.csv\n")

# -----------------------------------------------------------------------------
# Plot (clean slide-ready style)
# -----------------------------------------------------------------------------
p_delta <- ggplot(delta_df, aes(x = reorder(TestSite, delta_RMSE),
                                y = delta_RMSE,
                                fill = performance)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7) +
  coord_flip() +
  labs(
    title = "General vs Local AGC Model Performance",
    subtitle = expression(Delta ~ RMSE == RMSE[AllSite] - RMSE[Local]),
    x = "Test Site",
    y = expression(Delta ~ RMSE)
  ) +
  scale_fill_manual(values = c(
    "Local clearly better"   = "grey30",
    "AllSite clearly better" = "darkgreen",
    "Comparable"             = "grey70"
  )) +
  theme_minimal(base_size = 15) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_delta)

# Save plot
ggsave(file.path(result_dir, "agc_delta_rmse_barplot.png"),
       p_delta, width = 10, height = 6, dpi = 300)

ggsave(file.path(result_dir, "agc_delta_rmse_barplot.pdf"),
       p_delta, width = 10, height = 6)

cat("delta_RMSE bar plot saved (png + pdf)\n")