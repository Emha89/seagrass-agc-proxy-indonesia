# =============================================================================
# 18_agc_model.R
# Final AGC model of the proxy workflow, trained on the field point estimate
# (AGC_pred). Predictors are the 64 GSE bands, depth, distance to land and
# the seven proxy layers. The carbon index is the PREDICTED carbon index
# (step 12), the same layer that is used for mapping in GEE. The settings
# are those of Table S3 of the paper (grid search in 17_agc_tuning.R).
#
# The fit on the training samples is not an accuracy estimate. Accuracy comes
# from the validation scripts in R/validation/.
#
# Inputs   result/dataAGC_<year>.csv (built by 16_agc_data.R),
#          result/predicted_indexC_<year>.csv, result/predicted_tAGB_<year>.csv,
#          result/dataPA_<year>.csv
# Outputs  result/final_rf_model_AGC.rds, training_data_for_GEE_AGC.csv,
#          var_importance_rf_AGC.csv (Figure S3a), final_model_AGC_summary.csv,
#          predicted_tAGC_<year>.csv
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
# STEP 1: Build dataAGC_YYYY
#   16_agc_data.R builds dataAGC_YYYY from dataPA and predicted_tAGB, and
#   provides prepare_agc_training_data().
# -----------------------------------------------------------------------------
source(here("R", "workflow", "16_agc_data.R"))

# -----------------------------------------------------------------------------
# STEP 2: Put the PREDICTED carbon index into dataAGC_YYYY, then load the
#         training rows
#   dataPA_YYYY holds the field carbon index, which 16_agc_data.R copies into
#   dataAGC_YYYY. The AGC model uses the predicted carbon index instead:
#     carbon_index_obs = field-derived value (kept for checks only)
#     carbon_index     = carbon_index_pred from predicted_indexC_YYYY.csv,
#                        matched by gee_id within the same year
#   Running it twice gives the same result.
# -----------------------------------------------------------------------------
for (y in 2017:2023) {
  agc_path <- file.path(result_dir, sprintf("dataAGC_%d.csv", y))
  ci_path  <- file.path(result_dir, sprintf("predicted_indexC_%d.csv", y))
  if (!file.exists(agc_path)) next
  if (!file.exists(ci_path)) stop("", basename(ci_path), " not found, needed for ", basename(agc_path))

  a <- read_csv(agc_path, show_col_types = FALSE, guess_max = 100000)

  ci <- read_csv(ci_path, show_col_types = FALSE, guess_max = 100000) %>%
    transmute(gee_id_key  = as.character(gee_id),
              ci_pred_new = suppressWarnings(as.numeric(carbon_index_pred))) %>%
    distinct(gee_id_key, .keep_all = TRUE)

  # Keep the field value once (only if it was not saved before)
  if (!"carbon_index_obs" %in% names(a)) a$carbon_index_obs <- a$carbon_index

  a <- a %>%
    mutate(gee_id_key = as.character(gee_id)) %>%
    left_join(ci, by = "gee_id_key") %>%
    mutate(carbon_index = ci_pred_new) %>%
    select(-gee_id_key, -ci_pred_new, -any_of("carbon_index_pred"))

  write_csv(a, agc_path)
  cat(sprintf("   dataAGC_%d.csv: carbon_index = predicted (%d rows, %d without a predicted value)\n",
              y, nrow(a), sum(is.na(a$carbon_index))))
}

env_vars   <- c("depth", "distToLand")
gse_vars   <- paste0("GSE_A", sprintf("%02d", 0:63))
proxy_vars <- c("PA_prob", "tSPC_pred", "tAGB_pred", "carbon_index")
morph_pred <- c("P_mixed_short_plus_mono_short", "P_mixed_long", "P_mono_Ea")

agc_predictors <- c(gse_vars, env_vars, proxy_vars, morph_pred)
response_var   <- "AGC_pred"   # field point estimate

training_df <- prepare_agc_training_data() %>%
  filter(!is.na(AGC_pred), AGC_pred > 0,
         !is.na(AGC_low), !is.na(AGC_up),
         AGC_low >= 0, AGC_up > AGC_low)

missing_vars <- setdiff(c(response_var, agc_predictors), names(training_df))
if (length(missing_vars) > 0) {
  stop("Missing variables in dataAGC: ", paste(missing_vars, collapse = ", "))
}

training_df_final <- training_df %>%
  filter(!if_any(all_of(c(response_var, agc_predictors)), is.na))

cat(sprintf("AGC training rows: %d (%d removed for missing predictors)\n",
            nrow(training_df_final), nrow(training_df) - nrow(training_df_final)))
if ("loc" %in% names(training_df_final)) print(count(training_df_final, loc))


# -----------------------------------------------------------------------------
# STEP 3: Check that carbon_index is the PREDICTED carbon index
#   Each training value is matched (gee_id + value) against the predicted
#   carbon index files and against the field carbon index in dataPA.
# -----------------------------------------------------------------------------
read_pairs <- function(pattern, col) {
  map_dfr(2017:2023, function(y) {
    p <- file.path(result_dir, sprintf(pattern, y))
    if (!file.exists(p)) return(NULL)
    read_csv(p, show_col_types = FALSE, guess_max = 100000) %>%
      transmute(key = paste(as.character(gee_id),
                            round(suppressWarnings(as.numeric(.data[[col]])), 6)))
  }) %>% distinct()
}

ci_pred_keys <- read_pairs("predicted_indexC_%d.csv", "carbon_index_pred")
ci_obs_keys  <- read_pairs("dataPA_%d.csv", "carbon_index")
train_keys   <- paste(as.character(training_df_final$gee_id),
                      round(training_df_final$carbon_index, 6))

share_pred <- mean(train_keys %in% ci_pred_keys$key)
share_obs  <- mean(train_keys %in% ci_obs_keys$key)
cat(sprintf("carbon_index in the training rows: %.1f%% match the predicted values, %.1f%% match the field values\n",
            100 * share_pred, 100 * share_obs))

if (share_pred < 0.95) {
  stop("carbon_index is not the predicted carbon index. ",
       "Check STEP 2 of this script and the predicted_indexC_YYYY.csv files.")
}


# -----------------------------------------------------------------------------
# STEP 4: Settings (Table S3 of the paper)
# -----------------------------------------------------------------------------
model_params <- list(ntree = 900, mtry = 5, nodesize = 9, sample_fraction = 0.6)

grid_path <- file.path(result_dir, "grid_rf_best_AGC.csv")
if (file.exists(grid_path)) {
  g <- read_csv(grid_path, show_col_types = FALSE)
  same <- isTRUE(all.equal(c(g$ntree[1], g$mtry[1], g$nodesize[1], g$sample_fraction[1]),
                           unlist(model_params), check.attributes = FALSE))
  cat(sprintf("AGC hyperparameters: Table S3 values %s grid_rf_best_AGC.csv\n",
              if (same) "match" else "DIFFER FROM"))
  if (!same) print(g)
}
print(as_tibble(model_params))


# -----------------------------------------------------------------------------
# STEP 5: Final model
# -----------------------------------------------------------------------------
set.seed(42)
final_rf <- randomForest(
  formula    = as.formula(paste(response_var, "~", paste(agc_predictors, collapse = " + "))),
  data       = training_df_final[, c(response_var, agc_predictors)],
  ntree      = model_params$ntree,
  mtry       = model_params$mtry,
  nodesize   = model_params$nodesize,
  sampsize   = floor(model_params$sample_fraction * nrow(training_df_final)),
  replace    = FALSE,
  importance = TRUE
)

saveRDS(final_rf, file.path(result_dir, "final_rf_model_AGC.rds"))
cat("Saved final_rf_model_AGC.rds\n")


# -----------------------------------------------------------------------------
# STEP 6: Training table (target AGC_pred, as in gee/10_agc_model.js)
# -----------------------------------------------------------------------------
training_df_final$tAGC_pred <- predict(final_rf, newdata = training_df_final)

extra_cols   <- c("X", "Y", "location", "loc", "gee_id", "composition_id", "year")
cols_to_save <- c("AGC_pred", "AGC_low", "AGC_up", "AGC_CIwidt", agc_predictors,
                  "tAGC_pred", extra_cols)

write_csv(training_df_final %>% select(any_of(cols_to_save)),
          file.path(result_dir, "training_data_for_GEE_AGC.csv"))
cat("Saved training_data_for_GEE_AGC.csv (target AGC_pred)\n")


# -----------------------------------------------------------------------------
# STEP 7: Training fit, field uncertainty and variable importance
#   The training fit is NOT an accuracy estimate (see R/validation/).
#   SD_data = median of (AGC_up - AGC_low) / (2 x 1.96) over the training
#   rows, the typical uncertainty of a reference AGC value.
# -----------------------------------------------------------------------------
obs      <- training_df_final$AGC_pred
fit      <- training_df_final$tAGC_pred
fit_rmse <- sqrt(mean((obs - fit)^2))
fit_r2   <- cor(obs, fit)^2
sd_data  <- median((training_df_final$AGC_up - training_df_final$AGC_low) / (2 * 1.96))

cat(sprintf("\nAGC training fit (in-sample, not accuracy): RMSE %.2f | R² %.3f\n",
            fit_rmse, fit_r2))
cat(sprintf("SD_data (median SD of the reference AGC): %.3f gC m-2\n", sd_data))

var_imp_df <- as.data.frame(final_rf$importance) %>%
  tibble::rownames_to_column("variable") %>%
  arrange(desc(IncNodePurity)) %>%
  rename(importance = IncNodePurity)
write_csv(var_imp_df, file.path(result_dir, "var_importance_rf_AGC.csv"))

write_csv(tibble(model = "AGC", target = response_var, n_rows = nrow(training_df_final),
                 ntree = model_params$ntree, mtry = model_params$mtry,
                 nodesize = model_params$nodesize,
                 sample_fraction = model_params$sample_fraction,
                 train_fit_RMSE = fit_rmse, train_fit_R2 = fit_r2,
                 SD_data = sd_data,
                 carbon_index_match_predicted_pct = 100 * share_pred),
          file.path(result_dir, "final_model_AGC_summary.csv"))

p_varimp <- var_imp_df %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "grey") +
  coord_flip() +
  labs(title = "Top 20 variable importance (AGC, point-estimate target)",
       x = "Variable", y = "IncNodePurity") +
  theme_minimal(base_size = 13)
print(p_varimp)


# -----------------------------------------------------------------------------
# STEP 8: tAGC_pred at the field points of every year
#   Starts from predicted_tAGB_YYYY.csv, takes the predicted carbon index (and
#   the P_* columns if they are missing) from predicted_indexC_YYYY.csv, adds
#   tAGC_pred and writes predicted_tAGC_YYYY.csv. These files are not needed
#   for GEE.
# -----------------------------------------------------------------------------
for (yr in 2017:2024) {
  tagb_path <- file.path(result_dir, sprintf("predicted_tAGB_%d.csv", yr))
  idx_path  <- file.path(result_dir, sprintf("predicted_indexC_%d.csv", yr))
  if (!file.exists(tagb_path) || !file.exists(idx_path)) {
    cat(sprintf("   Year %d skipped (predicted_tAGB or predicted_indexC not found)\n", yr))
    next
  }

  df_yr <- read_csv(tagb_path, show_col_types = FALSE, guess_max = 100000) %>%
    select(-any_of(c("carbon_index", "carbon_index_pred")))

  need_morph <- setdiff(morph_pred, names(df_yr))
  idx <- read_csv(idx_path, show_col_types = FALSE, guess_max = 100000) %>%
    select(gee_id, carbon_index_pred, any_of(need_morph)) %>%
    distinct(gee_id, .keep_all = TRUE)

  df_yr <- df_yr %>%
    left_join(idx, by = "gee_id") %>%
    rename(carbon_index = carbon_index_pred)

  missing_yr <- setdiff(agc_predictors, names(df_yr))
  if (length(missing_yr) > 0) {
    warning(sprintf("Year %d skipped, missing predictors: %s", yr, paste(missing_yr, collapse = ", ")))
    next
  }

  n_before <- nrow(df_yr)
  df_yr <- df_yr %>%
    mutate(across(all_of(agc_predictors), ~ suppressWarnings(as.numeric(.x)))) %>%
    filter(if_all(all_of(agc_predictors), ~ !is.na(.x)))

  if (nrow(df_yr) == 0) {
    warning(sprintf("Year %d skipped, no rows with all predictors", yr))
    next
  }

  df_yr$tAGC_pred <- predict(final_rf, newdata = as.data.frame(df_yr[, agc_predictors]))

  write_csv(df_yr, file.path(result_dir, sprintf("predicted_tAGC_%d.csv", yr)))
  cat(sprintf("predicted_tAGC_%d.csv: %d rows (%d removed for missing predictors)\n",
              yr, nrow(df_yr), n_before - nrow(df_yr)))
}

cat("\nAGC model done. Next: the scripts in R/validation/\n")
