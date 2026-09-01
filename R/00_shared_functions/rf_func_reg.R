# =============================================================================
# rf_func_reg.R
# Shared regression helper functions for the Study 2 proxy-model pipeline
#
# Covers the general-purpose Random Forest regression toolkit used by the
# SPC, AGB, carbon index, and AGC stages: single-fit and Monte Carlo
# evaluation, hyperparameter tuning (mtry, ntree, and full grid search via
# ranger), result summarising, applying a trained model with a fixed
# predictor structure, and per-year performance evaluation from saved
# prediction files.
#
# Requires rf_func_class.R to be sourced first: build_rf_formula() is
# defined there (identical utility, used by both the classification and
# regression toolkits) and is not repeated here to avoid two diverging
# copies of the same function.
#
# Directory handling: result_dir arguments are left as relative paths on
# purpose -- adjust the calling script's working directory or pass an
# explicit path when sourcing this file.
#
# Dependencies: dplyr, randomForest, caret, tibble, purrr, readr (loaded
# below), plus ranger:: and Metrics:: used via explicit namespacing (ranger
# is only required if you call evaluate_rf_grid_ranger_regression()).
#
# Note: run_montecarlo_rf_regression() uses the `%||%` operator (from base R
# 4.4+, or from rlang/purrr on older R versions). Confirm this is available
# in your R version before sourcing.
# =============================================================================

library(dplyr)
library(randomForest)
library(caret)
library(tibble)
library(purrr)
library(readr)


# =============================================================================
# SECTION 1 -- Single-fit and Monte Carlo evaluation
# =============================================================================

# -----------------------------------------------------------------------------
# fit_rf_regression()
# Fit a single Random Forest regression model on one random train/test
# split and return the test-set RMSE.
#
# Input:  data (data frame already in memory), covars (predictor names),
#         response (numeric column name)
# Params: train_frac (split proportion), ntree, mtry (defaults to
#         sqrt(#covars) if not supplied)
# Output: a single RMSE value (numeric)
# -----------------------------------------------------------------------------
fit_rf_regression <- function(seed, data, covars, response,
                              train_frac = 0.7, ntree = 500, mtry = NULL) {
  set.seed(seed)
  train_idx <- sample(seq_len(nrow(data)), size = floor(train_frac * nrow(data)))
  train_data <- data[train_idx, ]
  test_data  <- data[-train_idx, ]
  
  if (is.null(mtry)) mtry <- floor(sqrt(length(covars)))
  formula <- build_rf_formula(covars, response)
  model <- randomForest(formula, data = train_data, ntree = ntree, mtry = mtry)
  
  preds <- predict(model, newdata = test_data)
  actuals <- test_data[[response]]
  return(sqrt(mean((preds - actuals)^2, na.rm = TRUE)))
}

# -----------------------------------------------------------------------------
# run_montecarlo_rf_regression()
# Repeat fit_rf_regression-style evaluation across n_iter random train/test
# splits, returning RMSE/MAE/R2 per iteration (via caret's metric
# functions), suitable for computing CIs afterwards with
# summarize_ci_regression().
#
# Input:  data, covars, response (numeric)
# Params: n_iter (Monte Carlo iterations), ntree, mtry (defaults to
#         sqrt(#covars) if not supplied), nodesize, sample_fraction
#         (fraction of train rows sampled per tree via sampsize), train_frac,
#         seed
# Output: tibble with one row per iteration: iteration, RMSE, MAE, R2
# -----------------------------------------------------------------------------
run_montecarlo_rf_regression <- function(data, covars, response,
                                         n_iter = 100,
                                         ntree = 500,
                                         mtry = NULL,
                                         nodesize = 5,
                                         sample_fraction = 0.8,
                                         train_frac = 0.7,
                                         seed = 42) {
  set.seed(seed)
  seeds <- sample(1:10000, n_iter)
  
  results <- purrr::map_dfr(seeds, function(s) {
    train_idx <- sample(seq_len(nrow(data)), size = floor(train_frac * nrow(data)))
    train <- data[train_idx, ]
    test  <- data[-train_idx, ] %>% filter(if_all(all_of(covars), ~ !is.na(.)))
    
    model <- randomForest(
      x = train[, covars],
      y = train[[response]],
      ntree = ntree,
      mtry = mtry %||% floor(sqrt(length(covars))),
      nodesize = nodesize,
      replace = FALSE,
      sampsize = floor(sample_fraction * nrow(train))
    )
    
    preds <- predict(model, newdata = test[, covars])
    actual <- test[[response]]
    
    tibble::tibble(
      iteration = s,
      RMSE = caret::RMSE(preds, actual),
      MAE  = caret::MAE(preds, actual),
      R2   = caret::R2(preds, actual)
    )
  })
  
  return(results)
}


# =============================================================================
# SECTION 2 -- Hyperparameter tuning
# =============================================================================

# -----------------------------------------------------------------------------
# evaluate_mtry_regression()
# Grid-search tuning over candidate mtry values (number of predictors
# sampled at each split). For each value, repeats n_iter train/test splits
# and reports the median RMSE.
#
# Input:  data, covars, response (numeric)
# Params: mtry_values (defaults to a small candidate set derived from
#         sqrt(#covars), #covars/3, and 1:#covars if not supplied), ntree,
#         n_iter, train_frac
# Output: tibble with one row per mtry value: mtry, rmse (median)
# -----------------------------------------------------------------------------
evaluate_mtry_regression <- function(data, covars, response,
                                     mtry_values = NULL, ntree = 500,
                                     n_iter = 10, train_frac = 0.7) {
  p <- length(covars)
  if (is.null(mtry_values)) mtry_values <- unique(c(floor(sqrt(p)), floor(p / 3), seq(1, p)))
  
  results <- purrr::map_dfr(mtry_values, function(m) {
    metrics <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      set.seed(i + 1000)
      idx <- sample(seq_len(nrow(data)), size = floor(train_frac * nrow(data)))
      train_data <- data[idx, ]
      test_data  <- data[-idx, ]
      
      model <- randomForest(
        formula = build_rf_formula(covars, response),
        data = train_data,
        ntree = ntree,
        mtry = m
      )
      
      preds <- predict(model, newdata = test_data)
      actual <- test_data[[response]]
      metrics[i] <- sqrt(mean((preds - actual)^2, na.rm = TRUE))
    }
    tibble(mtry = m, rmse = median(metrics))
  })
  return(results)
}

# -----------------------------------------------------------------------------
# evaluate_ntree_regression()
# Same idea as evaluate_mtry_regression(), tuning over candidate ntree
# (number of trees) values instead.
#
# Input:  data, covars, response (numeric)
# Params: ntree_values (default seq(100, 1000, 100)), n_iter, train_frac
# Output: tibble with one row per ntree value: ntree, rmse (median)
# -----------------------------------------------------------------------------
evaluate_ntree_regression <- function(data, covars, response,
                                      ntree_values = seq(100, 1000, 100),
                                      n_iter = 10, train_frac = 0.7) {
  results <- purrr::map_dfr(ntree_values, function(n) {
    metrics <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      set.seed(i + 2000)
      idx <- sample(seq_len(nrow(data)), size = floor(train_frac * nrow(data)))
      train_data <- data[idx, ]
      test_data  <- data[-idx, ]
      
      model <- randomForest(
        formula = build_rf_formula(covars, response),
        data = train_data,
        ntree = n
      )
      
      preds <- predict(model, newdata = test_data)
      actual <- test_data[[response]]
      metrics[i] <- sqrt(mean((preds - actual)^2, na.rm = TRUE))
    }
    tibble(ntree = n, rmse = median(metrics))
  })
  return(results)
}

# -----------------------------------------------------------------------------
# evaluate_rf_grid_ranger_regression()
# Full hyperparameter grid search (ntree x mtry x nodesize x
# sample_fraction) for regression, using the faster `ranger` engine instead
# of `randomForest`. For each grid cell, repeats n_iter train/test splits
# and records median/mean/SD RMSE plus a 95% percentile CI.
#
# Input:  data, covars, response (numeric)
# Params: ntree_values, mtry_values (defaults to a small candidate set
#         derived from #covars if not supplied), nodesize_values,
#         sample_fraction_values (each a vector of candidates -- all
#         combinations are tested), n_iter, train_frac, seed
# Output: list(grid_results = full grid with metrics, best_params = the
#         single row with lowest rmse_median)
# Requires the `ranger` package; stops with an error if it isn't installed.
# -----------------------------------------------------------------------------
evaluate_rf_grid_ranger_regression <- function(data, covars, response,
                                               ntree_values = c(300, 500, 700),
                                               mtry_values = NULL,
                                               nodesize_values = c(3, 5, 7),
                                               sample_fraction_values = c(0.6, 0.8, 1.0),
                                               n_iter = 5, train_frac = 0.7, seed = 42) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required. Please install it via install.packages('ranger').")
  }
  
  if (is.null(mtry_values)) {
    p <- length(covars)
    mtry_values <- unique(c(floor(sqrt(p)), floor(p / 3), seq(1, p)))
  }
  
  results <- expand.grid(
    ntree = ntree_values,
    mtry = mtry_values,
    nodesize = nodesize_values,
    sample_fraction = sample_fraction_values,
    stringsAsFactors = FALSE
  )
  
  results$rmse_median <- NA
  results$rmse_mean <- NA
  results$rmse_sd <- NA
  results$ci_lower <- NA
  results$ci_upper <- NA
  
  for (i in seq_len(nrow(results))) {
    row <- results[i, ]
    rmse_vec <- numeric(n_iter)
    for (j in seq_len(n_iter)) {
      set.seed(seed + j)
      idx <- sample(seq_len(nrow(data)), size = floor(train_frac * nrow(data)))
      train_data <- data[idx, ]
      test_data  <- data[-idx, ]
      
      model <- ranger::ranger(
        formula = build_rf_formula(covars, response),
        data = train_data,
        num.trees = row$ntree,
        mtry = row$mtry,
        min.node.size = row$nodesize,
        sample.fraction = row$sample_fraction,
        replace = FALSE,
        importance = "impurity",
        seed = seed + j
      )
      
      preds <- predict(model, data = test_data)$predictions
      actuals <- test_data[[response]]
      rmse_vec[j] <- sqrt(mean((preds - actuals)^2, na.rm = TRUE))
    }
    
    results$rmse_median[i] <- median(rmse_vec)
    results$rmse_mean[i] <- mean(rmse_vec)
    results$rmse_sd[i] <- sd(rmse_vec)
    results$ci_lower[i] <- quantile(rmse_vec, 0.025)
    results$ci_upper[i] <- quantile(rmse_vec, 0.975)
  }
  
  best <- results[which.min(results$rmse_median), ]
  return(list(grid_results = results, best_params = best))
}


# =============================================================================
# SECTION 3 -- Result summarising
# =============================================================================

# -----------------------------------------------------------------------------
# summarize_best_params_regression()
# Pretty-print the best row returned by evaluate_rf_grid_ranger_regression().
#
# Input:  grid_result_list, the list returned by
#         evaluate_rf_grid_ranger_regression() (uses its $best_params
#         element)
# Output: prints a formatted summary to the console; invisibly returns the
#         best_params row
# -----------------------------------------------------------------------------
summarize_best_params_regression <- function(grid_result_list) {
  best <- grid_result_list$best_params
  cat("Best Regression Parameters:\n")
  cat(sprintf("  ntree           : %d\n", best$ntree))
  cat(sprintf("  mtry            : %d\n", best$mtry))
  cat(sprintf("  nodesize        : %d\n", best$nodesize))
  cat(sprintf("  sample_fraction : %.2f\n", best$sample_fraction))
  cat(sprintf("  RMSE (median)   : %.3f\n", best$rmse_median))
  cat(sprintf("  95%% CI          : %.3f - %.3f\n", best$ci_lower, best$ci_upper))
  invisible(best)
}

# -----------------------------------------------------------------------------
# summarize_ci_regression()
# Compute mean, SD, and a 95% percentile CI (2.5%-97.5%) for RMSE, MAE, and
# R2, across Monte Carlo iterations.
#
# Input:  df, the tibble returned by run_montecarlo_rf_regression() (expects
#         columns: RMSE, MAE, R2)
# Output: a single-row summary data frame with mean/sd/lower/upper for each
#         metric
# -----------------------------------------------------------------------------
summarize_ci_regression <- function(df) {
  df %>% summarise(
    RMSE_mean  = mean(RMSE, na.rm = TRUE),
    RMSE_sd    = sd(RMSE, na.rm = TRUE),
    MAE_mean   = mean(MAE, na.rm = TRUE),
    MAE_sd     = sd(MAE, na.rm = TRUE),
    R2_mean    = mean(R2, na.rm = TRUE),
    R2_sd      = sd(R2, na.rm = TRUE),
    
    RMSE_lower = quantile(RMSE, 0.025, na.rm = TRUE),
    RMSE_upper = quantile(RMSE, 0.975, na.rm = TRUE),
    MAE_lower  = quantile(MAE, 0.025, na.rm = TRUE),
    MAE_upper  = quantile(MAE, 0.975, na.rm = TRUE),
    R2_lower   = quantile(R2, 0.025, na.rm = TRUE),
    R2_upper   = quantile(R2, 0.975, na.rm = TRUE)
  )
}


# =============================================================================
# SECTION 4 -- Applying a trained model / per-year performance evaluation
# =============================================================================

# -----------------------------------------------------------------------------
# apply_model_regression()
# Apply an already-trained model to new data, first coercing each predictor
# column to match the type (and, for factors, the exact levels) recorded in
# predictor_structure. This matters because a Random Forest fitted with a
# factor predictor expects the same set of levels at prediction time --
# mismatched levels (e.g. a level present in new data but absent from
# training) will cause predict() to fail or behave unexpectedly.
#
# Input:  df (new data to predict on), predictors (character vector of
#         predictor column names), model (an already-trained model object),
#         predictor_structure (a list with $classes[[var]] and
#         $levels[[var]] recorded from the training data -- typically
#         captured once when the model was originally fitted)
# Output: the predict() result (numeric vector of predictions)
# -----------------------------------------------------------------------------
apply_model_regression <- function(df, predictors, model, predictor_structure) {
  df_pred <- df %>% select(all_of(predictors))
  
  for (var in predictors) {
    cls <- predictor_structure$classes[[var]][1]
    if (cls == "numeric") {
      df_pred[[var]] <- as.numeric(df_pred[[var]])
    } else if (cls == "integer") {
      df_pred[[var]] <- as.integer(df_pred[[var]])
    } else if (cls == "factor") {
      lvls <- predictor_structure$levels[[var]]
      df_pred[[var]] <- factor(df_pred[[var]], levels = lvls)
    } else if (cls == "character") {
      df_pred[[var]] <- as.character(df_pred[[var]])
    }
  }
  
  predict(model, newdata = df_pred)
}

# -----------------------------------------------------------------------------
# eval_regression_performance()
# Read per-year prediction CSVs and compute RMSE, MAE, and R2 for each year.
# Generic across response variable name: for a given `response` (e.g.
# "tSPC"), the function looks for both the `response` column (actual value)
# and a `<response>_pred` column (predicted value) in the same file.
#
# Note: R2 here is computed manually as 1 - SS_res/SS_tot (see the nested
# r_squared() helper), which is not necessarily identical to caret::R2()
# used in run_montecarlo_rf_regression() above -- worth checking both give
# consistent numbers if you report both in the paper.
#
# Input:  <result_dir>/<pattern><year>.csv for each year in `years`
#         expected columns: <response>, <response>_pred
# Params: pattern (filename prefix before the year), response (base column
#         name, default "tSPC")
# Output: tibble with one row per year (for years where the file exists):
#         year, RMSE, MAE, R2
# -----------------------------------------------------------------------------
eval_regression_performance <- function(result_dir, pattern, years, response = "tSPC") {
  
  r_squared <- function(actual, predicted) {
    ss_res <- sum((actual - predicted)^2)
    ss_tot <- sum((actual - mean(actual))^2)
    1 - ss_res / ss_tot
  }
  
  results <- lapply(years, function(yr) {
    path <- file.path(result_dir, paste0(pattern, yr, ".csv"))
    if (!file.exists(path)) return(NULL)
    
    df <- read_csv(path, show_col_types = FALSE) %>%
      filter(!is.na(.data[[response]]), !is.na(.data[[paste0(response, "_pred")]]))
    
    tibble(
      year = yr,
      RMSE = Metrics::rmse(df[[response]], df[[paste0(response, "_pred")]]),
      MAE  = Metrics::mae(df[[response]], df[[paste0(response, "_pred")]]),
      R2   = r_squared(df[[response]], df[[paste0(response, "_pred")]])
    )
  }) %>% bind_rows()
  
  return(results)
}