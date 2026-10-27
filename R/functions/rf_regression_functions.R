# =============================================================================
# rf_regression_functions.R
# Random Forest helper functions for the regression models (percent cover,
# carbon index, AGB and AGC): grid search with ranger, and prediction with a
# fixed predictor structure. Sources rf_classification_functions.R for
# build_rf_formula().
# =============================================================================

library(dplyr)
library(randomForest)
library(caret)
library(tibble)
library(purrr)
library(readr)

# build_rf_formula() is defined in the classification functions
source(here::here("R", "functions", "rf_classification_functions.R"))


# =============================================================================
# SECTION 1 -- Hyperparameter tuning
# =============================================================================

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
# SECTION 2 -- Applying a trained model
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

