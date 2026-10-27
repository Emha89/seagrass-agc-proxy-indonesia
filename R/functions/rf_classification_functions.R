# =============================================================================
# rf_classification_functions.R
# Random Forest helper functions for the classification models (presence and
# leaf morphology): training data, formula, grid search with ranger, threshold
# evaluation. build_rf_formula() is also used by the regression models.
#
# Default folders point to result/ in the repository (here::here("result")).
# =============================================================================

library(dplyr)
library(randomForest)
library(caret)
library(tibble)
library(purrr)


# =============================================================================
# SECTION 1 -- Data preparation and formula building
# =============================================================================

# -----------------------------------------------------------------------------
# prepare_training_data()
# Combine per-year PA (occurrence) data into a single data frame, with
# optional class balancing.
#
# Input:  <input_dir>/dataPA_<year>.csv for each year in `years`
#         expected column: PA (0/1)
# Params: balance = FALSE (default) keeps the data as-is; TRUE downsamples
#         the majority class to match the minority class count (only
#         triggers if the PA=0 and PA=1 counts actually differ)
# Output: combined data frame, with an added `year` factor column
# -----------------------------------------------------------------------------
prepare_training_data <- function(
    input_dir = here::here("result"),
    years = c("2018", "2019", "2020", "2021", "2022", "2023"),
    balance = FALSE,
    seed = 42
) {
  set.seed(seed)
  
  # Combine data across years, if the file exists
  df_all <- purrr::map_dfr(years, function(yr) {
    path <- file.path(input_dir, paste0("dataPA_", yr, ".csv"))
    if (file.exists(path)) {
      readr::read_csv(path, show_col_types = FALSE) %>%
        dplyr::mutate(year = as.factor(yr))
    } else {
      warning(sprintf("File not found: %s", path))
      NULL
    }
  })
  
  # Keep only valid PA values and convert to factor
  df_all <- df_all %>%
    dplyr::filter(PA %in% c(0, 1)) %>%
    dplyr::mutate(PA = as.factor(PA))
  
  # Basic distribution info
  n_total <- nrow(df_all)
  n_pa1 <- sum(df_all$PA == 1, na.rm = TRUE)
  n_pa0 <- sum(df_all$PA == 0, na.rm = TRUE)
  
  cat(sprintf("Training data loaded: %d rows | PA=1: %d | PA=0: %d\n",
              n_total, n_pa1, n_pa0))
  
  # If balancing was not requested
  if (!balance) {
    cat("Using the data as-is (no balancing).\n")
    return(df_all)
  }
  
  # -------------------------------------
  # DYNAMIC BALANCING LOGIC
  # -------------------------------------
  if (n_pa0 > n_pa1) {
    cat("Balancing: downsampling PA=0 to match the PA=1 count.\n")
    df_balanced <- df_all %>%
      dplyr::group_by(PA) %>%
      dplyr::mutate(n_class = dplyr::n()) %>%
      dplyr::group_split() %>%
      {
        df_0 <- .[[which(levels(df_all$PA) == "0")]]
        df_1 <- .[[which(levels(df_all$PA) == "1")]]
        df_0_down <- df_0 %>% dplyr::sample_n(size = nrow(df_1), replace = FALSE)
        dplyr::bind_rows(df_0_down, df_1)
      } %>%
      dplyr::ungroup()
    
    cat(sprintf("Balancing complete: PA=1 and PA=0 each have %d samples (total %d)\n",
                n_pa1, nrow(df_balanced)))
    return(df_balanced)
    
  } else {
    cat("Distribution is already balanced, no balancing needed.\n")
    return(df_all)
  }
}

# -----------------------------------------------------------------------------
# build_rf_formula()
# Small utility: builds an R formula object ("response ~ covar1 + covar2 +
# ...") from a vector of covariate names and a response name. Used by every
# model-fitting function below.
#
# Input:  covars (character vector of predictor names), response (single
#         column name string)
# Output: a formula object; stops with an error if covars is empty/NULL
# -----------------------------------------------------------------------------
build_rf_formula <- function(covars, response) {
  if (length(covars) == 0 || is.null(covars)) stop("Covariate list is empty or NULL.")
  as.formula(paste(response, "~", paste(covars, collapse = " + ")))
}


# =============================================================================
# SECTION 2 -- Binary classification functions (e.g. PA / occurrence)
# =============================================================================

# -----------------------------------------------------------------------------
# evaluate_thresholds_classification()
# For each candidate probability threshold, repeat n_iter random train/test
# splits and RF fits, reporting mean and SD accuracy at that threshold.
# Used to choose an operating threshold (see plot_pa_threshold_accuracy in
# plot_functions.R for the companion plot).
#
# Input:  data, covars, response (binary 0/1), thresholds (numeric vector
#         of candidate cutoffs to test)
# Params: n_iter, ntree, mtry
# Output: tibble with one row per threshold: threshold, accuracy, accuracy_sd
# -----------------------------------------------------------------------------
evaluate_thresholds_classification <- function(data, covars, response, 
                                               thresholds, n_iter = 10, 
                                               ntree = 500, mtry = NULL) {
  results <- purrr::map_dfr(thresholds, function(thresh) {
    acc_vec <- numeric(n_iter)
    
    for (i in seq_len(n_iter)) {
      set.seed(i + 100)
      
      index <- caret::createDataPartition(data[[response]], p = 0.7, list = FALSE)
      train_data <- data[index, ]
      test_data  <- data[-index, ]
      
      model <- randomForest(
        formula = build_rf_formula(covars, response),
        data = train_data,
        ntree = ntree,
        mtry = mtry,  # mtry must be passed explicitly here
        importance = FALSE
      )
      
      probs <- predict(model, newdata = test_data, type = "prob")[, "1"]
      preds <- ifelse(probs > thresh, 1, 0)
      actuals <- as.numeric(as.character(test_data[[response]]))
      
      acc_vec[i] <- mean(preds == actuals, na.rm = TRUE)
    }
    
    tibble(
      threshold = thresh,
      accuracy = mean(acc_vec),
      accuracy_sd = sd(acc_vec)
    )
  })
  
  return(results)
}

# -----------------------------------------------------------------------------
# evaluate_rf_grid_ranger_classification()
# Full hyperparameter grid search (ntree x mtry x nodesize x
# sample_fraction) for binary classification, using the faster `ranger`
# engine instead of `randomForest`. For each grid cell, repeats n_iter
# train/test splits and records mean/SD accuracy plus a 95% percentile CI.
#
# Input:  data, covars, response (binary 0/1)
# Params: ntree_values, mtry_values, nodesize_values, sample_fraction_values
#         (each a vector of candidates -- all combinations are tested),
#         threshold, n_iter, train_frac, seed
# Output: list(grid_results = full grid with metrics, best_params = the
#         single row with highest accuracy_mean)
# Requires the `ranger` package; stops with an error if it isn't installed.
# -----------------------------------------------------------------------------
evaluate_rf_grid_ranger_classification <- function(data, covars, response,
                                                   ntree_values, mtry_values,
                                                   nodesize_values, sample_fraction_values,
                                                   threshold = 0.7, n_iter = 5,
                                                   train_frac = 0.7, seed = 42) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required for this function. Please install it.")
  }
  
  set.seed(seed)
  grid_results <- expand.grid(
    ntree = ntree_values,
    mtry = mtry_values,
    nodesize = nodesize_values,
    sample_fraction = sample_fraction_values,
    stringsAsFactors = FALSE
  )
  
  grid_results$accuracy_mean <- NA
  grid_results$accuracy_sd   <- NA
  grid_results$ci_lower <- NA
  grid_results$ci_upper <- NA
  
  for (i in seq_len(nrow(grid_results))) {
    ntree <- grid_results$ntree[i]
    mtry  <- grid_results$mtry[i]
    nodesize <- grid_results$nodesize[i]
    sampfrac <- grid_results$sample_fraction[i]
    
    acc_iter <- c()
    for (iter in 1:n_iter) {
      index <- caret::createDataPartition(data[[response]], p = train_frac, list = FALSE)
      train_data <- data[index, ]
      test_data  <- data[-index, ]
      
      model <- ranger::ranger(
        formula = build_rf_formula(covars, response),
        data = train_data,
        num.trees = ntree,
        mtry = mtry,
        min.node.size = nodesize,
        sample.fraction = sampfrac,
        probability = TRUE,
        classification = TRUE,
        seed = seed + iter
      )
      
      probs <- predict(model, data = test_data)$predictions[, "1"]
      preds <- ifelse(probs > threshold, 1, 0)
      actual <- as.numeric(as.character(test_data[[response]]))
      acc_iter <- c(acc_iter, mean(preds == actual, na.rm = TRUE))
    }
    
    grid_results$accuracy_mean[i] <- mean(acc_iter)
    grid_results$accuracy_sd[i]   <- sd(acc_iter)
    grid_results$ci_lower[i] <- quantile(acc_iter, 0.025)
    grid_results$ci_upper[i] <- quantile(acc_iter, 0.975)
  }
  
  best_idx <- which.max(grid_results$accuracy_mean)
  best_params <- grid_results[best_idx, ]
  
  return(list(grid_results = grid_results, best_params = best_params))
}

# =============================================================================
# SECTION 3 -- Multi-class classification functions (morphology, 3-class)
# These mirror the binary-classification functions in Section 2 above, but
# handle a multi-level factor response (e.g. leaf morphology) instead of a
# 0/1 response. They do not affect the PA binary-classification functions.
# =============================================================================

# -----------------------------------------------------------------------------
# evaluate_rf_grid_ranger_multiclass()
# Multi-class counterpart to evaluate_rf_grid_ranger_classification(): full
# hyperparameter grid search (ntree x mtry x nodesize x sample_fraction)
# using `ranger`, assigning each test case to the class with the highest
# predicted probability, and recording mean/SD accuracy plus a 95%
# percentile CI per grid cell.
#
# Input:  data, covars, response (multi-level factor/character column)
# Params: ntree_values, mtry_values, nodesize_values, sample_fraction_values
#         (each a vector of candidates -- all combinations are tested),
#         n_iter, train_frac, seed
# Output: list(grid_results = full grid with metrics, best_params = the
#         single row with highest accuracy_mean)
# Requires the `ranger` package; stops with an error if it isn't installed.
# -----------------------------------------------------------------------------
evaluate_rf_grid_ranger_multiclass <- function(
    data, covars, response,
    ntree_values, mtry_values,
    nodesize_values, sample_fraction_values,
    n_iter = 5,
    train_frac = 0.7,
    seed = 42
) {
  
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required. Please install it.")
  }
  
  set.seed(seed)
  
  # Grid of hyperparameter combinations
  grid_results <- expand.grid(
    ntree = ntree_values,
    mtry = mtry_values,
    nodesize = nodesize_values,
    sample_fraction = sample_fraction_values,
    stringsAsFactors = FALSE
  )
  
  # Storage for metrics
  grid_results$accuracy_mean <- NA
  grid_results$accuracy_sd   <- NA
  grid_results$ci_lower <- NA
  grid_results$ci_upper <- NA
  
  # ------------------------------------------------------------
  # Loop over all parameter combinations
  # ------------------------------------------------------------
  for (i in seq_len(nrow(grid_results))) {
    
    ntree    <- grid_results$ntree[i]
    mtry     <- grid_results$mtry[i]
    nodesize <- grid_results$nodesize[i]
    sampfrac <- grid_results$sample_fraction[i]
    
    acc_iter <- c()
    
    for (iter in 1:n_iter) {
      
      # Train/test split
      idx <- caret::createDataPartition(
        data[[response]],
        p = train_frac,
        list = FALSE
      )
      
      train_data <- data[idx, ]
      test_data  <- data[-idx, ]
      
      # Train ranger RF in probability mode
      model <- ranger::ranger(
        formula = build_rf_formula(covars, response),
        data = train_data,
        num.trees = ntree,
        mtry = mtry,
        min.node.size = nodesize,
        sample.fraction = sampfrac,
        probability = TRUE,
        classification = TRUE,
        seed = seed + iter
      )
      
      # --------------------------------------------------------
      # Multi-class prediction: predictions is a probability
      # matrix (n rows x K classes)
      # --------------------------------------------------------
      probs <- predict(model, data = test_data)$predictions
      
      # Predicted class = highest probability column
      pred_class <- colnames(probs)[max.col(probs)]
      
      # True class labels
      actual_class <- as.character(test_data[[response]])
      
      # Accuracy
      acc_iter <- c(acc_iter,
                    mean(pred_class == actual_class, na.rm = TRUE))
    }
    
    # Store tuning results
    grid_results$accuracy_mean[i] <- mean(acc_iter)
    grid_results$accuracy_sd[i]   <- sd(acc_iter)
    grid_results$ci_lower[i]      <- quantile(acc_iter, 0.025)
    grid_results$ci_upper[i]      <- quantile(acc_iter, 0.975)
  }
  
  # Best parameters (highest mean accuracy)
  best_idx <- which.max(grid_results$accuracy_mean)
  best_params <- grid_results[best_idx, ]
  
  return(list(
    grid_results = grid_results,
    best_params  = best_params
  ))
}
