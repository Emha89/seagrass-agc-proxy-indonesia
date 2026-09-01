# =============================================================================
# rf_func_class.R
# Shared classification helper functions for the Study 2 proxy-model pipeline
#
# Covers the general-purpose Random Forest classification toolkit used by
# both the PA (occurrence, binary) and morphology (3-class) stages: data
# preparation, formula building, single-fit and Monte Carlo evaluation,
# threshold/hyperparameter tuning (mtry, ntree, and full grid search via
# ranger), and result summarising.
#
# Directory handling: input_dir defaults are left as relative paths
# ("../result") on purpose -- adjust the calling script's working directory
# or pass an explicit input_dir when sourcing this file.
#
# Dependencies: dplyr, randomForest, caret, tibble, purrr (loaded below),
# plus readr:: and ranger:: used via explicit namespacing (ranger is only
# required if you call one of the *_ranger_* grid-search functions).
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
    input_dir = "../result",
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
# fit_rf_classification()
# Fit a single Random Forest classifier on one random train/test split and
# return the classification accuracy at a given probability threshold.
#
# Input:  data (data frame already in memory), covars (predictor names),
#         response (binary 0/1 column name)
# Params: train_frac (split proportion), ntree, threshold (probability
#         cutoff for the positive class), mtry (defaults to sqrt(#covars)
#         if not supplied)
# Output: a single accuracy value (numeric)
# -----------------------------------------------------------------------------
fit_rf_classification <- function(seed, data, covars, response,
                                  train_frac = 0.7, ntree = 500, threshold = 0.7, mtry = NULL) {
  set.seed(seed)
  train_idx <- caret::createDataPartition(data[[response]], p = train_frac, list = FALSE)
  train_data <- data[train_idx, ]
  test_data  <- data[-train_idx, ]
  
  train_data[[response]] <- as.factor(train_data[[response]])
  test_data[[response]]  <- as.factor(test_data[[response]])
  
  if (is.null(mtry)) mtry <- floor(sqrt(length(covars)))
  formula <- build_rf_formula(covars, response)
  model <- randomForest(formula, data = train_data, ntree = ntree, mtry = mtry)
  
  probs <- predict(model, newdata = test_data, type = "prob")[, "1"]
  preds <- ifelse(probs > threshold, 1, 0)
  actuals <- as.numeric(as.character(test_data[[response]]))
  return(mean(preds == actuals))
}

# -----------------------------------------------------------------------------
# run_montecarlo_rf_classification()
# Repeat fit_rf_classification-style evaluation across n_iter random
# train/test splits, returning full accuracy/precision/recall/F1 per
# iteration (via caret::confusionMatrix), suitable for computing CIs
# afterwards with summarize_ci_classification().
#
# Input:  data, covars, response (binary 0/1)
# Params: n_iter (Monte Carlo iterations), train_frac, ntree, mtry,
#         threshold (default 0.5 here -- note this differs from the 0.7
#         default used in the other classification functions below, worth
#         double-checking this is intentional), seed
# Output: tibble with one row per iteration: iteration, accuracy, precision,
#         recall, f1
# -----------------------------------------------------------------------------
run_montecarlo_rf_classification <- function(data, covars, response,
                                             n_iter = 100, train_frac = 0.7,
                                             ntree = 500, mtry = NULL,
                                             threshold = 0.5, seed = 42) {
  set.seed(seed)
  seeds <- sample(1:10000, n_iter)
  
  results <- purrr::map_dfr(seeds, function(s) {
    train_idx <- caret::createDataPartition(data[[response]], p = train_frac, list = FALSE)
    train_data <- data[train_idx, ]
    test_data  <- data[-train_idx, ] %>% filter(if_all(all_of(covars), ~ !is.na(.)))
    
    model <- randomForest(
      x = train_data[, covars],
      y = as.factor(train_data[[response]]),
      ntree = ntree,
      mtry = mtry,
      replace = FALSE
    )
    
    actual <- factor(test_data[[response]], levels = c(0, 1))
    probs  <- predict(model, newdata = test_data[, covars], type = "prob")[, "1"]
    preds  <- ifelse(probs >= threshold, 1, 0)
    preds  <- factor(preds, levels = c(0, 1))
    
    cm <- caret::confusionMatrix(preds, actual, positive = "1")
    tibble::tibble(
      iteration = s,
      accuracy  = as.numeric(cm$overall["Accuracy"]),
      precision = as.numeric(cm$byClass["Precision"]),
      recall    = as.numeric(cm$byClass["Recall"]),
      f1        = as.numeric(cm$byClass["F1"])
    )
  })
  return(results)
}

# -----------------------------------------------------------------------------
# evaluate_thresholds_classification()
# For each candidate probability threshold, repeat n_iter random train/test
# splits and RF fits, reporting mean and SD accuracy at that threshold.
# Used to choose an operating threshold (see plot_pa_threshold_accuracy in
# rf_func_vis.R for the companion plot).
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
# evaluate_mtry_classification()
# Grid-search tuning over candidate mtry values (number of predictors
# sampled at each split). For each value, repeats n_iter train/test splits
# and reports the median accuracy.
#
# Input:  data, covars, response (binary 0/1)
# Params: mtry_values (defaults to a small candidate set derived from
#         sqrt(#covars), #covars/3, and 1:#covars if not supplied), ntree,
#         threshold, n_iter, train_frac
# Output: tibble with one row per mtry value: mtry, accuracy (median)
# -----------------------------------------------------------------------------
evaluate_mtry_classification <- function(data, covars, response,
                                         mtry_values = NULL, ntree = 500,
                                         threshold = 0.7, n_iter = 10, train_frac = 0.7) {
  p <- length(covars)
  if (is.null(mtry_values)) mtry_values <- unique(c(floor(sqrt(p)), floor(p / 3), seq(1, p)))
  
  results <- purrr::map_dfr(mtry_values, function(m) {
    metrics <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      set.seed(i + 1000)
      idx <- caret::createDataPartition(data[[response]], p = train_frac, list = FALSE)
      train_data <- data[idx, ]
      test_data  <- data[-idx, ]
      
      train_data[[response]] <- as.factor(train_data[[response]])
      test_data[[response]]  <- as.factor(test_data[[response]])
      
      model <- randomForest(
        formula = build_rf_formula(covars, response),
        data = train_data,
        ntree = ntree,
        mtry = m
      )
      
      probs <- predict(model, newdata = test_data, type = "prob")[, "1"]
      preds <- ifelse(probs > threshold, 1, 0)
      actual <- as.numeric(as.character(test_data[[response]]))
      metrics[i] <- mean(preds == actual, na.rm = TRUE)
    }
    tibble(mtry = m, accuracy = median(metrics))
  })
  return(results)
}

# -----------------------------------------------------------------------------
# evaluate_ntree_classification()
# Same idea as evaluate_mtry_classification(), tuning over candidate ntree
# (number of trees) values instead.
#
# Input:  data, covars, response (binary 0/1)
# Params: ntree_values (default seq(100, 1000, 100)), threshold, n_iter,
#         train_frac
# Output: tibble with one row per ntree value: ntree, accuracy (median)
# -----------------------------------------------------------------------------
evaluate_ntree_classification <- function(data, covars, response,
                                          ntree_values = seq(100, 1000, 100),
                                          threshold = 0.7, n_iter = 10, train_frac = 0.7) {
  results <- purrr::map_dfr(ntree_values, function(n) {
    metrics <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      set.seed(i + 2000)
      idx <- caret::createDataPartition(data[[response]], p = train_frac, list = FALSE)
      train_data <- data[idx, ]
      test_data  <- data[-idx, ]
      
      train_data[[response]] <- as.factor(train_data[[response]])
      test_data[[response]]  <- as.factor(test_data[[response]])
      
      model <- randomForest(
        formula = build_rf_formula(covars, response),
        data = train_data,
        ntree = n
      )
      
      probs <- predict(model, newdata = test_data, type = "prob")[, "1"]
      preds <- ifelse(probs > threshold, 1, 0)
      actual <- as.numeric(as.character(test_data[[response]]))
      metrics[i] <- mean(preds == actual, na.rm = TRUE)
    }
    tibble(ntree = n, accuracy = median(metrics))
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

# -----------------------------------------------------------------------------
# summarize_best_params_classification()
# Pretty-print the best row returned by
# evaluate_rf_grid_ranger_classification() (or the multiclass version below).
#
# Input:  grid_result_list, the list returned by one of the grid-search
#         functions above (uses its $best_params element)
# Output: prints a formatted summary to the console; invisibly returns the
#         best_params row
# -----------------------------------------------------------------------------
summarize_best_params_classification <- function(grid_result_list) {
  best <- grid_result_list$best_params
  cat("Best Classification Parameters:\n")
  cat(sprintf("  ntree           : %d\n", best$ntree))
  cat(sprintf("  mtry            : %d\n", best$mtry))
  cat(sprintf("  nodesize        : %d\n", best$nodesize))
  cat(sprintf("  sample_fraction : %.2f\n", best$sample_fraction))
  cat(sprintf("  Accuracy (mean) : %.3f\n", best$accuracy_mean))
  cat(sprintf("  95%% CI          : %.3f - %.3f\n", best$ci_lower, best$ci_upper))
  invisible(best)
}

# -----------------------------------------------------------------------------
# summarize_ci_classification()
# Compute mean, SD, and a 95% percentile CI (2.5%-97.5%) for accuracy,
# precision, recall, and F1, across Monte Carlo iterations.
#
# Input:  df, the tibble returned by run_montecarlo_rf_classification()
#         (expects columns: accuracy, precision, recall, f1)
# Output: a single-row summary data frame with mean/sd/lower/upper for
#         each metric
# -----------------------------------------------------------------------------
summarize_ci_classification <- function(df) {
  df %>% summarise(
    accuracy_mean  = mean(accuracy, na.rm = TRUE),
    accuracy_sd    = sd(accuracy, na.rm = TRUE),
    precision_mean = mean(precision, na.rm = TRUE),
    precision_sd   = sd(precision, na.rm = TRUE),
    recall_mean    = mean(recall, na.rm = TRUE),
    recall_sd      = sd(recall, na.rm = TRUE),
    f1_mean        = mean(f1, na.rm = TRUE),
    f1_sd          = sd(f1, na.rm = TRUE),
    
    accuracy_lower  = quantile(accuracy, 0.025, na.rm = TRUE),
    accuracy_upper  = quantile(accuracy, 0.975, na.rm = TRUE),
    precision_lower = quantile(precision, 0.025, na.rm = TRUE),
    precision_upper = quantile(precision, 0.975, na.rm = TRUE),
    recall_lower    = quantile(recall, 0.025, na.rm = TRUE),
    recall_upper    = quantile(recall, 0.975, na.rm = TRUE),
    f1_lower        = quantile(f1, 0.025, na.rm = TRUE),
    f1_upper        = quantile(f1, 0.975, na.rm = TRUE)
  )
}


# =============================================================================
# SECTION 3 -- Multi-class classification functions (morphology, 3-class)
# These mirror the binary-classification functions in Section 2 above, but
# handle a multi-level factor response (e.g. leaf morphology) instead of a
# 0/1 response. They do not affect the PA binary-classification functions.
# =============================================================================

# -----------------------------------------------------------------------------
# run_montecarlo_rf_multiclass()
# Multi-class counterpart to run_montecarlo_rf_classification(): repeats
# n_iter random train/test splits, fitting a Random Forest classifier on a
# multi-level factor response, and reports overall accuracy and kappa per
# iteration.
#
# Input:  data, covars, response (multi-level factor/character column,
#         e.g. morph3)
# Params: n_iter, train_frac, ntree, mtry, seed
# Output: tibble with one row per iteration: iteration, accuracy, kappa
# -----------------------------------------------------------------------------
run_montecarlo_rf_multiclass <- function(data, covars, response,
                                         n_iter = 100, train_frac = 0.7,
                                         ntree = 500, mtry = NULL, seed = 42) {
  set.seed(seed)
  seeds <- sample(1:10000, n_iter)
  
  results <- purrr::map_dfr(seeds, function(s) {
    idx <- caret::createDataPartition(data[[response]], p = train_frac, list = FALSE)
    train_data <- data[idx, ]
    test_data  <- data[-idx, ] %>% filter(if_all(all_of(covars), ~ !is.na(.)))
    
    model <- randomForest(
      x = train_data[, covars],
      y = as.factor(train_data[[response]]),
      ntree = ntree,
      mtry = mtry,
      replace = FALSE
    )
    
    preds <- predict(model, newdata = test_data[, covars])
    actual <- factor(test_data[[response]])
    
    preds <- factor(preds, levels = levels(actual))  # force matching factor levels
    cm <- caret::confusionMatrix(preds, actual)
    
    tibble::tibble(
      iteration = s,
      accuracy  = as.numeric(cm$overall["Accuracy"]),
      kappa     = as.numeric(cm$overall["Kappa"])
    )
  })
  
  return(results)
}

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