# =============================================================================
# rf_func_vis.R
# Shared visualization functions for the Study 2 proxy-model pipeline
#
# Covers plotting functions for all six proxy stages:
#   - PA        seagrass occurrence probability
#   - MORPH3    leaf morphology probability (3-class)
#   - SPC/PCT   percent cover
#   - AGB       above-ground biomass
#   - AGC       above-ground carbon (final stage)
#
# Most PA / MORPH3 / SPC functions read pre-saved prediction CSVs from
# result_dir (default "../result"); most AGB / AGC functions instead take an
# already-loaded data frame (df) as their argument. This mixed pattern
# matches how each stage's evaluation script currently calls these
# functions and has not been changed here.
#
# Directory handling: result_dir defaults are left as relative paths
# ("../result") on purpose -- adjust the calling script's working directory
# or pass an explicit result_dir when sourcing this file.
#
# Dependencies: dplyr, readr, ggplot2, tidyr (loaded below), plus purrr::,
# tibble::, Metrics::, caret::, and tidyselect:: used via explicit
# namespacing in a few functions below (no need to library() those too).
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)


# =============================================================================
# SECTION 1 -- PA (occurrence probability) visualization functions
# =============================================================================

# -----------------------------------------------------------------------------
# plot_pa_performance_by_year()
# Grouped bar chart of PA (occurrence) model accuracy, precision, recall,
# and F1 score for each survey year.
#
# Input:  <result_dir>/pa_performance_by_year.csv
#         expected columns: year, accuracy, precision, recall, f1_score
# Output: ggplot bar chart, one group of bars per year
# -----------------------------------------------------------------------------
plot_pa_performance_by_year <- function(result_dir = "../result") {
  perf_df <- read_csv(file.path(result_dir, "pa_performance_by_year.csv"), show_col_types = FALSE)
  
  metrics_df <- perf_df %>%
    select(year, accuracy, precision, recall, f1_score) %>%
    pivot_longer(cols = c(accuracy, precision, recall, f1_score),
                 names_to = "metric", values_to = "value") %>%
    mutate(year_label = factor(year, levels = unique(year)))
  
  ggplot(metrics_df, aes(x = year_label, y = value, fill = metric)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = round(value, 2)), position = position_dodge(0.9), vjust = -0.5) +
    labs(title = "PA Model Performance per Year", x = "Year", y = "Metric Value") +
    scale_fill_manual(values = c(
      "accuracy" = "gray20",
      "precision" = "gray90",
      "recall" = "gray70",
      "f1_score" = "gray50"
    )) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_pa_prob_distribution()
# Histogram of predicted PA probability, split by the actual (true) class,
# pooled across all available years.
#
# Input:  <result_dir>/predicted_PA_<year>.csv for year in 2017:2023
#         expected columns: PA (0/1), PA_prob
# Output: ggplot overlapping histogram, coloured by actual class
# -----------------------------------------------------------------------------
plot_pa_prob_distribution <- function(result_dir = "../result") {
  years <- 2017:2023
  files <- paste0(result_dir, "/predicted_PA_", years, ".csv")
  prob_df <- lapply(files, function(path) {
    if (file.exists(path)) {
      read_csv(path, show_col_types = FALSE) %>%
        filter(!is.na(PA), !is.na(PA_prob)) %>%
        mutate(year = gsub("\\D", "", basename(path)))
    } else {
      NULL
    }
  }) %>% bind_rows()
  
  prob_df$PA <- factor(prob_df$PA, levels = c(0, 1), labels = c("Absent", "Present"))
  
  ggplot(prob_df, aes(x = PA_prob, fill = PA)) +
    geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
    labs(title = "Distribution of Predicted Probability by Actual Class",
         x = "Predicted Probability (PA_prob)", y = "Count") +
    scale_fill_manual(values = c("Absent" = "#999999", "Present" = "#1B9E77")) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_pa_spatial_map()
# Point map of predicted PA probability across Indonesia, pooled across all
# available years.
#
# Input:  <result_dir>/predicted_PA_<year>.csv for year in 2017:2023
#         expected columns: PA, PA_prob, lon, lat
# Output: ggplot map (yellow = low probability, dark green = high probability)
# -----------------------------------------------------------------------------
plot_pa_spatial_map <- function(result_dir = "../result") {
  years <- 2017:2023
  files <- paste0(result_dir, "/predicted_PA_", years, ".csv")
  prob_df <- lapply(files, function(path) {
    if (file.exists(path)) {
      read_csv(path, show_col_types = FALSE) %>%
        filter(!is.na(PA), !is.na(PA_prob)) %>%
        mutate(year = gsub("\\D", "", basename(path)))
    } else {
      NULL
    }
  }) %>% bind_rows()
  
  map_df <- prob_df %>%
    rename(Longitude = lon, Latitude = lat) %>%
    filter(!is.na(Longitude), !is.na(Latitude), !is.na(PA_prob))
  
  ggplot(map_df, aes(x = Longitude, y = Latitude, color = PA_prob)) +
    borders("world", xlim = c(94, 142), ylim = c(-11, 6), fill = "gray95", colour = "gray70") +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_gradient(name = "PA Probability", low = "yellow", high = "darkgreen") +
    coord_fixed(ratio = 1.3, xlim = c(94, 142), ylim = c(-11, 6)) +
    labs(title = "Predicted Seagrass Presence Probability (2017-2023)",
         x = "Longitude", y = "Latitude") +
    theme_minimal()
}

# -----------------------------------------------------------------------------
# plot_pa_threshold_accuracy()
# Line chart of classification accuracy across a range of probability
# thresholds, with the chosen operating threshold marked.
#
# Input:  <result_dir>/threshold_eval_PA.csv
#         expected columns: threshold, accuracy
# Output: ggplot line + point chart with a dashed vertical line at `threshold`
# -----------------------------------------------------------------------------
plot_pa_threshold_accuracy <- function(result_dir = "../result", threshold = 0.6) {
  threshold_df <- read_csv(file.path(result_dir, "threshold_eval_PA.csv"), show_col_types = FALSE)
  
  ggplot(threshold_df, aes(x = threshold, y = accuracy)) +
    geom_line(color = "grey", linewidth = 1) +
    geom_point(size = 2) +
    geom_text(aes(label = round(accuracy, 2)), vjust = -0.8, size = 3.5) +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(title = "Threshold vs Accuracy", x = "Threshold", y = "Accuracy") +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_rf_tuning_heatmap()
# Generic heatmap for visualizing a Random Forest hyperparameter grid
# search. Not stage-specific -- pass in whichever tuning results data frame
# you want to inspect.
#
# Input:  grid_results, a data frame already in memory (not read from file)
#         must contain the columns named by x_var, y_var, facet_var, fill_var
#         defaults assume columns: mtry, nodesize, ntree, rmse
# Output: ggplot tile heatmap, faceted by facet_var, cell values labelled
# -----------------------------------------------------------------------------
plot_rf_tuning_heatmap <- function(grid_results,
                                   x_var = "mtry",
                                   y_var = "nodesize",
                                   facet_var = "ntree",
                                   fill_var = "rmse",
                                   title = "RF Tuning Heatmap (RMSE)") {
  stopifnot(all(c(x_var, y_var, fill_var) %in% colnames(grid_results)))
  
  ggplot(grid_results, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[fill_var]])) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(.data[[fill_var]], 1)), color = "white", size = 3) +
    scale_fill_gradient(low = "steelblue", high = "red", name = toupper(fill_var)) +
    facet_wrap(vars(.data[[facet_var]])) +
    labs(title = title, x = x_var, y = y_var) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_pa_vs_tspc()
# Scatter plot of predicted PA probability against predicted percent cover
# (tSPC), coloured by the actual PA class. Used to sanity-check that high
# occurrence probability lines up with non-trivial cover.
#
# Input:  <result_dir>/predicted_PA_<year>.csv for year in 2017:2023
#         expected columns: PA, PA_prob, tSPC
# Output: ggplot scatter plot with reference lines at `threshold` (x) and
#         0.1 (y); emits a warning and returns NULL if tSPC is unavailable
# -----------------------------------------------------------------------------
plot_pa_vs_tspc <- function(result_dir = "../result", threshold = 0.6) {
  years <- 2017:2023
  files <- paste0(result_dir, "/predicted_PA_", years, ".csv")
  prob_df <- lapply(files, function(path) {
    if (file.exists(path)) {
      read_csv(path, show_col_types = FALSE) %>%
        filter(!is.na(PA), !is.na(PA_prob)) %>%
        mutate(year = gsub("\\D", "", basename(path)))
    } else {
      NULL
    }
  }) %>% bind_rows()
  
  df_plot <- prob_df %>% filter(!is.na(tSPC))
  
  if (nrow(df_plot) > 0) {
    ggplot(df_plot, aes(x = PA_prob, y = tSPC, color = factor(PA))) +
      geom_point(alpha = 0.5, size = 2) +
      scale_color_manual(values = c("blue", "darkgreen"), labels = c("PA = 0", "PA = 1")) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
      geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 0.1, linetype = "dotted", color = "red") +
      labs(title = "PA Probability vs tSPC",
           x = "PA Probability", y = "tSPC", color = "Actual PA") +
      theme_minimal()
  } else {
    warning("No data available for tSPC-based plot.")
  }
}

# -----------------------------------------------------------------------------
# plot_pa_vs_b2_stddev()
# Scatter plot of predicted PA probability against the B2 band standard
# deviation predictor, coloured by the actual PA class.
#
# Input:  <result_dir>/predicted_PA_<year>.csv for year in 2017:2023
#         expected columns: PA, PA_prob, B2_stddev
# Output: ggplot scatter plot with a dashed reference line at `threshold`;
#         emits a warning and returns NULL if B2_stddev is unavailable
# -----------------------------------------------------------------------------
plot_pa_vs_b2_stddev <- function(result_dir = "../result", threshold = 0.7) {
  years <- 2017:2023
  files <- paste0(result_dir, "/predicted_PA_", years, ".csv")
  prob_df <- lapply(files, function(path) {
    if (file.exists(path)) {
      read_csv(path, show_col_types = FALSE) %>%
        filter(!is.na(PA), !is.na(PA_prob)) %>%
        mutate(year = gsub("\\D", "", basename(path)))
    } else {
      NULL
    }
  }) %>% bind_rows()
  
  df_plot <- prob_df %>% filter(!is.na(B2_stddev))
  
  if (nrow(df_plot) > 0) {
    ggplot(df_plot, aes(x = PA_prob, y = B2_stddev, color = factor(PA))) +
      geom_point(alpha = 0.5, size = 2) +
      scale_color_manual(values = c("blue", "darkgreen"), labels = c("PA = 0", "PA = 1")) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
      geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
      labs(title = "PA Probability vs B2_stddev",
           x = "PA Probability", y = "B2_stddev", color = "Actual PA") +
      theme_minimal()
  } else {
    warning("No data available for B2_stddev-based plot.")
  }
}

# -----------------------------------------------------------------------------
# plot_pa_boxplot_yearly()
# Boxplot of predicted PA probability by actual class, for a single survey
# year. Also reports any readr parsing issues found in that year's file.
#
# Input:  <result_dir>/predicted_PA_<year>.csv
#         expected columns: PA, PA_prob
# Output: ggplot boxplot with a dashed reference line at 0.7; returns NULL
#         (with a warning) if the file is missing or only one class is present
# -----------------------------------------------------------------------------
plot_pa_boxplot_yearly <- function(year, result_dir = "../result") {
  file_path <- file.path(result_dir, paste0("predicted_PA_", year, ".csv"))
  
  if (!file.exists(file_path)) {
    warning(sprintf("File not found for year: %s", year))
    return(NULL)
  }
  
  df <- readr::read_csv(
    file_path,
    show_col_types = FALSE,
    guess_max = 100000,
    na = c("", "NA", "NaN", "null"),
    col_select = tidyselect::everything()
  )
  
  # Report any parsing issues found by readr
  pb <- readr::problems(df)
  if (nrow(pb) > 0) {
    message(sprintf("Parsing issues in year %s: %d rows", year, nrow(pb)))
    print(utils::head(pb, 5))
  }
  
  # Keep only rows with a valid PA and PA_prob
  df <- df %>%
    dplyr::filter(!is.na(PA), !is.na(PA_prob))
  
  # Only proceed if both classes (0 and 1) are present
  pa_classes <- unique(df$PA)
  if (length(pa_classes) < 2) {
    warning(sprintf("Year %s: only one class (%s) found -- boxplot skipped.", year, pa_classes))
    return(NULL)
  }
  
  df <- df %>%
    dplyr::mutate(PA = factor(PA, levels = c(0, 1), labels = c("Absent", "Present")))
  
  p <- ggplot(df, aes(x = PA, y = PA_prob, fill = PA)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    scale_fill_manual(values = c("Absent" = "#D95F02", "Present" = "#1B9E77")) +
    labs(
      title = paste("Predicted PA Probability by Class - Year", year),
      x = "Actual Class", y = "Predicted PA Probability"
    ) +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
  
  return(p)
}

# -----------------------------------------------------------------------------
# plot_pa_prob_distribution_optimal()
# Finds the classification threshold that maximises F1 score or Youden's J
# statistic for a single survey year, then plots the predicted-probability
# distribution with that optimal threshold marked.
#
# Input:  <result_dir>/predicted_PA_<year>.csv
#         expected columns: PA, PA_prob
# Params: method = "f1" (default) or "youden", selects which statistic to
#         maximise when searching thresholds 0.10-0.90 in steps of 0.01
# Output: ggplot histogram with a dashed vertical line at the optimal
#         threshold; returns NULL if the file is missing or empty
# -----------------------------------------------------------------------------
plot_pa_prob_distribution_optimal <- function(year, result_dir = "../result", method = "f1") {
  file_path <- file.path(result_dir, paste0("predicted_PA_", year, ".csv"))
  if (!file.exists(file_path)) {
    warning(sprintf("File not found for year: %s", year))
    return(NULL)
  }
  
  df <- readr::read_csv(
    file_path,
    show_col_types = FALSE,
    guess_max = 100000,
    na = c("", "NA", "NaN", "null"),
    col_select = c(PA, PA_prob)
  ) %>%
    dplyr::filter(!is.na(PA), !is.na(PA_prob)) %>%
    dplyr::mutate(PA = as.integer(as.character(PA)))
  
  if (nrow(df) == 0) return(NULL)
  
  # Search for the optimal threshold
  thresholds <- seq(0.1, 0.9, by = 0.01)
  metrics <- purrr::map_dfr(thresholds, function(th) {
    pred <- ifelse(df$PA_prob >= th, 1, 0)
    TP <- sum(pred == 1 & df$PA == 1)
    FP <- sum(pred == 1 & df$PA == 0)
    FN <- sum(pred == 0 & df$PA == 1)
    TN <- sum(pred == 0 & df$PA == 0)
    
    precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
    recall    <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
    f1        <- ifelse(!is.na(precision) && !is.na(recall) && (precision+recall)>0,
                        2 * precision * recall / (precision + recall), NA)
    sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
    spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
    youden <- ifelse(!is.na(sens) & !is.na(spec), sens + spec - 1, NA)
    
    tibble::tibble(threshold = th, f1 = f1, youden = youden)
  })
  
  if (method == "f1") {
    best_row <- metrics[which.max(metrics$f1), ]
  } else {
    best_row <- metrics[which.max(metrics$youden), ]
  }
  best_th <- best_row$threshold
  
  message(sprintf("Year %s: optimal threshold = %.2f (method: %s)", year, best_th, method))
  
  # Plot the probability distribution
  df$PA <- factor(df$PA, levels = c(0, 1), labels = c("Absent", "Present"))
  
  p <- ggplot(df, aes(x = PA_prob, fill = PA)) +
    geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6) +
    geom_vline(xintercept = best_th, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = best_th, y = 50, label = paste0("Optimal th = ", round(best_th, 2)),
             angle = 90, vjust = -0.5, color = "red") +
    scale_fill_manual(values = c("Absent" = "#999999", "Present" = "#1B9E77")) +
    labs(
      title = paste("Distribution of Predicted PA Probability - Year", year),
      x = "Predicted Probability", y = "Count", fill = "Actual Class"
    ) +
    theme_minimal(base_size = 14)
  
  return(p)
}


# =============================================================================
# SECTION 2 -- Morphology (MORPH3, 3-class) visualization functions
# =============================================================================

# -----------------------------------------------------------------------------
# plot_morph3_performance_by_year()
# Grouped bar chart of morphology classification accuracy, macro-F1, and
# kappa for each survey year.
#
# Input:  <result_dir>/morph3_performance_by_year.csv
#         expected columns: year, accuracy, macro_f1, kappa
# Output: ggplot bar chart, one group of bars per year
# -----------------------------------------------------------------------------
plot_morph3_performance_by_year <- function(result_dir = "../result") {
  
  perf_df <- read_csv(
    file.path(result_dir, "morph3_performance_by_year.csv"),
    show_col_types = FALSE
  )
  
  metrics_df <- perf_df %>%
    select(year, accuracy, macro_f1, kappa) %>%
    pivot_longer(cols = c(accuracy, macro_f1, kappa),
                 names_to = "metric",
                 values_to = "value") %>%
    filter(!is.na(value))
  
  ggplot(metrics_df, aes(x = year, y = value, fill = metric)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = round(value, 2)),
              position = position_dodge(0.9),
              vjust = -0.5,
              size = 3.5) +
    labs(
      title = "Morphology (morph3) Classification Performance per Year",
      x = "Year",
      y = "Metric Value"
    ) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_morph3_prob_distribution()
# Histogram of predicted class probability for each of the three morphology
# classes, faceted by predicted class and coloured by the true class.
# Predictions and ground truth are joined on gee_id.
#
# Input:  <result_dir>/predicted_MORPH3_probs_<year>.csv (columns: gee_id,
#         P_mixed_short_plus_mono_short, P_mixed_long, P_mono_Ea, ...)
#         <result_dir>/dataCARBON_<year>.csv (columns: gee_id, morph3)
#         for year in 2017:2024
# Output: ggplot histogram faceted by predicted class, with manuscript-
#         friendly class labels (short-leaf / mixed-leaf / long-leaf)
# -----------------------------------------------------------------------------
plot_morph3_prob_distribution <- function(result_dir = "../result") {
  
  years <- 2017:2024
  
  df_all <- purrr::map_dfr(years, function(yr) {
    
    file_path <- file.path(result_dir,
                           paste0("predicted_MORPH3_probs_", yr, ".csv"))
    
    truth_path <- file.path(result_dir,
                            paste0("dataCARBON_", yr, ".csv"))
    
    if (!file.exists(file_path) || !file.exists(truth_path)) return(NULL)
    
    df_pred <- read_csv(file_path, show_col_types = FALSE)
    df_truth <- read_csv(truth_path, show_col_types = FALSE) %>%
      select(gee_id, morph3)
    
    inner_join(df_truth, df_pred, by = "gee_id") %>%
      filter(!is.na(morph3)) %>%
      mutate(year = yr)
  })
  
  prob_cols <- grep("^P_", names(df_all), value = TRUE)
  
  df_long <- df_all %>%
    pivot_longer(cols = all_of(prob_cols),
                 names_to = "class",
                 values_to = "probability")
  
  # Relabel classes for manuscript consistency
  df_long <- df_long %>%
    mutate(
      morph3 = recode(morph3,
                      "mixed_short_plus_mono_short" = "short-leaf",
                      "mixed_long"                  = "mixed-leaf",
                      "mono_Ea"                     = "long-leaf"),
      class = recode(class,
                     "P_mixed_short_plus_mono_short" = "short-leaf",
                     "P_mixed_long"                  = "mixed-leaf",
                     "P_mono_Ea"                     = "long-leaf")
    )
  
  ggplot(df_long, aes(x = probability, fill = morph3)) +
    geom_histogram(binwidth = 0.05, alpha = 0.6) +
    facet_wrap(~ class, scales = "free_y") +
    labs(
      x = "Predicted Probability",
      y = "Count"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# plot_morph3_spatial_map()
# Point map of one morphology class's predicted probability across
# Indonesia, pooled across all available years.
#
# Input:  <result_dir>/predicted_MORPH3_probs_<year>.csv for year in
#         2017:2024; expected columns: lon, lat, and the column named by
#         target_class (default "P_mono_Ea")
# Output: ggplot map (yellow = low probability, dark red = high probability);
#         stops with an error if target_class is not present in the data
# -----------------------------------------------------------------------------
plot_morph3_spatial_map <- function(result_dir = "../result",
                                    target_class = "P_mono_Ea") {
  
  years <- 2017:2024
  
  df_all <- purrr::map_dfr(years, function(yr) {
    
    file_path <- file.path(result_dir,
                           paste0("predicted_MORPH3_probs_", yr, ".csv"))
    
    if (!file.exists(file_path)) return(NULL)
    
    read_csv(file_path, show_col_types = FALSE) %>%
      filter(!is.na(lat), !is.na(lon)) %>%
      mutate(year = yr)
  })
  
  if (!target_class %in% names(df_all)) {
    stop("Requested class probability column not found.")
  }
  
  ggplot(df_all, aes(x = lon, y = lat, color = .data[[target_class]])) +
    geom_point(size = 2.2, alpha = 0.7) +
    coord_fixed() +
    scale_color_gradient(low = "yellow", high = "darkred") +
    labs(
      title = paste("Spatial Distribution of", target_class),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal()
}

# -----------------------------------------------------------------------------
# plot_morph3_boxplot_yearly()
# Boxplot of predicted class probability by the true morphology class, for
# a single survey year, faceted by predicted class. Predictions and ground
# truth are joined on gee_id.
#
# Input:  <result_dir>/predicted_MORPH3_probs_<year>.csv
#         <result_dir>/dataCARBON_<year>.csv (column: gee_id, morph3)
# Output: ggplot boxplot faceted by predicted class, with manuscript-
#         friendly class labels; returns NULL if either file is missing
# -----------------------------------------------------------------------------
plot_morph3_boxplot_yearly <- function(year, result_dir = "../result") {
  
  file_path <- file.path(result_dir,
                         paste0("predicted_MORPH3_probs_", year, ".csv"))
  
  truth_path <- file.path(result_dir,
                          paste0("dataCARBON_", year, ".csv"))
  
  if (!file.exists(file_path) || !file.exists(truth_path)) return(NULL)
  
  df_pred <- read_csv(file_path, show_col_types = FALSE)
  df_truth <- read_csv(truth_path, show_col_types = FALSE) %>%
    select(gee_id, morph3)
  
  df <- inner_join(df_truth, df_pred, by = "gee_id") %>%
    filter(!is.na(morph3))
  
  prob_cols <- grep("^P_", names(df), value = TRUE)
  
  df_long <- df %>%
    pivot_longer(cols = all_of(prob_cols),
                 names_to = "class",
                 values_to = "probability")
  
  # Relabel classes for manuscript consistency
  df_long <- df_long %>%
    mutate(
      morph3 = recode(morph3,
                      "mixed_short_plus_mono_short" = "short-leaf",
                      "mixed_long"                  = "mixed-leaf",
                      "mono_Ea"                     = "long-leaf"),
      class = recode(class,
                     "P_mixed_short_plus_mono_short" = "short-leaf",
                     "P_mixed_long"                  = "mixed-leaf",
                     "P_mono_Ea"                     = "long-leaf")
    )
  
  ggplot(df_long, aes(x = morph3, y = probability, fill = morph3)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
    facet_wrap(~ class) +
    labs(
      title = paste("Morphology Probability Boxplot - Year", year),
      x = "True Morphology Class",
      y = "Predicted Probability"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
}


# =============================================================================
# SECTION 3 -- SPC / percent cover (PCT) visualization functions
# =============================================================================

# -----------------------------------------------------------------------------
# plot_regression_performance_by_year()
# Grouped bar chart of SPC (percent cover) RMSE and MAE for each survey year.
#
# Input:  <result_dir>/pct_performance_by_year.csv
#         expected columns: year, RMSE, MAE
# Output: ggplot bar chart, one group of bars per year
# -----------------------------------------------------------------------------
plot_regression_performance_by_year <- function(result_dir = "../result") {
  df <- read_csv(file.path(result_dir, "pct_performance_by_year.csv"), show_col_types = FALSE)
  
  perf_long <- df %>%
    mutate(year = as.character(year)) %>%
    pivot_longer(cols = c(RMSE, MAE), names_to = "metric", values_to = "value") %>%
    mutate(year = factor(year, levels = unique(year)))
  
  ggplot(perf_long, aes(x = year, y = value, fill = metric)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = round(value, 2)), position = position_dodge(0.9), vjust = -0.5) +
    labs(title = "PCT Model Performance per Year", y = "Metric Value", x = "Year") +
    scale_fill_manual(values = c("RMSE" = "#555555", "MAE" = "#BBBBBB")) +
    theme_minimal()
}

# -----------------------------------------------------------------------------
# plot_regression_performance_by_location()
# Horizontal bar chart of SPC RMSE, MAE, and R2 for each survey location,
# sorted by metric value. Non-finite values (e.g. Inf from a degenerate
# location) are dropped before plotting.
#
# Input:  <result_dir>/pct_performance_by_location.csv
#         expected columns: loc, RMSE, MAE, R2
# Output: ggplot horizontal bar chart
# -----------------------------------------------------------------------------
plot_regression_performance_by_location <- function(result_dir = "../result") {
  df <- read_csv(file.path(result_dir, "pct_performance_by_location.csv"), show_col_types = FALSE)
  
  perf_long <- df %>%
    pivot_longer(cols = c(RMSE, MAE, R2), names_to = "metric", values_to = "value") %>%
    filter(!is.na(value), is.finite(value))  # drop NA and Inf values
  
  ggplot(perf_long, aes(x = reorder(loc, value), y = value, fill = metric)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = round(value, 2)), position = position_dodge(0.9), vjust = -0.5, size = 3.5) +
    coord_flip() +
    labs(title = "PCT Model Performance by Location", x = "Location", y = "Metric Value") +
    scale_fill_manual(values = c("RMSE" = "#888888", "MAE" = "#BBBBBB", "R2" = "#4CAF50")) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_regression_boxplot_actual_vs_pred_bin()
# Boxplot of actual tSPC values, grouped into bins of predicted tSPC value,
# pooled across all available years. Values outside (0, 100] are dropped.
#
# Input:  <result_dir>/predicted_tSPC_<year>.csv for year in 2017:2023
#         expected columns: tSPC, tSPC_pred
# Output: ggplot boxplot, one box per 10-unit predicted-value bin
# -----------------------------------------------------------------------------
plot_regression_boxplot_actual_vs_pred_bin <- function(result_dir = "../result") {
  years <- 2017:2023
  files <- paste0(result_dir, "/predicted_tSPC_", years, ".csv")
  
  scatter_df <- lapply(files, function(path) {
    if (file.exists(path)) {
      read_csv(path, show_col_types = FALSE) %>%
        filter(
          !is.na(tSPC), !is.na(tSPC_pred),
          tSPC > 0, tSPC <= 100,
          tSPC_pred > 0, tSPC_pred <= 100
        ) %>%
        mutate(year = gsub("\\D", "", basename(path)))
    } else {
      NULL
    }
  }) %>% bind_rows()
  
  scatter_df <- scatter_df %>%
    mutate(Pred_Bin = cut(tSPC_pred, breaks = seq(0, 100, by = 10), include.lowest = TRUE))
  
  ggplot(scatter_df, aes(x = Pred_Bin, y = tSPC)) +
    geom_boxplot(fill = "gray85") +
    labs(
      title = "Actual tSPC by Predicted Interval",
      x = "Predicted tSPC Bin (%)",
      y = "Actual tSPC (%)"
    ) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_regression_error_by_actual_bin()
# Boxplot of absolute prediction error, grouped into bins of actual tSPC
# value, pooled across all available years. Values outside (0, 100] are
# dropped.
#
# Input:  <result_dir>/predicted_tSPC_<year>.csv for year in 2017:2023
#         expected columns: tSPC, tSPC_pred
# Output: ggplot boxplot, one box per 10-unit actual-value bin
# -----------------------------------------------------------------------------
plot_regression_error_by_actual_bin <- function(result_dir = "../result") {
  years <- 2017:2023
  files <- paste0(result_dir, "/predicted_tSPC_", years, ".csv")
  
  error_df <- lapply(files, function(path) {
    if (file.exists(path)) {
      read_csv(path, show_col_types = FALSE) %>%
        filter(
          !is.na(tSPC), !is.na(tSPC_pred),
          tSPC > 0, tSPC <= 100,
          tSPC_pred > 0, tSPC_pred <= 100
        ) %>%
        mutate(year = gsub("\\D", "", basename(path)))
    } else {
      NULL
    }
  }) %>% bind_rows()
  
  error_df <- error_df %>%
    mutate(
      Error = abs(tSPC - tSPC_pred),
      Actual_Bin = cut(tSPC, breaks = seq(0, 100, by = 10), include.lowest = TRUE)
    )
  
  ggplot(error_df, aes(x = Actual_Bin, y = Error)) +
    geom_boxplot(fill = "gray85") +
    labs(
      title = "Absolute Error by Actual tSPC Interval",
      x = "Actual tSPC Bin (%)",
      y = "Absolute Error"
    ) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_regression_spatial_map()
# Point map of predicted tSPC across Indonesia, pooled across all available
# years.
#
# Input:  <result_dir>/predicted_tSPC_<year>.csv for year in 2017:2024
#         expected columns: tSPC_pred, lon, lat
# Output: ggplot map (yellow = low cover, dark green = high cover); returns
#         NULL (with a warning) if no rows have valid coordinates
# -----------------------------------------------------------------------------
plot_regression_spatial_map <- function(result_dir = "../result") {
  years <- 2017:2024
  files <- paste0(result_dir, "/predicted_tSPC_", years, ".csv")
  
  map_df <- lapply(files, function(path) {
    if (file.exists(path)) {
      read_csv(path, show_col_types = FALSE) %>%
        filter(!is.na(tSPC_pred), !is.na(lon), !is.na(lat)) %>%
        rename(Longitude = lon, Latitude = lat)
    } else {
      NULL
    }
  }) %>% bind_rows()
  
  if (nrow(map_df) == 0) {
    warning("No data with valid coordinates found for spatial plot.")
    return(NULL)
  }
  
  ggplot(map_df, aes(x = Longitude, y = Latitude, color = tSPC_pred)) +
    borders("world", xlim = c(94, 142), ylim = c(-11, 6), fill = "gray95", colour = "gray70") +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_gradient(name = "tSPC Prediction", low = "yellow", high = "darkgreen") +
    coord_fixed(ratio = 1.3, xlim = c(94, 142), ylim = c(-11, 6)) +
    labs(title = "Predicted tSPC Map (2017-2024)", x = "Longitude", y = "Latitude") +
    theme_minimal()
}


# =============================================================================
# SECTION 4 -- AGB (above-ground biomass) visualization functions
#
# Unlike the PA / MORPH3 / SPC functions above, these take an already-loaded
# data frame (df) rather than reading CSVs from result_dir directly.
# =============================================================================

# -----------------------------------------------------------------------------
# plot_regression_performance_by_year_agb()
# Grouped bar chart of AGB RMSE and MAE for each survey year, computed
# directly from df.
#
# Input:  df with columns AGB_pred (reference/proxy AGB), tAGB_pred
#         (model-predicted AGB), year
# Output: ggplot bar chart, one group of bars per year
# -----------------------------------------------------------------------------
plot_regression_performance_by_year_agb <- function(df) {
  library(ggplot2)
  library(dplyr)
  library(Metrics)
  library(tidyr)
  
  # Compute metrics per year
  summary_df <- df %>%
    filter(!is.na(AGB_pred), !is.na(tAGB_pred), !is.na(year)) %>%
    group_by(year) %>%
    summarise(
      RMSE = Metrics::rmse(AGB_pred, tAGB_pred),
      MAE  = Metrics::mae(AGB_pred, tAGB_pred),
      R2   = caret::R2(tAGB_pred, AGB_pred),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(MAE, RMSE), names_to = "metric", values_to = "value")
  
  # Plot
  ggplot(summary_df, aes(x = factor(year), y = value, fill = metric)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = round(value, 2)), position = position_dodge(0.9),
              vjust = -0.5, size = 3.5) +
    labs(
      title = "AGB Model Performance per Year",
      x = "Year",
      y = "Metric Value"
    ) +
    scale_fill_manual(values = c("MAE" = "gray60", "RMSE" = "gray30")) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_regression_performance_by_location_agb()
# Horizontal bar chart of AGB RMSE, MAE, and R2 for each survey location,
# sorted by metric value.
#
# Input:  df with columns loc, AGB_pred, tAGB_pred
# Output: ggplot horizontal bar chart
# -----------------------------------------------------------------------------
plot_regression_performance_by_location_agb <- function(df) {
  perf_long <- df %>%
    filter(!is.na(loc), !is.na(AGB_pred), !is.na(tAGB_pred)) %>%
    group_by(loc) %>%
    summarise(
      RMSE = Metrics::rmse(AGB_pred, tAGB_pred),
      MAE  = Metrics::mae(AGB_pred, tAGB_pred),
      R2   = caret::R2(tAGB_pred, AGB_pred),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(RMSE, MAE, R2), names_to = "metric", values_to = "value") %>%
    filter(!is.na(value), is.finite(value))
  
  ggplot(perf_long, aes(x = reorder(loc, value), y = value, fill = metric)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = round(value, 2)), position = position_dodge(0.9), vjust = -0.5, size = 3.5) +
    coord_flip() +
    labs(title = "AGB Model Performance by Location", x = "Location", y = "Metric Value") +
    scale_fill_manual(values = c("RMSE" = "#59A14F", "MAE" = "#EDC948", "R2" = "#B07AA1")) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_regression_boxplot_actual_vs_pred_bin_agb()
# Boxplot of actual AGB values (AGB_pred), grouped into bins of predicted
# AGB value (tAGB_pred).
#
# Input:  df with columns AGB_pred, tAGB_pred
# Output: ggplot boxplot, one box per predicted-value bin (bin edges chosen
#         automatically via pretty())
# -----------------------------------------------------------------------------
plot_regression_boxplot_actual_vs_pred_bin_agb <- function(df) {
  df <- df %>%
    filter(!is.na(AGB_pred), !is.na(tAGB_pred)) %>%
    mutate(Pred_Bin = cut(
      tAGB_pred,
      breaks = pretty(tAGB_pred, n = 10),
      include.lowest = TRUE
    ))
  
  ggplot(df, aes(x = Pred_Bin, y = AGB_pred)) +
    geom_boxplot(fill = "gray85", color = "black") +
    labs(
      title = "Actual AGB vs Predicted AGB Bins",
      x = "Predicted AGB (tAGB_pred) Bins",
      y = "Actual AGB (AGB_pred)"
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -----------------------------------------------------------------------------
# plot_regression_spatial_map_agb()
# Point map of predicted AGB across Indonesia.
#
# Input:  df with columns tAGB_pred, lon, lat
# Output: ggplot map (cyan = low biomass, dark blue = high biomass); returns
#         NULL (with a warning) if no rows have valid coordinates
# -----------------------------------------------------------------------------
plot_regression_spatial_map_agb <- function(df) {
  df <- df %>%
    filter(!is.na(tAGB_pred), !is.na(lon), !is.na(lat)) %>%
    rename(Longitude = lon, Latitude = lat)
  
  if (nrow(df) == 0) {
    warning("No valid spatial data for plotting.")
    return(NULL)
  }
  
  ggplot(df, aes(x = Longitude, y = Latitude, color = tAGB_pred)) +
    borders("world", xlim = c(94, 142), ylim = c(-11, 6), fill = "gray95", colour = "gray70") +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_gradient(name = "AGB Prediction", low = "cyan", high = "darkblue") +
    coord_fixed(ratio = 1.3, xlim = c(94, 142), ylim = c(-11, 6)) +
    labs(title = "Predicted AGB Map (2017-2023)", x = "Longitude", y = "Latitude") +
    theme_minimal()
}

# -----------------------------------------------------------------------------
# plot_regression_error_by_actual_bin_agb()
# Boxplot of absolute prediction error, grouped into bins of actual AGB
# value (AGB_pred).
#
# Input:  df with columns AGB_pred, tAGB_pred
# Output: ggplot boxplot, one box per actual-value bin (bin edges chosen
#         automatically via pretty())
# -----------------------------------------------------------------------------
plot_regression_error_by_actual_bin_agb <- function(df) {
  df <- df %>%
    filter(!is.na(AGB_pred), !is.na(tAGB_pred)) %>%
    mutate(
      Error = abs(AGB_pred - tAGB_pred),
      Actual_Bin = cut(AGB_pred, breaks = pretty(AGB_pred, n = 10), include.lowest = TRUE)
    )
  
  ggplot(df, aes(x = Actual_Bin, y = Error)) +
    geom_boxplot(fill = "gray85") +
    labs(
      title = "Absolute Error by Actual AGB Bin",
      x = "Actual AGB Bin",
      y = "Absolute Error"
    ) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_agb_rmse_with_range_by_location()
# Horizontal bar chart of AGB RMSE per location, annotated with each
# location's predicted AGB range (min-max) as a text label.
#
# Input:  df with columns loc, AGB_pred, tAGB_pred
# Output: ggplot horizontal bar chart, sorted by RMSE
# -----------------------------------------------------------------------------
plot_agb_rmse_with_range_by_location <- function(df) {
  # Compute RMSE per location
  perf <- df %>%
    filter(!is.na(loc), !is.na(AGB_pred), !is.na(tAGB_pred)) %>%
    group_by(loc) %>%
    summarise(
      RMSE = Metrics::rmse(AGB_pred, tAGB_pred),
      AGB_min = min(tAGB_pred, na.rm = TRUE),
      AGB_max = max(tAGB_pred, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(perf, aes(x = reorder(loc, RMSE), y = RMSE)) +
    geom_col(fill = "gray85") +
    geom_text(aes(label = round(RMSE, 1)), hjust = -0.1, size = 3.5) +
    geom_text(aes(y = 0, label = paste0("(", round(AGB_min), "-", round(AGB_max), ")")),
              hjust = -0.1, size = 3.2, color = "gray30") +
    coord_flip() +
    labs(
      title = "AGB RMSE by Location with Predicted AGB Range",
      x = "Location",
      y = "RMSE (tAGB_pred)",
      subtitle = "Text label shows predicted AGB range (min-max)"
    ) +
    theme_minimal(base_size = 14)
}


# =============================================================================
# SECTION 5 -- AGC (above-ground carbon, final stage) visualization functions
#
# Structurally identical to the AGB functions in Section 4 above; only the
# column names differ (AGC_pred / tAGC_pred instead of AGB_pred / tAGB_pred).
# =============================================================================

# -----------------------------------------------------------------------------
# plot_regression_performance_by_year_agc()
# Grouped bar chart of AGC RMSE and MAE for each survey year, computed
# directly from df.
#
# Input:  df with columns AGC_pred (reference/proxy AGC), tAGC_pred
#         (model-predicted AGC), year
# Output: ggplot bar chart, one group of bars per year
# -----------------------------------------------------------------------------
plot_regression_performance_by_year_agc <- function(df) {
  summary_df <- df %>%
    filter(!is.na(AGC_pred), !is.na(tAGC_pred), !is.na(year)) %>%
    group_by(year) %>%
    summarise(
      MAE = mae(AGC_pred, tAGC_pred),
      RMSE = rmse(AGC_pred, tAGC_pred),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(MAE, RMSE), names_to = "metric", values_to = "value")
  
  ggplot(summary_df, aes(x = factor(year), y = value, fill = metric)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = round(value, 2)), position = position_dodge(0.9),
              vjust = -0.5, size = 3.5) +
    labs(
      title = "AGC Model Performance per Year",
      x = "Year",
      y = "Metric Value"
    ) +
    scale_fill_manual(values = c("MAE" = "gray60", "RMSE" = "gray30")) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_regression_performance_by_location_agc()
# Horizontal bar chart of AGC RMSE, MAE, and R2 for each survey location,
# sorted by metric value.
#
# Input:  df with columns loc, AGC_pred, tAGC_pred
# Output: ggplot horizontal bar chart
# -----------------------------------------------------------------------------
plot_regression_performance_by_location_agc <- function(df) {
  perf_long <- df %>%
    filter(!is.na(loc), !is.na(AGC_pred), !is.na(tAGC_pred)) %>%
    group_by(loc) %>%
    summarise(
      RMSE = rmse(AGC_pred, tAGC_pred),
      MAE  = mae(AGC_pred, tAGC_pred),
      R2   = R2(tAGC_pred, AGC_pred),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(RMSE, MAE, R2), names_to = "metric", values_to = "value") %>%
    filter(!is.na(value), is.finite(value))
  
  ggplot(perf_long, aes(x = reorder(loc, value), y = value, fill = metric)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = round(value, 2)), position = position_dodge(0.9),
              vjust = -0.5, size = 3.5) +
    coord_flip() +
    labs(title = "AGC Model Performance by Location", x = "Location", y = "Metric Value") +
    scale_fill_manual(values = c("RMSE" = "#59A14F", "MAE" = "#EDC948", "R2" = "#B07AA1")) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_regression_boxplot_actual_vs_pred_bin_agc()
# Boxplot of actual AGC values (AGC_pred), grouped into bins of predicted
# AGC value (tAGC_pred).
#
# Input:  df with columns AGC_pred, tAGC_pred
# Output: ggplot boxplot, one box per predicted-value bin (bin edges chosen
#         automatically via pretty())
# -----------------------------------------------------------------------------
plot_regression_boxplot_actual_vs_pred_bin_agc <- function(df) {
  df <- df %>%
    filter(!is.na(AGC_pred), !is.na(tAGC_pred)) %>%
    mutate(Pred_Bin = cut(
      tAGC_pred,
      breaks = pretty(tAGC_pred, n = 10),
      include.lowest = TRUE
    ))
  
  ggplot(df, aes(x = Pred_Bin, y = AGC_pred)) +
    geom_boxplot(fill = "gray85", color = "black") +
    labs(
      title = "Actual AGC vs Predicted AGC Bins",
      x = "Predicted AGC (tAGC_pred) Bins",
      y = "Actual AGC (AGC_pred)"
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -----------------------------------------------------------------------------
# plot_regression_spatial_map_agc()
# Point map of predicted AGC across Indonesia.
#
# Input:  df with columns tAGC_pred, lon, lat
# Output: ggplot map (cyan = low carbon, dark blue = high carbon); returns
#         NULL (with a warning) if no rows have valid coordinates
# -----------------------------------------------------------------------------
plot_regression_spatial_map_agc <- function(df) {
  df <- df %>%
    filter(!is.na(tAGC_pred), !is.na(lon), !is.na(lat)) %>%
    rename(Longitude = lon, Latitude = lat)
  
  if (nrow(df) == 0) {
    warning("No valid spatial data for plotting.")
    return(NULL)
  }
  
  ggplot(df, aes(x = Longitude, y = Latitude, color = tAGC_pred)) +
    borders("world", xlim = c(94, 142), ylim = c(-11, 6),
            fill = "gray95", colour = "gray70") +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_gradient(name = "AGC Prediction", low = "cyan", high = "darkblue") +
    coord_fixed(ratio = 1.3, xlim = c(94, 142), ylim = c(-11, 6)) +
    labs(title = "Predicted AGC Map (2017-2023)", x = "Longitude", y = "Latitude") +
    theme_minimal()
}

# -----------------------------------------------------------------------------
# plot_regression_error_by_actual_bin_agc()
# Boxplot of absolute prediction error, grouped into bins of actual AGC
# value (AGC_pred).
#
# Input:  df with columns AGC_pred, tAGC_pred
# Output: ggplot boxplot, one box per actual-value bin (bin edges chosen
#         automatically via pretty())
# -----------------------------------------------------------------------------
plot_regression_error_by_actual_bin_agc <- function(df) {
  df <- df %>%
    filter(!is.na(AGC_pred), !is.na(tAGC_pred)) %>%
    mutate(
      Error = abs(AGC_pred - tAGC_pred),
      Actual_Bin = cut(AGC_pred, breaks = pretty(AGC_pred, n = 10), include.lowest = TRUE)
    )
  
  ggplot(df, aes(x = Actual_Bin, y = Error)) +
    geom_boxplot(fill = "gray85") +
    labs(
      title = "Absolute Error by Actual AGC Bin",
      x = "Actual AGC Bin",
      y = "Absolute Error"
    ) +
    theme_minimal(base_size = 14)
}

# -----------------------------------------------------------------------------
# plot_agc_rmse_with_range_by_location()
# Horizontal bar chart of AGC RMSE per location, annotated with each
# location's predicted AGC range (min-max) as a text label.
#
# Input:  df with columns loc, AGC_pred, tAGC_pred
# Output: ggplot horizontal bar chart, sorted by RMSE
# -----------------------------------------------------------------------------
plot_agc_rmse_with_range_by_location <- function(df) {
  perf <- df %>%
    filter(!is.na(loc), !is.na(AGC_pred), !is.na(tAGC_pred)) %>%
    group_by(loc) %>%
    summarise(
      RMSE = rmse(AGC_pred, tAGC_pred),
      AGC_min = min(tAGC_pred, na.rm = TRUE),
      AGC_max = max(tAGC_pred, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(perf, aes(x = reorder(loc, RMSE), y = RMSE)) +
    geom_col(fill = "gray85") +
    geom_text(aes(label = round(RMSE, 1)), hjust = -0.1, size = 3.5) +
    geom_text(aes(y = 0, label = paste0("(", round(AGC_min), "-", round(AGC_max), ")")),
              hjust = -0.1, size = 3.2, color = "gray30") +
    coord_flip() +
    labs(
      title = "AGC RMSE by Location with Predicted AGC Range",
      x = "Location",
      y = "RMSE (tAGC_pred)",
      subtitle = "Text label shows predicted AGC range (min-max)"
    ) +
    theme_minimal(base_size = 14)
}