# =============================================================================
# plot_functions.R
# Plots used by the presence model scripts (02 and 04 in R/workflow/).
#
# Default folders point to result/ in the repository (here::here("result")).
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)


# =============================================================================
# PA (presence) plots
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
plot_pa_performance_by_year <- function(result_dir = here::here("result")) {
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
plot_pa_prob_distribution <- function(result_dir = here::here("result")) {
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
plot_pa_spatial_map <- function(result_dir = here::here("result")) {
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
plot_pa_threshold_accuracy <- function(result_dir = here::here("result"), threshold = 0.6) {
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
plot_pa_vs_tspc <- function(result_dir = here::here("result"), threshold = 0.6) {
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
# plot_pa_boxplot_yearly()
# Boxplot of predicted PA probability by actual class, for a single survey
# year. Also reports any readr parsing issues found in that year's file.
#
# Input:  <result_dir>/predicted_PA_<year>.csv
#         expected columns: PA, PA_prob
# Output: ggplot boxplot with a dashed reference line at 0.7; returns NULL
#         (with a warning) if the file is missing or only one class is present
# -----------------------------------------------------------------------------
plot_pa_boxplot_yearly <- function(year, result_dir = here::here("result")) {
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
plot_pa_prob_distribution_optimal <- function(year, result_dir = here::here("result"), method = "f1") {
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


