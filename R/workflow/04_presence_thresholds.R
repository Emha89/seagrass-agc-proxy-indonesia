# =============================================================================
# 04_presence_thresholds.R
# Accuracy of the presence model in each year, and the probability threshold
# with the highest F1-score in each year (values shown in Figure S1 of the
# paper). Each year compares the field labels with the probabilities predicted
# from the GSE layer of that year.
#
# Inputs   result/dataPA_<year>.csv, result/predicted_PA_<year>.csv
# Outputs  result/pa_performance_by_year.csv and plots
# =============================================================================

library(dplyr)
library(readr)
library(purrr)
library(caret)
library(ggplot2)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "functions", "plot_functions.R"))  # visualization functions

data_dir <- here("result")
data_years <- as.character(2017:2024)
out_csv <- file.path(data_dir, "pa_performance_by_year.csv")

# -----------------------------------------------------------------------------
# Helper Function: Safe division
# -----------------------------------------------------------------------------
safe_div <- function(x, y) ifelse(y == 0, NA, x / y)

evaluate_metrics <- function(actual, predicted) {
  TP <- sum(predicted == 1 & actual == 1, na.rm = TRUE)
  TN <- sum(predicted == 0 & actual == 0, na.rm = TRUE)
  FP <- sum(predicted == 1 & actual == 0, na.rm = TRUE)
  FN <- sum(predicted == 0 & actual == 1, na.rm = TRUE)
  
  acc <- safe_div(TP + TN, TP + TN + FP + FN)
  precision <- safe_div(TP, TP + FP)
  recall <- safe_div(TP, TP + FN)
  f1 <- ifelse(!is.na(precision) && !is.na(recall) && (precision + recall) > 0,
               2 * precision * recall / (precision + recall), NA)
  
  data.frame(
    accuracy = round(acc, 3),
    precision = round(precision, 3),
    recall = round(recall, 3),
    f1_score = round(f1, 3),
    TP = TP, FP = FP, FN = FN, TN = TN
  )
}

# -----------------------------------------------------------------------------
# Block 1: Evaluate Per-Year Metrics
# -----------------------------------------------------------------------------
eval_results <- map_dfr(data_years, function(yr) {
  pred_path <- file.path(data_dir, paste0("predicted_PA_", yr, ".csv"))
  truth_path <- file.path(data_dir, paste0("dataPA_", yr, ".csv"))
  
  if (!file.exists(pred_path) || !file.exists(truth_path)) {
    message(sprintf("Skipping year %s (missing file)", yr))
    return(NULL)
  }
  
  df_pred <- read_csv(pred_path, show_col_types = FALSE) %>%
    select(gee_id, PA_pred, PA_prob)
  
  df_truth <- read_csv(truth_path, show_col_types = FALSE) %>%
    select(gee_id, PA)
  
  df_eval <- inner_join(df_truth, df_pred, by = "gee_id") %>%
    filter(!is.na(PA), !is.na(PA_pred))
  
  if (nrow(df_eval) == 0) return(NULL)
  
  actual <- as.integer(as.character(df_eval$PA))
  predicted <- as.integer(df_eval$PA_pred)
  
  met <- evaluate_metrics(actual, predicted)
  met$year <- yr
  met$n <- nrow(df_eval)
  met
})

cat("Columns in evaluation results:\n")
print(names(eval_results))

write_csv(eval_results, out_csv)
cat(sprintf("Saved evaluation results to: %s\n", out_csv))
print(eval_results)

# -----------------------------------------------------------------------------
# Block 2: Visualization
# -----------------------------------------------------------------------------
# Accuracy, precision, recall, f1 per year
cat("\nPlotting PA performance by year...\n")
print(plot_pa_performance_by_year())

# Distribution of predicted probabilities by class
cat("\nPlotting PA probability distribution...\n")
print(plot_pa_prob_distribution())

# Spatial prediction map
cat("\nPlotting spatial map of PA predictions...\n")
print(plot_pa_spatial_map())

# Boxplot per year
for (yr in data_years) {
  cat(sprintf("\nPlotting boxplot for year %s...\n", yr))
  p <- plot_pa_boxplot_yearly(yr)
  if (!is.null(p)) print(p)
}

# Distribution per year with optimal threshold (F1-based)
for (yr in data_years) {
  cat(sprintf("\nPlotting optimal threshold distribution for year %s...\n", yr))
  p <- plot_pa_prob_distribution_optimal(yr, method = "f1")
  if (!is.null(p)) print(p)
}

# Additional: Plot PA vs tSPC or B2_stddev if available
cat("\nPlotting PA vs tSPC...\n")
print(plot_pa_vs_tspc())

cat("\nAll visualizations generated.\n")