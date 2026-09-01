# =============================================================================
# rf_evalModel_MORPH3.R
#
# Evaluate morphology (morph3) classification predictions per year.
# Output files are designed for supplementary tables and figures in the
# paper.
#
# Pipeline position: stage 2 of 6, final step. Follows
# rf_applyModel_MORPH3.R.
#
# FIX APPLIED: the original version compared morph3_pred (derived from
# P_-prefixed probability column names, e.g. "P_mono_Ea") directly against
# morph3 (the unprefixed ground-truth label, e.g. "mono_Ea"). Since these
# strings never match, factor(predicted, levels = levels(factor(actual)))
# turned every predicted value into NA, which caret::confusionMatrix()
# does not accept. Fixed by stripping the "P_" prefix before comparison --
# a one-line change that does not affect how any metric is computed, only
# how the predicted class label is spelled.
#
# Inputs:
#   - predicted_MORPH3_probs_<year>.csv  (from rf_applyModel_MORPH3.R)
#   - dataCARBON_<year>.csv              (ground-truth morph3 label)
#
# Outputs (written to result/):
#   1) morph3_performance_by_year.csv     accuracy, macro-F1, kappa per year
#   2) morph3_confusion_by_year.csv       full confusion matrices per year
#   3) morph3_classwise_metrics_by_year.csv  precision/recall/F1 per class per year
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

library(dplyr)
library(readr)
library(purrr)
library(caret)
library(ggplot2)
library(tibble)
library(here)   # install.packages("here") if you don't have it yet

source(here("R", "00_shared_functions", "rf_func_vis.R"))

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
data_dir <- here("result")
data_years <- as.character(2017:2024)

# Output file paths (paper supplementary)
out_perf <- file.path(data_dir, "morph3_performance_by_year.csv")
out_cm   <- file.path(data_dir, "morph3_confusion_by_year.csv")
out_cls  <- file.path(data_dir, "morph3_classwise_metrics_by_year.csv")

# -----------------------------------------------------------------------------
# Helper Function: Multi-class Evaluation Metrics
# -----------------------------------------------------------------------------
evaluate_metrics_multiclass <- function(actual, predicted) {
  
  actual <- factor(actual)
  predicted <- factor(predicted, levels = levels(actual))
  
  cm <- caret::confusionMatrix(predicted, actual)
  
  # Overall metrics
  acc <- as.numeric(cm$overall["Accuracy"])
  kappa <- as.numeric(cm$overall["Kappa"])
  
  # Macro-F1 (mean of class-wise F1)
  f1_per_class <- cm$byClass[, "F1"]
  macro_f1 <- mean(f1_per_class, na.rm = TRUE)
  
  return(data.frame(
    accuracy = round(acc, 3),
    kappa    = round(kappa, 3),
    macro_f1 = round(macro_f1, 3)
  ))
}

# -----------------------------------------------------------------------------
# STORAGE OBJECTS FOR SUPPLEMENTARY EXPORT
# -----------------------------------------------------------------------------
confusion_list <- list()
classwise_list <- list()

# -----------------------------------------------------------------------------
# Block 1: Evaluate Per-Year Morphology Metrics
# -----------------------------------------------------------------------------
eval_results <- map_dfr(data_years, function(yr) {
  
  cat(sprintf("\nEvaluating Morph3 predictions for year: %s\n", yr))
  
  pred_path  <- file.path(data_dir, paste0("predicted_MORPH3_probs_", yr, ".csv"))
  truth_path <- file.path(data_dir, paste0("dataCARBON_", yr, ".csv"))
  
  # Skip year if files are missing
  if (!file.exists(pred_path) || !file.exists(truth_path)) {
    message(sprintf("Skipping year %s (missing file)", yr))
    return(NULL)
  }
  
  # ------------------------------------------------------------
  # Load predictions (probabilities)
  # ------------------------------------------------------------
  df_pred <- read_csv(pred_path, show_col_types = FALSE)
  
  # ------------------------------------------------------------
  # Load ground-truth morphology labels
  # ------------------------------------------------------------
  df_truth <- read_csv(truth_path, show_col_types = FALSE) %>%
    select(gee_id, morph3)
  
  # Merge truth + predicted probabilities
  df_eval <- inner_join(df_truth, df_pred, by = "gee_id") %>%
    filter(!is.na(morph3))
  
  if (nrow(df_eval) == 0) return(NULL)
  
  # ------------------------------------------------------------
  # Predicted class = highest probability among P_* columns.
  # The "P_" prefix is stripped so the predicted label matches the
  # unprefixed morph3 ground-truth spelling (see fix note at top).
  # ------------------------------------------------------------
  prob_cols <- grep("^P_", names(df_eval), value = TRUE)
  
  df_eval$morph3_pred <- sub(
    "^P_", "",
    colnames(df_eval[, prob_cols])[max.col(df_eval[, prob_cols])]
  )
  
  actual    <- df_eval$morph3
  predicted <- df_eval$morph3_pred
  
  # ------------------------------------------------------------
  # SAFETY CHECK: Confusion matrix requires >=2 classes
  # ------------------------------------------------------------
  n_classes <- length(unique(actual))
  
  if (n_classes < 2) {
    warning(sprintf(
      "Year %s skipped: only one morphology class (%s) found.",
      yr, unique(actual)
    ))
    
    # Still return basic info for record
    return(data.frame(
      accuracy = NA,
      kappa = NA,
      macro_f1 = NA,
      year = yr,
      n = nrow(df_eval)
    ))
  }
  
  # ------------------------------------------------------------
  # Confusion Matrix (Supplementary Table Sx)
  # ------------------------------------------------------------
  cm <- caret::confusionMatrix(
    factor(predicted, levels = levels(factor(actual))),
    factor(actual)
  )
  
  cm_df <- as.data.frame(cm$table)
  names(cm_df) <- c("Predicted", "Reference", "Count")
  
  cm_df <- cm_df %>%
    mutate(
      year = yr,
      n_test = nrow(df_eval)
    )
  
  confusion_list[[yr]] <<- cm_df
  
  # ------------------------------------------------------------
  # Class-wise Precision/Recall/F1 (Supplementary Table Sx)
  # ------------------------------------------------------------
  byclass <- as.data.frame(cm$byClass)
  
  classwise_df <- byclass %>%
    rownames_to_column("Class") %>%
    select(
      Class,
      Precision,
      Recall,
      F1
    ) %>%
    mutate(
      year = yr,
      n_test = nrow(df_eval)
    )
  
  classwise_list[[yr]] <<- classwise_df
  
  # ------------------------------------------------------------
  # Overall yearly metrics (Accuracy, Macro-F1, Kappa)
  # ------------------------------------------------------------
  met <- evaluate_metrics_multiclass(actual, predicted)
  
  met$year <- yr
  met$n <- nrow(df_eval)
  
  return(met)
})

# -----------------------------------------------------------------------------
# Export 1: Overall yearly performance
# -----------------------------------------------------------------------------
write_csv(eval_results, out_perf)

cat("\nSaved Morphology performance summary:\n")
cat(out_perf, "\n")

print(eval_results)

# -----------------------------------------------------------------------------
# Export 2: Yearly confusion matrices (Table Sx)
# -----------------------------------------------------------------------------
confusion_df <- bind_rows(confusion_list)

write_csv(confusion_df, out_cm)

cat("\nSaved yearly confusion matrices:\n")
cat(out_cm, "\n")

# -----------------------------------------------------------------------------
# Export 3: Yearly class-wise metrics (Table Sx)
# -----------------------------------------------------------------------------
classwise_df <- bind_rows(classwise_list)

write_csv(classwise_df, out_cls)

cat("\nSaved yearly class-wise Precision/Recall/F1:\n")
cat(out_cls, "\n")

# -----------------------------------------------------------------------------
# Block 2: Visualizations (Appendix Figures)
# -----------------------------------------------------------------------------

cat("\nPlotting Morph3 performance by year...\n")
print(plot_morph3_performance_by_year())

cat("\nPlotting Morph3 probability distributions...\n")
print(plot_morph3_prob_distribution())

cat("\nPlotting Morph3 spatial probability map (default mono_Ea)...\n")
print(plot_morph3_spatial_map(target_class = "P_mono_Ea"))

# Boxplots per year (probabilities by true class)
for (yr in data_years) {
  cat(sprintf("\nPlotting Morph3 probability boxplots for year %s...\n", yr))
  p <- plot_morph3_boxplot_yearly(yr)
  if (!is.null(p)) print(p)
}

cat("\nAll Morph3 evaluation outputs generated successfully.\n")