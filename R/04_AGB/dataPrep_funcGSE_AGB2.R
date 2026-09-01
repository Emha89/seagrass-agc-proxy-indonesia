# =============================================================================
# dataPrep_func_AGB.R
# Merge AGB reference data (dataPA) with predicted SPC (tSPC_pred), and
# prepare AGB training data with optional CI-width and outlier filtering.
#
# Pipeline position: stage 4 of 6. Reads dataPA_<year>.csv (from
# 01_occurrence_PA -- must already include the AGB_pred reference column
# from the original GT label file) and predicted_tSPC_<year>.csv (from
# 03_SPC). Note: years covered are 2017-2023 only (no 2024), consistent
# with the SPC training range.
#
# Provides two alternative training-data-prep functions -- the next script
# in this folder presumably calls one of them:
#   - prepare_training_data_agb_withFilter(): applies both a CI-width
#     threshold (top 25% widest confidence intervals excluded) and a
#     multi-feature outlier filter
#   - prepare_training_data_agb_noFilter(): only keeps valid AGB_pred > 0,
#     no CI-width or outlier filtering
#
# OPEN QUESTION: detect_outliers_agb() combines the top-N GSE feature
# outlier checks with AND (&), so a row is only flagged as an outlier if
# it is outside the IQR bounds on every one of the top N features
# simultaneously. This is a very strict definition (few rows will ever
# qualify) -- confirm whether this is intentional or whether OR (|) was
# meant instead.
#
# Column names compositio and AGB_CIwidt appear to be 10-character
# truncated field names (a common ESRI shapefile field-name limit
# artifact from an earlier export step), not typos -- kept as-is.
#
# Data availability: this script's inputs are not included in this
# repository.
# =============================================================================

# -----------------------------------------------------------------------------
# LIBRARIES
# -----------------------------------------------------------------------------
library(dplyr)
library(readr)
library(randomForest)
library(purrr)
library(here)   # install.packages("here") if you don't have it yet

# -----------------------------------------------------------------------------
# SETUP: Paths and Output Directory
# -----------------------------------------------------------------------------
output_dir <- here("result")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# SPC PREDICTION PATHS TO MERGE (BASED ON SPC PIPELINE)
# -----------------------------------------------------------------------------
paths <- list(
  "2017" = here("result", "predicted_tSPC_2017.csv"),
  "2018" = here("result", "predicted_tSPC_2018.csv"),
  "2019" = here("result", "predicted_tSPC_2019.csv"),
  "2020" = here("result", "predicted_tSPC_2020.csv"),
  "2021" = here("result", "predicted_tSPC_2021.csv"),
  "2022" = here("result", "predicted_tSPC_2022.csv"),
  "2023" = here("result", "predicted_tSPC_2023.csv")
)

# -----------------------------------------------------------------------------
# FUNCTION: Merge Predicted SPC into AGB Reference (PA + tSPC from SPC file)
# -----------------------------------------------------------------------------
merge_and_export_agb <- function(spc_path, year_str) {
  # Load AGB reference (dataPA) as the main basis
  pa_path <- file.path(output_dir, paste0("dataPA_", year_str, ".csv"))
  if (!file.exists(pa_path)) {
    cat(sprintf("File dataPA_%s.csv not found. Skipping.\n", year_str))
    return(NULL)
  }
  
  # Load all columns from dataPA (main data basis)
  df_pa <- read_csv(pa_path, show_col_types = FALSE)
  
  # Load predicted SPC (only the additional columns needed)
  df_spc <- read_csv(spc_path, show_col_types = FALSE) %>%
    select(gee_id, PA_prob, PA_pred, tSPC_pred)
  
  # Join by gee_id (dataPA as the main reference)
  df <- df_pa %>%
    left_join(df_spc, by = "gee_id")
  
  # Result summary
  n_missing <- sum(is.na(df$tSPC_pred))
  cat(sprintf("AGB data for year %s merged. Rows without tSPC prediction: %d\n",
              year_str, n_missing))
  
  # Save the full merged result
  write_csv(df, file.path(output_dir, paste0("dataAGB_", year_str, ".csv")))
}

# -----------------------------------------------------------------------------
# EXECUTE MERGE PER YEAR
# -----------------------------------------------------------------------------
purrr::iwalk(paths, merge_and_export_agb)

# -----------------------------------------------------------------------------
# FUNCTION: Detect outliers by top-N GSE features (IQR logic)
# -----------------------------------------------------------------------------
detect_outliers_agb <- function(df, top_n_features = 10, factor = 0.5) {
  gse_features <- df %>%
    select(starts_with("GSE_")) %>%
    select(where(is.numeric)) %>%
    colnames()
  
  if (length(gse_features) < 5) {
    stop("Not enough GSE_* features found in the dataset.")
  }
  
  rf_model <- randomForest(
    x = df[, gse_features],
    y = df$AGB_pred,
    ntree = 300,
    importance = TRUE
  )
  
  importance_df <- as.data.frame(randomForest::importance(rf_model)) %>%
    tibble::rownames_to_column("feature") %>%
    arrange(desc(`%IncMSE`))
  
  top_features <- head(importance_df$feature, top_n_features)
  
  # NOTE: outlier_flag uses AND across all top features -- see the open
  # question at the top of this file.
  df$outlier_flag <- TRUE
  for (f in top_features) {
    Q1 <- quantile(df[[f]], 0.25, na.rm = TRUE)
    Q3 <- quantile(df[[f]], 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    lower <- Q1 - factor * IQR_val
    upper <- Q3 + factor * IQR_val
    df$outlier_flag <- df$outlier_flag & (df[[f]] < lower | df[[f]] > upper)
  }
  
  return(df)
}

# -----------------------------------------------------------------------------
# FUNCTION: Prepare Combined AGB Training Data (WITH filtering CI + Outlier)
# -----------------------------------------------------------------------------
prepare_training_data_agb_withFilter <- function(
    input_dir = output_dir,
    years = c("2017", "2018", "2019", "2020", "2021", "2022", "2023")
) {
  df_all <- purrr::map_dfr(years, function(yr) {
    path <- file.path(input_dir, paste0("dataAGB_", yr, ".csv"))
    if (file.exists(path)) {
      read_csv(path, show_col_types = FALSE) %>%
        mutate(year = as.factor(yr)) %>%
        mutate(across(
          .cols = -c(gee_id, loc, compositio),
          .fns = ~ suppressWarnings(as.numeric(.))
        ))
    } else {
      message(sprintf("Missing file: %s", path))
      NULL
    }
  })
  
  df_all <- df_all %>%
    filter(!is.na(AGB_pred), AGB_pred > 0, !is.na(AGB_CIwidt))
  
  ci_threshold <- quantile(df_all$AGB_CIwidt, 0.75, na.rm = TRUE)
  df_all <- df_all %>% filter(AGB_CIwidt <= ci_threshold)
  
  df_flagged <- detect_outliers_agb(df_all, factor = 0.5)
  
  df_clean <- df_flagged %>%
    filter(!outlier_flag) %>%
    select(-outlier_flag)
  
  cat(sprintf("Data with CI & outlier filtering: %d rows remaining\n", nrow(df_clean)))
  
  return(df_clean)
}

# -----------------------------------------------------------------------------
# FUNCTION: Prepare Combined AGB Training Data (WITHOUT filtering)
# -----------------------------------------------------------------------------
prepare_training_data_agb_noFilter <- function(
    input_dir = output_dir,
    years = c("2017", "2018", "2019", "2020", "2021", "2022", "2023")
) {
  df_all <- purrr::map_dfr(years, function(yr) {
    path <- file.path(input_dir, paste0("dataAGB_", yr, ".csv"))
    if (file.exists(path)) {
      read_csv(path, show_col_types = FALSE) %>%
        mutate(year = as.factor(yr)) %>%
        mutate(across(
          .cols = -c(gee_id, loc, compositio),
          .fns = ~ suppressWarnings(as.numeric(.))
        ))
    } else {
      message(sprintf("Missing file: %s", path))
      NULL
    }
  })
  
  # only keep valid AGB_pred > 0
  df_all <- df_all %>%
    filter(!is.na(AGB_pred), AGB_pred > 0)
  
  cat(sprintf("Data without CI & outlier filtering: %d rows remaining\n", nrow(df_all)))
  
  return(df_all)
}