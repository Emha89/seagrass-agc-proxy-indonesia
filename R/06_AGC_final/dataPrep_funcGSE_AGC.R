# =============================================================================
# dataPrep_funcGSE_AGC2.R
# Merge AGC reference data (dataPA) with predicted tAGB, and prepare AGC
# training data with optional CI-width and outlier filtering.
#
# CONFIRMED VERSION: rf_model_AGC3.R sources this file by name
# (dataPrep_funcGSE_AGC2.R), which is the direct evidence that this is
# the data-prep version actually used -- not the alternative draft that
# explicitly re-merges carbon_index and morph3 P_* columns from
# dataCARBON_<year>.csv / predicted_MORPH3_probs_<year>.csv.
#
# This version instead only CHECKS whether carbon_index and the morph3
# P_* columns are already present in dataPA_<year>.csv, and warns if not
# -- it does not merge them in itself. This means it depends entirely on
# dataPA_<year>.csv already being the pass-2 output of
# dataPrep_func_PA.R (i.e. that script must have been re-run after
# 02_morphology/rf_applyModel_MORPH3.R and
# 05_carbon_index/carbon_scaling_join_SP_finalDell.R have both produced
# their outputs). See the pipeline ordering note in dataPrep_func_PA.R
# for the full multi-pass cycle.
#
# Pipeline position: stage 6 of 6 (final). Reads dataPA_<year>.csv (from
# 01_occurrence_PA, pass 2) and predicted_tAGB_<year>.csv (from 04_AGB).
#
# Provides two alternative training-data-prep functions -- the model
# script calls prepare_training_data_agc_noFilter() (see the equivalent
# AGB functions in 04_AGB/dataPrep_funcGSE_AGB2.R for the same
# withFilter / noFilter pattern, including the same open question about
# the outlier AND-logic in detect_outliers_agc()).
#
# Column names compositio and AGC_CIwidt appear to be 10-character
# truncated field names (a common ESRI shapefile field-name limit
# artifact), not typos -- kept as-is.
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
# PREDICTED tAGB PATHS TO MERGE
# -----------------------------------------------------------------------------
paths <- list(
  "2017" = here("result", "predicted_tAGB_2017.csv"),
  "2018" = here("result", "predicted_tAGB_2018.csv"),
  "2019" = here("result", "predicted_tAGB_2019.csv"),
  "2020" = here("result", "predicted_tAGB_2020.csv"),
  "2021" = here("result", "predicted_tAGB_2021.csv"),
  "2022" = here("result", "predicted_tAGB_2022.csv"),
  "2023" = here("result", "predicted_tAGB_2023.csv")
)

# -----------------------------------------------------------------------------
# FUNCTION: Merge Predicted tAGB into AGC Reference
# Checks (does not merge) Morph3 probability predictors for Stage-B AGC
# regression, and checks (does not merge) carbon_index -- both are
# expected to already be present via dataPA_<year>.csv's pass-2 origin.
# -----------------------------------------------------------------------------
merge_and_export_agc <- function(tagb_path, year_str) {
  
  pa_path    <- file.path(output_dir, paste0("dataPA_", year_str, ".csv"))
  out_path <- file.path(output_dir, paste0("dataAGC_", year_str, ".csv"))
  
  # ---------------------------------------------------------------------------
  # Validate required files
  # ---------------------------------------------------------------------------
  if (!file.exists(pa_path)) {
    cat(sprintf("File dataPA_%s.csv not found. Skipping.\n", year_str))
    return(NULL)
  }
  
  if (!file.exists(tagb_path)) {
    cat(sprintf("File predicted_tAGB_%s.csv not found. Skipping.\n", year_str))
    return(NULL)
  }
  
  # ---------------------------------------------------------------------------
  # Load PA reference data (already includes carbon_index, if pass-2)
  # ---------------------------------------------------------------------------
  df_pa <- read_csv(pa_path, show_col_types = FALSE)
  
  if ("carbon_index" %in% names(df_pa)) {
    cat(sprintf("carbon_index already present in dataPA_%s.csv\n", year_str))
  } else {
    warning(sprintf("carbon_index NOT found in dataPA_%s.csv\n", year_str))
  }
  
  # ---------------------------------------------------------------------------
  # Load predicted tAGB outputs
  # ---------------------------------------------------------------------------
  df_tagb <- read_csv(tagb_path, show_col_types = FALSE) %>%
    select(gee_id, PA_prob, PA_pred, tSPC_pred, tAGB_pred)
  
  # ---------------------------------------------------------------------------
  # Merge PA + predicted tAGB
  # ---------------------------------------------------------------------------
  df <- df_pa %>%
    left_join(df_tagb, by = "gee_id")
  
  # ---------------------------------------------------------------------------
  # Check Morph3 predictors already exist in dataPA_YYYY
  # ---------------------------------------------------------------------------
  morph_cols <- c("P_mixed_short_plus_mono_short",
                  "P_mixed_long",
                  "P_mono_Ea")
  
  missing_morph <- setdiff(morph_cols, names(df_pa))
  
  if (length(missing_morph) == 0) {
    cat(sprintf("Morph3 probabilities already present in dataPA_%s.csv\n", year_str))
  } else {
    warning(sprintf("Missing Morph3 columns in dataPA_%s.csv: %s\n",
                    year_str, paste(missing_morph, collapse = ", ")))
  }
  
  # ---------------------------------------------------------------------------
  # Report missing predictions
  # ---------------------------------------------------------------------------
  n_missing <- sum(is.na(df$tAGB_pred))
  
  cat(sprintf("dataAGC_%s created. Missing tAGB predictions: %d\n",
              year_str, n_missing))
  
  # ---------------------------------------------------------------------------
  # Save output
  # ---------------------------------------------------------------------------
  write_csv(df, out_path)
  cat(sprintf("Saved: %s\n", out_path))
}

# -----------------------------------------------------------------------------
# EXECUTE MERGE PER YEAR
# -----------------------------------------------------------------------------
purrr::iwalk(paths, merge_and_export_agc)

# -----------------------------------------------------------------------------
# FUNCTION: Detect outliers by top-N GSE features (IQR logic)
#
# NOTE: outlier_flag uses AND across all top features -- a row is only
# flagged as an outlier if it is outside the IQR bounds on every one of
# the top N features simultaneously (a very strict definition). Same open
# question as detect_outliers_agb() in 04_AGB/dataPrep_funcGSE_AGB2.R --
# confirm whether this is intentional or whether OR (|) was meant instead.
# -----------------------------------------------------------------------------
detect_outliers_agc <- function(df, top_n_features = 10, factor = 0.5) {
  gse_features <- df %>%
    select(starts_with("GSE_")) %>%
    select(where(is.numeric)) %>%
    colnames()
  
  if (length(gse_features) < 5) {
    stop("Not enough GSE_* features found in the dataset.")
  }
  
  rf_model <- randomForest(
    x = df[, gse_features],
    y = df$AGC_pred,
    ntree = 300,
    importance = TRUE
  )
  
  importance_df <- as.data.frame(randomForest::importance(rf_model)) %>%
    tibble::rownames_to_column("feature") %>%
    arrange(desc(`%IncMSE`))
  
  top_features <- head(importance_df$feature, top_n_features)
  
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
# FUNCTION: Prepare Combined AGC Training Data (WITH filtering CI + Outlier)
# -----------------------------------------------------------------------------
prepare_training_data_agc_withFilter <- function(
    input_dir = output_dir,
    years = c("2017", "2018", "2019", "2020", "2021", "2022", "2023")
) {
  df_all <- purrr::map_dfr(years, function(yr) {
    path <- file.path(input_dir, paste0("dataAGC_", yr, ".csv"))
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
    filter(!is.na(AGC_pred), AGC_pred > 0, !is.na(AGC_CIwidt))
  
  ci_threshold <- quantile(df_all$AGC_CIwidt, 0.75, na.rm = TRUE)
  df_all <- df_all %>% filter(AGC_CIwidt <= ci_threshold)
  
  df_flagged <- detect_outliers_agc(df_all, factor = 0.5)
  
  df_clean <- df_flagged %>%
    filter(!outlier_flag) %>%
    select(-outlier_flag)
  
  cat(sprintf("Data AGC with CI & outlier filtering: %d rows remaining\n", nrow(df_clean)))
  
  return(df_clean)
}

# -----------------------------------------------------------------------------
# FUNCTION: Prepare Combined AGC Training Data (WITHOUT filtering)
# -----------------------------------------------------------------------------
prepare_training_data_agc_noFilter <- function(
    input_dir = output_dir,
    years = c("2017", "2018", "2019", "2020", "2021", "2022", "2023")
) {
  df_all <- purrr::map_dfr(years, function(yr) {
    path <- file.path(input_dir, paste0("dataAGC_", yr, ".csv"))
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
    filter(!is.na(AGC_pred), AGC_pred > 0)
  
  cat(sprintf("Data AGC without CI & outlier filtering: %d rows remaining\n", nrow(df_all)))
  
  return(df_all)
}