# =============================================================================
# 16_agc_data.R
# Joins the predicted AGB to the yearly tables (dataAGC_<year>.csv) and
# defines prepare_agc_training_data(), which keeps the samples with a valid
# AGC point estimate. Sourced by 17_agc_tuning.R and 18_agc_model.R.
#
# Inputs   result/dataPA_<year>.csv (second run of step 01),
#          result/predicted_tAGB_<year>.csv
# Outputs  result/dataAGC_<year>.csv
# =============================================================================

# -----------------------------------------------------------------------------
# LIBRARIES
# -----------------------------------------------------------------------------
library(dplyr)
library(readr)
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
# Checks that carbon_index and the morphology probabilities are present
# (second run of step 01). 18_agc_model.R then replaces the field carbon
# index with the predicted carbon index.
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
  # Load PA reference data (second run of step 01)
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
# FUNCTION: Combined AGC training data (all years)
# -----------------------------------------------------------------------------
prepare_agc_training_data <- function(
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