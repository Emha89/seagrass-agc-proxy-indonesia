# =============================================================================
# 13_agb_data.R
# Joins the predicted percent cover to the yearly tables (dataAGB_<year>.csv)
# and defines prepare_agb_training_data(), which keeps the samples with a
# valid AGB point estimate. Sourced by 14_agb_tuning.R and 15_agb_model.R.
#
# Inputs   result/dataPA_<year>.csv, result/predicted_tSPC_<year>.csv
# Outputs  result/dataAGB_<year>.csv
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
# FUNCTION: Combined AGB training data (all years)
# -----------------------------------------------------------------------------
prepare_agb_training_data <- function(
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