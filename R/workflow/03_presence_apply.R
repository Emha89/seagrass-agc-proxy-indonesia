# =============================================================================
# 03_presence_apply.R
# Presence probability (PA_prob) at the field points of every year.
#
# Run this script twice: after step 02, and again after step 08 so that the
# leaf morphology probabilities (P_*) are added to the output.
#
# Inputs   data/gse_points_<year>.csv, data/study_regions.shp,
#          result/final_rf_model_PA.rds, result/rf_predictor_structure_PA.rds,
#          result/predicted_MORPH3_probs_<year>.csv (second run)
# Outputs  result/predicted_PA_<year>.csv
# =============================================================================

library(dplyr)
library(readr)
library(randomForest)
library(sf)
library(purrr)
library(here)   # install.packages("here") if you don't have it yet

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
`%||%` <- function(x, y) if (!is.null(x)) x else y

model_path <- here("result", "final_rf_model_PA.rds")
structure_path <- here("result", "rf_predictor_structure_PA.rds")
aoi_path <- here("data", "study_regions.shp")

output_dir <- here("result")
dir.create(output_dir, showWarnings = FALSE)

paths <- list(
  "2017" = here("data", "gse_points_2017.csv"),
  "2018" = here("data", "gse_points_2018.csv"),
  "2019" = here("data", "gse_points_2019.csv"),
  "2020" = here("data", "gse_points_2020.csv"),
  "2021" = here("data", "gse_points_2021.csv"),
  "2022" = here("data", "gse_points_2022.csv"),
  "2023" = here("data", "gse_points_2023.csv"),
  "2024" = here("data", "gse_points_2024.csv")
)

best_threshold <- 0.7

# -----------------------------------------------------------------------------
# Load model & predictor structure
# -----------------------------------------------------------------------------
model <- readRDS(model_path)
predictor_structure <- readRDS(structure_path)
pa_predictors <- names(predictor_structure$classes)

# -----------------------------------------------------------------------------
# Load AOI shapefile
# -----------------------------------------------------------------------------
if (!file.exists(aoi_path)) stop("AOI shapefile not found.")

aoi_sf <- st_read(aoi_path, quiet = TRUE) %>% st_transform(crs = 4326)

# -----------------------------------------------------------------------------
# Helper: coerce predictor columns to match the classes/factor levels
# recorded when the model was trained (predictor_structure). This mirrors
# apply_model_regression() in rf_regression_functions.R, but is kept as a separate local
# function here because this script needs classification probabilities
# (predict(..., type = "prob")) rather than a plain regression prediction.
# -----------------------------------------------------------------------------
coerce_types <- function(df, classes_list, levels_list) {
  for (nm in names(classes_list)) {
    cls <- classes_list[[nm]][1]
    if (!nm %in% names(df)) next
    if (cls == "numeric") {
      df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    } else if (cls == "integer") {
      df[[nm]] <- suppressWarnings(as.integer(df[[nm]]))
    } else if (cls == "factor") {
      lvls <- levels_list[[nm]] %||% unique(df[[nm]])
      df[[nm]] <- factor(df[[nm]], levels = lvls)
    } else if (cls == "character") {
      df[[nm]] <- as.character(df[[nm]])
    }
  }
  df
}

# -----------------------------------------------------------------------------
# Run prediction per year
# -----------------------------------------------------------------------------
for (yr in names(paths)) {
  cat(sprintf("\nPredicting PA for year: %s\n", yr))
  
  input_path <- paths[[yr]]
  if (!file.exists(input_path)) {
    warning(sprintf("File not found: %s", input_path))
    next
  }
  
  df_raw <- read_csv(input_path, show_col_types = FALSE, guess_max = 100000)
  
  # Check xcoord/ycoord
  if (!all(c("xcoord", "ycoord") %in% names(df_raw))) {
    warning(sprintf("%s: xcoord/ycoord not found. Skipping shapefile filtering.", yr))
    next
  }
  
  # Spatial filter -- remove = FALSE keeps xcoord/ycoord as plain columns
  # alongside the geometry, so they survive st_drop_geometry() below
  df_sf <- df_raw %>%
    filter(!is.na(xcoord), !is.na(ycoord)) %>%
    st_as_sf(coords = c("xcoord", "ycoord"), crs = 4326, remove = FALSE)
  
  df_sf_filtered <- st_join(df_sf, aoi_sf, left = FALSE)
  df_filtered <- df_sf_filtered %>% st_drop_geometry()
  
  cat(sprintf("%s: %d rows inside AOI.\n", yr, nrow(df_filtered)))
  
  if (nrow(df_filtered) == 0) {
    warning(sprintf("%s: No rows within AOI. Skipped.", yr))
    next
  }
  
  # Ensure all required predictors exist
  missing_cols <- setdiff(pa_predictors, names(df_filtered))
  for (mc in missing_cols) df_filtered[[mc]] <- NA
  
  df_pred <- df_filtered %>% select(any_of(pa_predictors))
  df_pred <- coerce_types(df_pred,
                          predictor_structure$classes,
                          predictor_structure$levels)
  
  keep_idx <- complete.cases(df_pred)
  if (sum(keep_idx) == 0) {
    warning(sprintf("%s: No complete rows after NA filtering. Skipped.", yr))
    next
  }
  
  df_pred_cc <- df_pred[keep_idx, ]
  df_final <- df_filtered[keep_idx, ]
  
  probs <- predict(model, newdata = df_pred_cc, type = "prob")[, "1"]
  df_final$PA_prob <- probs
  df_final$PA_pred <- as.integer(probs > best_threshold)
  
  # -------------------------------------------------------------------------
  # Join Morph3 class probability predictors into PA output
  # This ensures predicted_PA_YYYY.csv contains P_* columns for SPC & AGC
  # models downstream (second run, after step 08).
  # -------------------------------------------------------------------------
  morph_path <- file.path(output_dir,
                          paste0("predicted_MORPH3_probs_", yr, ".csv"))
  
  if (file.exists(morph_path)) {
    
    df_morph <- read_csv(morph_path, show_col_types = FALSE) %>%
      select(gee_id, starts_with("P_"))
    
    df_final <- df_final %>%
      left_join(df_morph, by = "gee_id")
    
    cat(sprintf("Morph3 probabilities added into PA output for year %s\n", yr))
    
  } else {
    
    warning(sprintf("Morph3 probability file not found for year %s -- P_* not added", yr))
    
  }
  
  # -------------------------------------------------------------------------
  # Save output (PA + Morph3 probabilities). xcoord/ycoord are renamed to
  # lon/lat here to match the plotting functions in plot_functions.R.
  # -------------------------------------------------------------------------
  df_final <- df_final %>% rename(lon = xcoord, lat = ycoord)
  
  out_path <- file.path(output_dir, paste0("predicted_PA_", yr, ".csv"))
  write_csv(df_final, out_path)
  
  cat(sprintf("Saved predictions -> %s (%d rows)\n", out_path, nrow(df_final)))
}

cat("\nAll predictions completed.\n")