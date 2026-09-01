# =============================================================================
# rf_applyModel_PA.R
# Apply the trained RF PA model to yearly predictor grids, producing
# predicted_PA_<year>.csv across the full study extent (not just labelled
# GT points).
#
# Pipeline position: stage 1 of 6, after rf_model_PA.R. Its output is
# consumed by rf_evalModel_PA.R and by the PA plotting functions in
# rf_func_vis.R (which expect columns PA_prob, and lon/lat for the spatial
# map -- see the fix note below).
#
# FIX APPLIED: the original version created the sf geometry with
# st_as_sf(..., coords = c("xcoord","ycoord")), which by default drops
# xcoord/ycoord as plain columns, and then dropped the geometry entirely
# with st_drop_geometry() -- so the saved CSV ended up with no coordinate
# columns at all. plot_pa_spatial_map() in rf_func_vis.R reads this exact
# file and expects lon/lat columns, so this has been fixed here by adding
# remove = FALSE (keeps xcoord/ycoord as regular columns) and renaming them
# to lon/lat before saving.
#
# OPEN QUESTION (could not verify without the data): some PA plotting
# functions in rf_func_vis.R (e.g. plot_pa_prob_distribution,
# plot_pa_boxplot_yearly) expect a ground-truth PA column, not just
# PA_prob/PA_pred. Confirm whether your GSE_training_<year>_CSV.csv files
# already contain a PA column for at least some rows -- if not, those
# specific plotting functions won't work against this script's output.
#
# Data availability: this script's inputs (yearly GSE training CSVs, AOI
# shapefile, and the model/predictor-structure files from rf_model_PA.R)
# are not included in this repository.
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
aoi_path <- here("data", "R2_case_study.shp")

output_dir <- here("result")
dir.create(output_dir, showWarnings = FALSE)

paths <- list(
  "2017" = here("data", "GSE_training_2017_CSV.csv"),
  "2018" = here("data", "GSE_training_2018_CSV.csv"),
  "2019" = here("data", "GSE_training_2019_CSV.csv"),
  "2020" = here("data", "GSE_training_2020_CSV.csv"),
  "2021" = here("data", "GSE_training_2021_CSV.csv"),
  "2022" = here("data", "GSE_training_2022_CSV.csv"),
  "2023" = here("data", "GSE_training_2023_CSV.csv"),
  "2024" = here("data", "GSE_training_2024_CSV.csv")
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
# apply_model_regression() in rf_func_reg.R, but is kept as a separate local
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
  # models downstream (produced by the 02_morphology stage -- see the
  # pipeline ordering note in dataPrep_func_PA.R for the same dependency).
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
  # lon/lat here to match what the plotting functions in rf_func_vis.R
  # expect (see the fix note at the top of this file).
  # -------------------------------------------------------------------------
  df_final <- df_final %>% rename(lon = xcoord, lat = ycoord)
  
  out_path <- file.path(output_dir, paste0("predicted_PA_", yr, ".csv"))
  write_csv(df_final, out_path)
  
  cat(sprintf("Saved predictions -> %s (%d rows)\n", out_path, nrow(df_final)))
}

cat("\nAll predictions completed.\n")