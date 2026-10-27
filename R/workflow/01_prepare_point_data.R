# =============================================================================
# 01_prepare_point_data.R
# Joins the field labels with the GSE, depth and distance-to-land values at
# each field point, and writes one table per survey year.
#
# Run this script twice:
#   first run   before step 02 (base tables)
#   second run  after step 08, to add the field carbon index (step 06) and
#               the leaf morphology probabilities (step 08)
#
# Inputs  (data/)    field_points.csv, gse_points_<year>.csv, study_regions.shp
# Outputs (result/)  dataPA_<year>.csv, summary_loc_year_counts.csv,
#                    data_all_PAsamples_QGIS.csv
# =============================================================================

# -----------------------------------------------------------------------------
# LOAD LIBRARIES
# -----------------------------------------------------------------------------
library(dplyr)
library(readr)
library(sf)
library(purrr)
library(tidyr)
library(here)   # install.packages("here") if you don't have it yet

# -----------------------------------------------------------------------------
# SETUP: Paths and Output Directory
# -----------------------------------------------------------------------------
# Inputs are read from data/ and outputs are written to result/ (see README).
output_dir <- here("result")
label_path <- here("data", "field_points.csv")

# The shapefile needs its sibling files (.shx, .dbf, .prj) in data/.
shapefile_path <- here("data", "study_regions.shp")

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

# -----------------------------------------------------------------------------
# EARLY DIAGNOSTIC: Filter GT data by shapefile before assessing
# -----------------------------------------------------------------------------

cat("Initial diagnostics on GT label data (after AOI filter)...\n")

aoi_sf <- st_read(shapefile_path, quiet = TRUE) %>%
  st_transform(crs = 4326)

gt_df <- read_csv(label_path, show_col_types = FALSE) %>%
  mutate(year = as.integer(year)) %>%
  distinct(gee_id, .keep_all = TRUE)

if (!all(c("xcoord", "ycoord") %in% names(gt_df))) {
  stop("Columns 'xcoord' and 'ycoord' not found in GT data.")
}

gt_sf <- gt_df %>%
  filter(!is.na(xcoord), !is.na(ycoord)) %>%
  st_as_sf(coords = c("xcoord", "ycoord"), crs = 4326)

gt_filtered <- st_join(gt_sf, aoi_sf, left = FALSE) %>%
  st_drop_geometry()

cat(sprintf("GT samples after AOI filtering: %d rows retained.\n", nrow(gt_filtered)))

duplicate_gee_ids <- gt_filtered %>%
  group_by(gee_id) %>%
  filter(n() > 1) %>%
  ungroup()

if (nrow(duplicate_gee_ids) > 0) {
  cat(sprintf("Found %d duplicate gee_id entries in filtered GT data.\n", nrow(duplicate_gee_ids)))
  write_csv(duplicate_gee_ids, file.path(output_dir, "duplicate_GT_gee_ids_filtered.csv"))
} else {
  cat("No duplicate gee_id entries in filtered GT data.\n")
}

gt_years <- gt_filtered %>%
  summarise(
    n_missing = sum(is.na(year)),
    years = paste(sort(unique(na.omit(year))), collapse = ", ")
  )

cat("Year info:\n")
cat(sprintf("- Missing year entries: %d\n", gt_years$n_missing))
cat(sprintf("- Years present: %s\n", gt_years$years))

if ("loc" %in% names(gt_filtered)) {
  gt_locs <- sort(unique(gt_filtered$loc))
  cat("Unique locations in GT data:\n")
  cat(paste("-", gt_locs, collapse = "\n"), "\n")
} else {
  cat("Column 'loc' not found in GT data.\n")
}

cat("GT assessment complete.\n\n")

# -----------------------------------------------------------------------------
# ANALYSIS: Per-Year Summary of Filtered GT Data
# -----------------------------------------------------------------------------

cat("Per-year GT data summary:\n")

gt_filtered_valid <- gt_filtered %>% filter(!is.na(year))
unique_years <- sort(unique(gt_filtered_valid$year))
summary_list <- list()

for (yr in unique_years) {
  gt_year <- gt_filtered_valid %>% filter(year == yr)
  locations <- sort(unique(gt_year$loc))
  n_rows <- nrow(gt_year)
  n_dup <- gt_year %>% group_by(gee_id) %>% filter(n() > 1) %>% nrow()
  
  cat("--------------------------------------------------\n")
  cat(sprintf("Year: %d\n", yr))
  cat(sprintf("Total samples: %d\n", n_rows))
  cat("Locations:\n")
  cat(paste(" -", locations, collapse = "\n"), "\n")
  cat(ifelse(n_dup > 0, sprintf("Duplicates: %d\n", n_dup), "No duplicate gee_id entries.\n"))
  
  summary_list[[as.character(yr)]] <- data.frame(
    year = yr,
    n_samples = n_rows,
    n_locations = length(locations),
    locations = paste(locations, collapse = "; "),
    n_duplicates = n_dup
  )
}

gt_summary_df <- bind_rows(summary_list)

# -----------------------------------------------------------------------------
# DEFINE PREDICTORS (for NA filtering)
# -----------------------------------------------------------------------------
env_vars <- c("distToLand", "rugosity", "slope", "depth", "carbon_index")
wave_vars <- c("elevation", "mean_wave_period", "sig_wave_height")
gse_bands <- paste0("GSE_A", sprintf("%02d", 0:63))
pa_predictors <- c(gse_bands, env_vars, wave_vars)
response_var <- "PA"

# -----------------------------------------------------------------------------
# MERGE GT + IMAGE DATA PER YEAR
# -----------------------------------------------------------------------------

cat("Joining GT with image data by gee_id per year...\n")

merged_yearly_list <- list()

for (yr in unique_years) {
  gt_year <- gt_filtered_valid %>% filter(year == yr)
  image_path <- paths[[as.character(yr)]]
  
  if (is.null(image_path) || !file.exists(image_path)) {
    cat(sprintf("Year %d: Image file not found. Skipping.\n", yr))
    next
  }
  
  image_df <- read_csv(image_path, show_col_types = FALSE)
  gse_cols  <- grep("^GSE_", names(image_df), value = TRUE)
  keep_cols <- unique(c("gee_id", env_vars[env_vars != "carbon_index"],
                        wave_vars, gse_cols, "xcoord", "ycoord"))
  image_df <- image_df %>% select(any_of(keep_cols))
  
  # -------------------------------------------------------------------------
  # Add carbon_index if available (step 06, second run only)
  # -------------------------------------------------------------------------
  carbon_path <- file.path(output_dir, paste0("dataCARBON_", yr, ".csv"))
  if (file.exists(carbon_path)) {
    carbon_df <- read_csv(carbon_path, show_col_types = FALSE) %>%
      select(gee_id, carbon_index)
    
    image_df <- left_join(image_df, carbon_df, by = "gee_id")
    cat(sprintf("Year %d: carbon_index merged (%d rows)\n", yr, nrow(carbon_df)))
  } else {
    cat(sprintf("Year %d: carbon_index file not found, skipping merge.\n", yr))
  }
  
  # -------------------------------------------------------------------------
  # Merge GT + image predictors
  # -------------------------------------------------------------------------
  merged_df <- inner_join(gt_year, image_df, by = "gee_id") %>%
    distinct(gee_id, .keep_all = TRUE) %>%
    rename(year_gt = year)
  
  # -------------------------------------------------------------------------
  # Join morphology probabilities if available (step 08, second run only).
  # Adds: P_mixed_short_plus_mono_short, P_mixed_long, P_mono_Ea
  # These are stored permanently in dataPA_<year>.csv once available, for
  # later use as AGC-stage predictors.
  # -------------------------------------------------------------------------
  morph_path <- file.path(output_dir, paste0("predicted_MORPH3_probs_", yr, ".csv"))
  
  if (file.exists(morph_path)) {
    
    morph_df <- read_csv(morph_path, show_col_types = FALSE) %>%
      select(gee_id, starts_with("P_"))
    
    merged_df <- merged_df %>%
      left_join(morph_df, by = "gee_id")
    
    cat(sprintf("Year %d: Morph3 probabilities merged (%d rows)\n",
                yr, nrow(morph_df)))
    
  } else {
    cat(sprintf("Year %d: Morph3 probability file not found. Skipped.\n", yr))
  }
  
  # -------------------------------------------------------------------------
  # Store result
  # -------------------------------------------------------------------------
  cat(sprintf("Year %d: GT=%d | Matches=%d\n",
              yr, nrow(gt_year), nrow(merged_df)))
  
  merged_yearly_list[[as.character(yr)]] <- merged_df
}

cat("All GT + image data successfully merged.\n\n")

# -----------------------------------------------------------------------------
# DIAGNOSTIC: NA counts per predictor column (per year)
# -----------------------------------------------------------------------------

cat("Checking NA counts per predictor per year...\n")

for (yr in names(merged_yearly_list)) {
  df <- merged_yearly_list[[yr]]
  
  na_summary <- df %>%
    summarise(across(all_of(pa_predictors), ~ sum(is.na(.)), .names = "na_{.col}")) %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "n_NA") %>%
    filter(n_NA > 0) %>%
    arrange(desc(n_NA))
  
  cat(sprintf("Year: %s -- Predictors with NA values:\n", yr))
  print(na_summary, n = Inf)
  cat("--------------------------------------------------\n")
}

cat("Completed NA diagnostics for all years.\n\n")

# -----------------------------------------------------------------------------
# EXPORT PA DATA PER YEAR
# -----------------------------------------------------------------------------
remove_NA_predictors <- FALSE

cat("Exporting PA summary per year...\n")

for (yr in names(merged_yearly_list)) {
  df <- merged_yearly_list[[yr]]
  
  if (!"PA" %in% names(df)) {
    cat(sprintf("Year %s: 'PA' column not found. Skipped.\n", yr))
    next
  }
  
  df <- df %>%
    mutate(PA = as.numeric(PA))
  
  if (remove_NA_predictors) {
    df <- df %>% drop_na(all_of(pa_predictors))
    cat("Rows with NA in predictors removed.\n")
  } else {
    cat("Rows with NA in predictors retained (no filtering).\n")
  }
  
  n_total <- nrow(df)
  n_pa1 <- sum(df$PA == 1, na.rm = TRUE)
  n_pa0 <- sum(df$PA == 0, na.rm = TRUE)
  
  cat("--------------------------------------------------\n")
  cat(sprintf("Year: %s | Total (after NA drop): %d | PA=1: %d | PA=0: %d\n", yr, n_total, n_pa1, n_pa0))
  
  output_path <- file.path(output_dir, paste0("dataPA_", yr, ".csv"))
  write_csv(df, output_path)
  cat(sprintf("Saved: %s\n", output_path))
}

cat("All per-year files successfully saved.\n\n")

# -----------------------------------------------------------------------------
# COMBINE ALL YEARS INTO ONE DATASET
# -----------------------------------------------------------------------------

cat("Combining all yearly dataPA_<year>.csv files...\n")

all_years <- names(merged_yearly_list)
all_files <- file.path(output_dir, paste0("dataPA_", all_years, ".csv"))
existing_files <- all_files[file.exists(all_files)]

data_list <- map(existing_files, ~ read_csv(.x, show_col_types = FALSE))
data_all <- bind_rows(data_list)

cat(sprintf("Total combined rows: %d\n\n", nrow(data_all)))

# -----------------------------------------------------------------------------
# SUMMARY: TOTALS AND PER LOCATION/YEAR STATS
# -----------------------------------------------------------------------------

cat("Computing total and per-location/year summaries...\n")

data_all <- data_all %>%
  mutate(
    PA = as.numeric(PA),
    year_gt = as.integer(year_gt),
    OBJECTID = as.character(OBJECTID)
  )

n_total <- nrow(data_all)
n_pa1 <- sum(data_all$PA == 1, na.rm = TRUE)
n_pa0 <- sum(data_all$PA == 0, na.rm = TRUE)

cat("=====================================\n")
cat("TOTAL SAMPLES (All Years Combined)\n")
cat(sprintf("Total: %d | PA=1: %d | PA=0: %d\n", n_total, n_pa1, n_pa0))
cat("=====================================\n\n")

# tSPC, AGB_pred and AGC_pred come from the field label file (field_points.csv).
required_cols <- c("loc", "year_gt", "PA", "tSPC", "AGB_pred", "AGC_pred")
missing_cols <- setdiff(required_cols, names(data_all))
if (length(missing_cols) > 0) {
  warning(paste("Missing columns in combined data:", paste(missing_cols, collapse = ", ")))
}

data_all <- data_all %>%
  mutate(
    tSPC = as.numeric(tSPC),
    AGB_pred = as.numeric(AGB_pred),
    AGC_pred = as.numeric(AGC_pred)
  )

summary_by_loc_year <- data_all %>%
  group_by(loc, year_gt) %>%
  summarise(
    n_total = n(),
    n_pa1 = sum(PA == 1, na.rm = TRUE),
    n_pa0 = sum(PA == 0, na.rm = TRUE),
    n_tSPC = sum(!is.na(tSPC)),
    n_AGB = sum(!is.na(AGB_pred)),
    n_AGC = sum(!is.na(AGC_pred)),
    .groups = "drop"
  ) %>%
  arrange(loc, year_gt)

print(summary_by_loc_year)

write_csv(summary_by_loc_year, file.path(output_dir, "summary_loc_year_counts.csv"))
cat("Summary table saved as 'summary_loc_year_counts.csv'\n")

# -----------------------------------------------------------------------------
# EXPORT: DATASET READY FOR QGIS VISUALIZATION
# -----------------------------------------------------------------------------

cat("Exporting combined dataset for QGIS visualization...\n")

data_qgis <- data_all %>%
  filter(!is.na(xcoord), !is.na(ycoord)) %>%
  mutate(
    xcoord = as.numeric(xcoord),
    ycoord = as.numeric(ycoord)
  )

output_qgis_path <- file.path(output_dir, "data_all_PAsamples_QGIS.csv")
write_csv(data_qgis, output_qgis_path)

cat("QGIS-ready dataset saved: 'data_all_PAsamples_QGIS.csv'\n")
cat(sprintf("Total points ready for QGIS: %d\n", nrow(data_qgis)))