# =============================================================================
# dataPrep_mergeSpecies.R
# Attach seagrass species composition data to dataPA_<year>.csv, saving the
# result as dataSP_<year>.csv.
#
# Pipeline position: stage 5 of 6 (prerequisite). Needs dataPA_<year>.csv
# to exist (pass 1 of 01_occurrence_PA/dataPrep_func_PA.R is enough --
# only the "compositio" column is used here, which doesn't change between
# pass 1 and pass 2). Must run before
# carbon_scaling_join_SP_finalDell.R in this same folder, which reads
# dataSP_<year>.csv as one of its two inputs.
#
# Extracted from a larger combined script that also included exploratory
# per-species spectral band plots and a species presence/dominance summary
# report -- neither of those feeds into any other script in the adopted
# pipeline, so they are not included here (kept in the local working copy
# only, consistent with the other exploratory-work-not-included material
# noted in the README).
#
# Data availability: this script's inputs (dataPA_<year>.csv and
# sg_compos.csv, the raw species composition survey data) are not
# included in this repository.
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(here)   # install.packages("here") if you don't have it yet

# -----------------------------------------------------------------------------
# Define paths
# -----------------------------------------------------------------------------
dataPA_dir <- here("result")
seagrass_path <- here("data", "sg_compos.csv")
output_dir <- dataPA_dir

# -----------------------------------------------------------------------------
# Load seagrass species composition data
# -----------------------------------------------------------------------------
seagrass_df <- read_csv(seagrass_path, show_col_types = FALSE)
cat("Seagrass composition data loaded.\n")

# -----------------------------------------------------------------------------
# Merge seagrass data with dataPA_<year>.csv files (2017-2024)
# -----------------------------------------------------------------------------
years <- 2017:2024

for (yr in years) {
  
  cat("Processing year:", yr, "\n")
  dataPA_path <- file.path(dataPA_dir, paste0("dataPA_", yr, ".csv"))
  
  if (!file.exists(dataPA_path)) {
    cat("File not found:", dataPA_path, "- Skipping.\n")
    next
  }
  
  dataPA_df <- read_csv(dataPA_path, show_col_types = FALSE)
  
  if (!"compositio" %in% names(dataPA_df)) {
    cat("Column 'compositio' missing in:", dataPA_path, "- Skipping.\n")
    next
  }
  
  merged_df <- left_join(dataPA_df, seagrass_df, by = "compositio") %>%
    filter(!is.na(compositio))
  
  cat(sprintf("Year %d: Merged %d -> Filtered: %d (non-NA compositio)\n",
              yr, nrow(dataPA_df), nrow(merged_df)))
  
  output_path <- file.path(output_dir, paste0("dataSP_", yr, ".csv"))
  
  # Remove duplicate columns produced by the merge (e.g. .y), and drop
  # the .x suffix from the remaining duplicated names
  merged_df <- merged_df %>%
    select(-matches("\\.y$")) %>%
    rename_with(~str_replace(., "\\.x$", ""), ends_with(".x"))
  
  write_csv(merged_df, output_path)
  cat("Saved:", output_path, "\n\n")
}

cat("Merging complete for all available years.\n")