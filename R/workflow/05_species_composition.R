# =============================================================================
# 05_species_composition.R
# Adds the species percent cover of each field point to the yearly tables.
#
# Inputs   result/dataPA_<year>.csv (first run of step 01),
#          data/species_composition.csv
# Outputs  result/dataSP_<year>.csv
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(here)   # install.packages("here") if you don't have it yet

# -----------------------------------------------------------------------------
# Define paths
# -----------------------------------------------------------------------------
dataPA_dir <- here("result")
seagrass_path <- here("data", "species_composition.csv")
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