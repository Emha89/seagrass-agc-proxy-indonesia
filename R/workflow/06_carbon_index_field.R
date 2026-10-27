# =============================================================================
# 06_carbon_index_field.R
# Field carbon index of each point from its species percent cover, and the
# three leaf morphology classes (morph3) used in step 07.
#
#   carbon_base    = sum(species cover x species weight) / 100
#   species_factor = 1 + max_bonus x (n_present - 1) / (n_total - 1)
#   carbon_index   = min(1, carbon_base x species_factor)
#
# The species weights (mean carbon per shoot) are given in Table S2 of the paper.
#
# Inputs   result/dataPA_<year>.csv, result/dataSP_<year>.csv
# Outputs  result/dataCARBON_<year>.csv
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(here)   # install.packages("here") if you don't have it yet

# -----------------------------------------------------------------------------
# 1. PATHS
# -----------------------------------------------------------------------------
result_dir <- here("result")  # all PA & SP files stored here
years <- 2017:2024

# -----------------------------------------------------------------------------
# 2. Species weights (mean carbon per shoot, Table S2 of the paper)
# -----------------------------------------------------------------------------
carbon_species <- c(
  Ea_SPC = 0.829,
  Th_SPC = 0.321,
  Cr_SPC = 0.152,
  Cs_SPC = 0.190,
  Si_SPC = 0.053,
  Hu_SPC = 0.014,
  Ho_SPC = 0.017,
  Hp_SPC = 0.031,
  Tc_SPC = 0.108,
  Hm_SPC = 0.000,
  Hs_SPC = 0.000,
  Hd_SPC = 0.000
)

# maximum richness bonus (1 = up to 2x scaling when all species present)
max_bonus <- 1

# -----------------------------------------------------------------------------
# 3. Morphology 3-class definition (prepared for next RF scripts)
# -----------------------------------------------------------------------------
morph_levels <- c(
  "mixed_short_plus_mono_short",
  "mixed_long",
  "mono_Ea"
)

# -----------------------------------------------------------------------------
# 4. Loop through each year
# -----------------------------------------------------------------------------
for (yr in years) {
  
  cat(sprintf("\nProcessing year %d...\n", yr))
  
  pa_path  <- file.path(result_dir, sprintf("dataPA_%d.csv", yr))
  sp_path  <- file.path(result_dir, sprintf("dataSP_%d.csv", yr))
  out_path <- file.path(result_dir, sprintf("dataCARBON_%d.csv", yr))
  
  # Skip if required files do not exist
  if (!file.exists(pa_path) || !file.exists(sp_path)) {
    cat(sprintf("Missing file for %d -- skipped.\n", yr))
    next
  }
  
  # -----------------------------------------------------------------------------
  # 5. Load data
  # -----------------------------------------------------------------------------
  pa_df <- read_csv(pa_path, show_col_types = FALSE)
  sp_df <- read_csv(sp_path, show_col_types = FALSE)
  
  # Normalize column names
  names(pa_df) <- str_trim(names(pa_df))
  names(sp_df) <- str_trim(names(sp_df))
  
  # Automatically detect species columns based on *_SPC pattern
  species_cols <- grep(
    "^(Ea|Th|Cr|Cs|Si|Hu|Ho|Hp|Tc|Hm|Hs|Hd)_SPC$",
    names(sp_df),
    value = TRUE
  )
  
  if (length(species_cols) == 0) {
    cat(sprintf("No species columns found for %d. Skipped.\n", yr))
    next
  }
  
  cat(sprintf("Detected species columns for %d: %s\n",
              yr, paste(species_cols, collapse = ", ")))
  
  # Total number of species available in the data
  n_total_species <- length(species_cols)
  
  # -----------------------------------------------------------------------------
  # 6. Merge PA + SP datasets by gee_id
  # (sg_morpho is kept so that morph3 can be created)
  # -----------------------------------------------------------------------------
  merged_df <- pa_df %>%
    left_join(
      sp_df %>% select(any_of(c("gee_id", "sg_morpho", species_cols))),
      by = "gee_id"
    )
  
  # Replace NA values in species columns with 0
  merged_df[species_cols][is.na(merged_df[species_cols])] <- 0
  
  # -----------------------------------------------------------------------------
  # 7. Create simplified morphology class (morph3)
  # This is saved here so later RF scripts can directly use it.
  # -----------------------------------------------------------------------------
  if ("sg_morpho" %in% names(merged_df)) {
    
    merged_df <- merged_df %>%
      mutate(
        sg_morpho = as.character(sg_morpho),
        
        # 3-class morphology reclassification
        morph3 = case_when(
          
          # Rule 1: Mixed-long canopy remains its own class
          sg_morpho == "mixed_long" ~ "mixed_long",
          
          # Rule 2: Mono meadow with Enhalus presence becomes mono_Ea
          sg_morpho == "mono" & Ea_SPC > 0 ~ "mono_Ea",
          
          # Rule 3: Remaining mono or mixed_short become short-canopy class
          sg_morpho %in% c("mixed_short", "mono") ~
            "mixed_short_plus_mono_short",
          
          # Rule 4: Other labels ignored
          TRUE ~ NA_character_
        ),
        
        morph3 = factor(morph3, levels = morph_levels)
      )
    
  } else {
    cat("sg_morpho column not found -> morph3 not created.\n")
  }
  
  # -----------------------------------------------------------------------------
  # 8. Calculate Carbon Index
  # -----------------------------------------------------------------------------
  merged_df <- merged_df %>%
    mutate(
      
      # Base carbon (species cover x carbon content)
      carbon_base = rowSums(
        sweep(across(all_of(species_cols)),
              2, carbon_species[species_cols], `*`) / 100,
        na.rm = TRUE
      ),
      
      # Number of species present (cover > 0)
      n_species_present = rowSums(
        across(all_of(species_cols)) > 0,
        na.rm = TRUE
      ),
      
      # Species richness factor (scaled to total available species)
      species_factor = ifelse(
        n_total_species > 1,
        1 + max_bonus * ((n_species_present - 1) /
                           (n_total_species - 1)),
        1
      ),
      
      # Final carbon value after scaling
      carbon_scaled = carbon_base * species_factor,
      
      # Final carbon index (0-1)
      carbon_index = pmin(1, carbon_scaled)
    )
  
  # -----------------------------------------------------------------------------
  # 9. Save yearly output (Carbon Index + Morph3 ready for next scripts)
  # -----------------------------------------------------------------------------
  write_csv(merged_df, out_path)
  
  cat(sprintf("Saved result to: %s\n", out_path))
}

cat("\nAll Carbon Index files successfully created with morph3 included.\n")