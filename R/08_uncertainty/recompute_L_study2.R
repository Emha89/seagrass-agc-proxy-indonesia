# =============================================================================
# recompute_L_study2.R
# Study 2 -- recompute the spatial correlation length L
# used in SD_C,adj = SD_total x L x sqrt(pi x A)
#
# CONTEXT:
#   - The published draft uses L = 30 m, justified as "seagrass patch size".
#   - Legendre (1993): spatial autocorrelation concerns correlation in the
#     model's RESIDUALS, not patterns in the raw measured object. So the
#     correct L is the correlation length of PREDICTION ERROR, not patch size.
#   - Study 2 already has LOLO-CV across 6 regions, so residuals with known
#     coordinates already exist. This gives a real but COARSE picture (only
#     6 discrete points), enough to show L >> 30 m, but not enough to fit a
#     smooth, well-resolved variogram range on its own.
#   - This script does two things: (1) a quick variogram from the existing
#     6-region LOLO-CV residuals as a sanity check, and (2) a finer spatial
#     block cross-validation (blockCV) on the underlying training points to
#     get a larger, better-distributed set of out-of-sample residuals for a
#     more defensible variogram fit.
#
# OPEN ITEM: you also have recompute_L_v2.R and an older "recompute L.R" on
# disk. This file matches the version already confirmed for the project --
# if v2 differs meaningfully, upload it too and we'll compare.
#
# EDIT THIS: only the column names in Section 0 need checking against your
# actual data file -- the paths below now resolve automatically via here(),
# based on wherever this repository sits on your computer.
# =============================================================================

library(dplyr)
library(readr)
library(sf)
library(gstat)
library(blockCV)
library(ranger)
library(ggplot2)
library(here)   # install.packages("here") if you don't have it yet

set.seed(123)

# -----------------------------------------------------------------------------
# 0. Paths and column names
# -----------------------------------------------------------------------------
base_dir     <- here()   # finds this repository's own root folder automatically,
                          # no matter where it's cloned on disk
train_path   <- file.path(base_dir, "data", "carbon_scaling_join_SP_finalDell.csv")  # training data with coords + AGC + predictors
lolo_resid_path <- file.path(base_dir, "result", "agc_lolo_residuals_by_region.csv") # if you already saved per-region LOLO residuals

coord_x_col  <- "xcoord"
coord_y_col  <- "ycoord"
response_col <- "AGC_pred"          # observed/reference AGC
predictor_cols <- c(paste0("GSE_A", sprintf("%02d", 0:63)),
                     "carbon_index", "tSPC_pred", "tAGB_pred",
                     "P_mixed_short_plus_mono_short", "P_mixed_long", "P_mono_Ea")

out_dir <- file.path(base_dir, "result_L_recompute")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# PART A -- Quick check using the existing 6-region LOLO-CV residuals
# =============================================================================
# Expected columns in lolo_resid_path: region, xcoord, ycoord, observed, predicted
# (one row per held-out point, from your existing LOLO-CV run)

cat("=== PART A: variogram from existing LOLO-CV residuals (6 regions) ===\n")

if (file.exists(lolo_resid_path)) {

  lolo_df <- read_csv(lolo_resid_path, show_col_types = FALSE) %>%
    mutate(residual = observed - predicted) %>%
    filter(!is.na(xcoord), !is.na(ycoord), !is.na(residual))

  cat(sprintf("Loaded %d LOLO-CV residuals across %d regions.\n",
              nrow(lolo_df), n_distinct(lolo_df$region)))

  lolo_sf <- st_as_sf(lolo_df, coords = c("xcoord", "ycoord"), crs = 4326) %>%
    st_transform(3857)  # metres-based CRS for distance calculations

  v_lolo <- variogram(residual ~ 1, data = lolo_sf)
  print(v_lolo)

  # Fit a spherical model to extract the range (= L)
  fit_lolo <- tryCatch(
    fit.variogram(v_lolo, vgm(psill = var(lolo_df$residual, na.rm = TRUE),
                               model = "Sph", range = 50000, nugget = 0)),
    error = function(e) { message("Fit failed: ", e$message); NULL }
  )

  if (!is.null(fit_lolo)) {
    cat("\nFitted variogram model (Part A, coarse):\n")
    print(fit_lolo)
    L_partA <- fit_lolo$range[2]
    cat(sprintf("\nEstimated L from 6-region LOLO-CV residuals: %.0f m\n", L_partA))
  } else {
    cat("\nCould not fit a variogram model with only 6 regions -- expected,")
    cat("\nthis part is a sanity check only. See Part B for the fuller estimate.\n")
    L_partA <- NA
  }

  p1 <- ggplot(v_lolo, aes(x = dist, y = gamma)) +
    geom_point(size = 3) +
    labs(title = "Variogram: Study 2 LOLO-CV residuals (6 regions)",
         x = "Distance (m)", y = "Semivariance") +
    theme_minimal(base_size = 13)
  print(p1)
  ggsave(file.path(out_dir, "variogram_partA_lolo_6region.png"), p1, width = 7, height = 5)

} else {
  cat(sprintf("File not found: %s\n", lolo_resid_path))
  cat("Skipping Part A. If you have per-region LOLO-CV predictions saved\n")
  cat("elsewhere, point lolo_resid_path at that file instead.\n")
  L_partA <- NA
}

# =============================================================================
# PART B -- Spatial block cross-validation on the training data for a finer,
#           better-distributed set of out-of-sample residuals
# =============================================================================
cat("\n=== PART B: spatial block CV on training data (finer resolution) ===\n")

train_df <- read_csv(train_path, show_col_types = FALSE) %>%
  filter(!is.na(.data[[coord_x_col]]), !is.na(.data[[coord_y_col]]),
         !is.na(.data[[response_col]]))

cat(sprintf("Loaded %d training records for spatial block CV.\n", nrow(train_df)))

train_sf <- st_as_sf(train_df, coords = c(coord_x_col, coord_y_col), crs = 4326)

# -----------------------------------------------------------------------------
# B1. Define spatial blocks. Start with 20 km; adjust based on how many folds
# you get -- too few points per fold makes RF unstable, too large a block
# size approaches LOLO-CV again.
# -----------------------------------------------------------------------------
block_size_m <- 20000  # 20 km, EDIT if needed

sb <- cv_spatial(
  x = train_sf,
  size = block_size_m,
  k = 10,
  selection = "random",
  iteration = 50,
  progress = FALSE
)

cat(sprintf("Created %d spatial folds using %.0f km blocks.\n",
            length(unique(sb$folds_ids)), block_size_m / 1000))

train_df$fold <- sb$folds_ids

# -----------------------------------------------------------------------------
# B2. Run Random Forest across folds, collect out-of-sample predictions
# -----------------------------------------------------------------------------
rf_formula <- as.formula(paste(response_col, "~", paste(predictor_cols, collapse = " + ")))

results_list <- list()

for (f in sort(unique(train_df$fold))) {

  train_fold <- train_df %>% filter(fold != f)
  test_fold  <- train_df %>% filter(fold == f)

  if (nrow(test_fold) < 5) next  # skip folds that are too small to be meaningful

  rf_fit <- ranger(
    formula = rf_formula,
    data = train_fold,
    num.trees = 500,
    respect.unordered.factors = TRUE
  )

  preds <- predict(rf_fit, data = test_fold)$predictions

  results_list[[as.character(f)]] <- tibble(
    fold = f,
    xcoord = test_fold[[coord_x_col]],
    ycoord = test_fold[[coord_y_col]],
    observed = test_fold[[response_col]],
    predicted = preds,
    residual = test_fold[[response_col]] - preds
  )

  cat(sprintf("Fold %d: n=%d, RMSE=%.2f\n",
              f, nrow(test_fold),
              sqrt(mean((test_fold[[response_col]] - preds)^2))))
}

block_resid_df <- bind_rows(results_list)
write_csv(block_resid_df, file.path(out_dir, "study2_blockCV_residuals.csv"))

overall_rmse <- sqrt(mean(block_resid_df$residual^2, na.rm = TRUE))
cat(sprintf("\nOverall spatial block CV RMSE: %.2f gC m-2\n", overall_rmse))
cat("This RMSE is the honest SD_model input -- compare it to the internal\n")
cat("Monte Carlo RMSE currently used in the thesis.\n")

# -----------------------------------------------------------------------------
# B3. Fit the variogram on the out-of-sample residuals
# -----------------------------------------------------------------------------
block_resid_sf <- st_as_sf(block_resid_df, coords = c("xcoord", "ycoord"), crs = 4326) %>%
  st_transform(3857)

v_block <- variogram(residual ~ 1, data = block_resid_sf, cutoff = 200000, width = 5000)
print(v_block)

fit_block <- fit.variogram(
  v_block,
  vgm(psill = var(block_resid_df$residual, na.rm = TRUE),
      model = "Sph", range = 50000, nugget = 0)
)

cat("\nFitted variogram model (Part B, spatial block CV residuals):\n")
print(fit_block)
L_partB <- fit_block$range[2]
cat(sprintf("\nEstimated L from spatial block CV residuals: %.0f m (%.1f km)\n",
            L_partB, L_partB / 1000))

p2 <- ggplot(v_block, aes(x = dist, y = gamma)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(se = FALSE, colour = "steelblue") +
  labs(title = "Variogram: Study 2 spatial block CV residuals",
       subtitle = sprintf("Fitted range (L) = %.0f m", L_partB),
       x = "Distance (m)", y = "Semivariance") +
  theme_minimal(base_size = 13)
print(p2)
ggsave(file.path(out_dir, "variogram_partB_blockCV.png"), p2, width = 7, height = 5)

# =============================================================================
# PART C -- Recompute the uncertainty formula with the new L
# =============================================================================
cat("\n=== PART C: recompute SD_C,adj and N_effective with the new L ===\n")

A_study2 <- NA  # EDIT: total study area in m2 for Study 2's 6 regions, if this
                # correction is applied at the Study 2 level (check your own Eq. 3-5)

L_old <- 30
L_new <- L_partB  # use the finer, better-resolved estimate

cat(sprintf("L (old, patch size assumption):        %d m\n", L_old))
cat(sprintf("L (new, from spatial block CV residuals): %.0f m\n", L_new))
cat(sprintf("Ratio (new / old): %.1fx\n", L_new / L_old))
cat("\nApply this L to your Eq. 3-5 uncertainty formula (SD_C,adj = SD_total x L x sqrt(pi*A))\n")
cat("using SD_model = the spatial block CV RMSE reported above, not the internal Monte Carlo RMSE.\n")

# Save a summary table
summary_out <- tibble(
  parameter = c("L_old_patch_size_m", "L_partA_lolo_6region_m", "L_partB_blockCV_m",
                "RMSE_internal_MC", "RMSE_spatial_blockCV"),
  value = c(L_old, L_partA, L_partB, NA, overall_rmse)
)
write_csv(summary_out, file.path(out_dir, "L_recompute_summary_study2.csv"))
print(summary_out)

cat("\nScript complete. Outputs saved to:", out_dir, "\n")
