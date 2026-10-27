# =============================================================================# 02_agc_pixel_lolo_transfer.R
# Leak-free validation of the AGC models in three designs. In every repeat,
# fold and pair, all five proxy models are refitted with the training samples
# only, and all quadrats of one 10 m pixel stay on the same side of a split.
#
#   PART A  Unseen regions (LOLO-CV): each region is held out in turn.
#   PART B  Single-region transfer: models of one region are applied to each
#           of the other five regions.
#   PART C  Tables and figures for PARTS A and B.
#   PART D  Pixel hold-out: 30% of the pixels of each region are held out,
#           ten repeats.
#   PART E  Accuracy of the percent cover model (samples with seagrass present)
#           and of the AGB model (target AGB_pred), random rows vs
#           pixel-grouped, ten repeats.
#
# Every design compares the proxy workflow, a GSE-only model and a mean-value
# baseline, and reports the bias of the regional mean AGC. The AGC and AGB
# models are trained on the field point estimates (AGC_pred, AGB_pred). Every
# AGC model uses the settings of Table S3 (900 trees, mtry 5, node size 9,
# sampling fraction 0.6).
#
# Used in the paper: Table 3, Figure 4, Tables S1, S6, S9 and S10.
# 03_agc_blocks_surveys_crossfit.R reuses STEPS 0 to 6 of this script.
#
# Model structure (same predictors as the scripts in R/workflow/)
#   PA      : depth + distToLand + GSE             -> PA_prob
#   MORPH3  : depth + distToLand + GSE             -> P_* (3 class probabilities)
#   SPC     : depth + distToLand + GSE + P_*       -> tSPC_pred
#   CINDEX  : depth + distToLand + GSE             -> carbon_index (predicted)
#   AGB     : depth + distToLand + GSE             -> tAGB_pred   (target AGB_pred)
#   AGC     : GSE + env + all proxies above        -> AGC         (target AGC_pred)
#
# Inputs (result/)
#   training_data_for_GEE_PA.csv, training_data_for_GEE_MORPH3_probs.csv,
#   training_data_for_GEE_PCT.csv, training_data_for_GEE_CINDEX.csv,
#   dataAGB_<year>.csv, dataAGC_<year>.csv, dataPA_<year>.csv,
#   predicted_indexC_<year>.csv, grid_rf_best_{PA,MORPH3,PCT,CINDEX,AGB}.csv
#
# Outputs (result/, suffix _pointTarget)
#   PART A  agc_lolocv_pointTarget.csv, agc_lolocv_table_pointTarget.csv
#   PART B  agc_cross_site_pointTarget.csv, agc_cross_site_matrix_pointTarget.csv
#   PART C  agc_lolocv_pointTarget.png, agc_cross_site_heatmap_pointTarget.png
#   PART D  agc_nestedMC_iterations_pointTarget.csv,
#           agc_nestedMC_summary_pointTarget.csv, agc_nestedMC_boxplot_pointTarget.png
#   PART E  proxy_stage_MC_iterations_pointTarget.csv,
#           proxy_stage_MC_summary_pointTarget.csv
#   checkpoint_pointTarget/   one .rds per finished fold or repeat
#
# How to run
#   1. Run once with test_mode <- TRUE (1 fold or repeat per part, 50 trees).
#      Outputs get "_TEST".
#   2. Set test_mode <- FALSE and run again. Finished folds and repeats are
#      saved as checkpoints, so a stopped run can be started again.
# =============================================================================


# -----------------------------------------------------------------------------
# STEP 0: Libraries, paths and run settings
# -----------------------------------------------------------------------------
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(randomForest)
library(ggplot2)

# Paths: open seagrass-agc-proxy-indonesia.Rproj in RStudio (or start R in the
# repository root), so that here() returns the repository folder.
library(here)

result_dir <- here("result")   # input tables and outputs (not distributed, see README)

# TRUE  = quick check (1 fold or iteration per part, 50 trees)
# FALSE = full run
test_mode <- TRUE

run_lolo       <- TRUE   # PART A
run_cross_site <- TRUE   # PART B
run_nested_mc  <- TRUE   # PART D
run_proxy_mc   <- TRUE   # PART E

min_samples <- 10        # minimum rows for a region or a proxy model
seed_base   <- 42
n_iter_mc   <- 10        # PART D and E iterations
test_frac   <- 0.3       # PART D and E: share of pixels (or rows) for testing
agc_years   <- 2017:2023 # years with AGC field data

# AGC settings (Table S3 of the paper), fixed for every AGC model
agc_hp <- list(ntree = 900, mtry = 5, nodesize = 9, sample_fraction = 0.6)
if (test_mode) agc_hp$ntree <- 50

suffix   <- if (test_mode) "_pointTarget_TEST" else "_pointTarget"
ckpt_dir <- file.path(result_dir, paste0("checkpoint", suffix))
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

t_start <- Sys.time()
cat(sprintf("Point-target validation started (%s mode)\n",
            if (test_mode) "TEST" else "FULL"))


# -----------------------------------------------------------------------------
# STEP 1: Predictor definitions
# -----------------------------------------------------------------------------
env_vars  <- c("depth", "distToLand")
gse_vars  <- paste0("GSE_A", sprintf("%02d", 0:63))
base_vars <- c(env_vars, gse_vars)            # PA, MORPH3, CINDEX and AGB models

morph_classes <- c("mixed_short_plus_mono_short", "mixed_long", "mono_Ea")
morph_pred    <- paste0("P_", morph_classes)

spc_vars <- c(env_vars, gse_vars, morph_pred) # SPC model

agc_pred_vars       <- c("PA_prob", "tSPC_pred", "tAGB_pred", "carbon_index")
agc_predictors      <- c(gse_vars, env_vars, agc_pred_vars, morph_pred)
gse_only_predictors <- c(gse_vars, env_vars)  # baseline without proxies

agc_response <- "AGC_pred"   # field point estimate of AGC
agb_response <- "AGB_pred"   # field point estimate of AGB


# -----------------------------------------------------------------------------
# STEP 2: Location lookup (gee_id -> loc) from the annual dataPA files
# -----------------------------------------------------------------------------
pa_files <- list.files(result_dir, pattern = "^dataPA_[0-9]{4}\\.csv$",
                       full.names = TRUE)
if (length(pa_files) == 0) stop("No dataPA_YYYY.csv files found in ", result_dir)

loc_lookup <- map_dfr(pa_files, function(p) {
  read_csv(p, show_col_types = FALSE,
           col_types = cols(.default = col_character()),
           col_select = c(gee_id, loc))
}) %>%
  filter(!is.na(gee_id), !is.na(loc)) %>%
  distinct(gee_id, loc)

conflicts <- loc_lookup %>% count(gee_id) %>% filter(n > 1)
if (nrow(conflicts) > 0) {
  print(head(conflicts))
  stop("Some gee_id values belong to more than one location. Check dataPA files.")
}
cat(sprintf("Location lookup ready: %d points\n", nrow(loc_lookup)))


# -----------------------------------------------------------------------------
# STEP 3: Load the training rows of every stage
# -----------------------------------------------------------------------------

# Pixel key: samples in the same 10 m pixel and year have identical GSE values.
# Values are rounded to 5 decimals to avoid tiny differences from CSV export.
make_pix <- function(df) {
  m <- round(as.matrix(df[, gse_vars]), 5)
  do.call(paste, c(as.data.frame(m), sep = "|"))
}

# 3a-3d. Proxy stages that do not change (same files as before)
load_stage <- function(file, response, predictors, stage) {
  path <- file.path(result_dir, file)
  if (!file.exists(path)) stop(sprintf("[%s] File not found: %s", stage, path))

  df <- read_csv(path, show_col_types = FALSE, guess_max = 100000)

  needed  <- c("gee_id", response, predictors)
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) {
    stop(sprintf("[%s] Missing columns: %s", stage, paste(missing, collapse = ", ")))
  }

  df <- df %>%
    select(all_of(needed)) %>%
    mutate(gee_id = as.character(gee_id)) %>%
    left_join(loc_lookup, by = "gee_id") %>%
    filter(!is.na(loc)) %>%
    filter(if_all(all_of(c(response, predictors)), ~ !is.na(.x)))

  df$pix <- make_pix(df)
  cat(sprintf("[%s] %d training rows loaded\n", stage, nrow(df)))
  df
}

pa_df <- load_stage("training_data_for_GEE_PA.csv", "PA", base_vars, "PA") %>%
  mutate(PA = factor(as.character(PA), levels = c("0", "1")))

morph_df <- load_stage("training_data_for_GEE_MORPH3_probs.csv", "morph3",
                       base_vars, "MORPH3")
unknown_class <- setdiff(unique(as.character(morph_df$morph3)), morph_classes)
if (length(unknown_class) > 0) {
  stop("[MORPH3] Unexpected class names: ", paste(unknown_class, collapse = ", "))
}
morph_df$morph3 <- factor(as.character(morph_df$morph3), levels = morph_classes)

spc_df <- load_stage("training_data_for_GEE_PCT.csv", "tSPC", spc_vars, "SPC")
ci_df  <- load_stage("training_data_for_GEE_CINDEX.csv", "carbon_index", base_vars, "CINDEX")

# 3e. AGB rows with the point estimate as target.
#     Rows with a valid point estimate and 95% interval.
agb_df <- map_dfr(agc_years, function(y) {
  p <- file.path(result_dir, sprintf("dataAGB_%d.csv", y))
  if (!file.exists(p)) {
    cat(sprintf("   [AGB] %s not found, year skipped\n", basename(p)))
    return(NULL)
  }
  read_csv(p, show_col_types = FALSE, guess_max = 100000) %>%
    select(any_of(c("gee_id", "AGB_pred", "AGB_low", "AGB_up", base_vars))) %>%
    mutate(gee_id = as.character(gee_id),
           across(-gee_id, ~ suppressWarnings(as.numeric(.x))))
})

missing_agb <- setdiff(c("gee_id", "AGB_pred", "AGB_low", "AGB_up", base_vars), names(agb_df))
if (length(missing_agb) > 0) {
  stop("[AGB] Missing columns in dataAGB files: ", paste(missing_agb, collapse = ", "))
}

agb_df <- agb_df %>%
  filter(!is.na(AGB_pred), AGB_pred > 0,
         !is.na(AGB_low), !is.na(AGB_up), AGB_low >= 0, AGB_up > AGB_low) %>%
  filter(if_all(all_of(base_vars), ~ !is.na(.x))) %>%
  left_join(loc_lookup, by = "gee_id") %>%
  filter(!is.na(loc)) %>%
  select(gee_id, loc, all_of(agb_response), all_of(base_vars))
agb_df$pix <- make_pix(agb_df)
cat(sprintf("[AGB] %d training rows loaded (target %s)\n", nrow(agb_df), agb_response))

# 3f. AGC rows with the point estimate as target.
#     Rows with a valid point estimate and 95% interval, and with all
#     predictors (3,676 rows in the paper). The in-sample proxy and
#     carbon index columns are read ONLY for this filter. They are not used as
#     predictors: every proxy is rebuilt in each fold.
insample_proxies <- c("PA_prob", "tSPC_pred", "tAGB_pred", morph_pred)

read_agc_year <- function(y) {
  agc_path <- file.path(result_dir, sprintf("dataAGC_%d.csv", y))
  pa_path  <- file.path(result_dir, sprintf("dataPA_%d.csv", y))
  ci_path  <- file.path(result_dir, sprintf("predicted_indexC_%d.csv", y))
  if (!all(file.exists(c(agc_path, pa_path, ci_path)))) return(NULL)

  a <- read_csv(agc_path, show_col_types = FALSE, guess_max = 100000) %>%
    select(any_of(c("gee_id", "loc", "AGC_pred", "AGC_low", "AGC_up",
                    base_vars, insample_proxies))) %>%
    mutate(gee_id = as.character(gee_id), loc = as.character(loc),
           across(-c(gee_id, loc), ~ suppressWarnings(as.numeric(.x))))

  ci_obs <- read_csv(pa_path, show_col_types = FALSE, guess_max = 100000) %>%
    transmute(gee_id = as.character(gee_id),
              ci_obs = suppressWarnings(as.numeric(carbon_index))) %>%
    distinct(gee_id, .keep_all = TRUE)

  ci_pred <- read_csv(ci_path, show_col_types = FALSE, guess_max = 100000) %>%
    transmute(gee_id = as.character(gee_id),
              ci_pred = suppressWarnings(as.numeric(carbon_index_pred))) %>%
    distinct(gee_id, .keep_all = TRUE)

  a %>%
    left_join(ci_obs,  by = "gee_id") %>%
    left_join(ci_pred, by = "gee_id")
}

agc_df <- map_dfr(agc_years, read_agc_year)

needed_agc  <- c("gee_id", "loc", "AGC_pred", "AGC_low", "AGC_up",
                 base_vars, insample_proxies, "ci_obs", "ci_pred")
missing_agc <- setdiff(needed_agc, names(agc_df))
if (length(missing_agc) > 0) {
  stop("[AGC] Missing columns: ", paste(missing_agc, collapse = ", "))
}

agc_df <- agc_df %>%
  filter(!is.na(loc), !is.na(AGC_pred), AGC_pred > 0,
         !is.na(AGC_low), !is.na(AGC_up), AGC_low >= 0, AGC_up > AGC_low) %>%
  filter(if_all(all_of(c(base_vars, insample_proxies, "ci_obs", "ci_pred")),
                ~ !is.na(.x))) %>%
  transmute(gee_id, loc, AGC_pred, across(all_of(base_vars)),
            carbon_index_obs = ci_obs)   # field carbon index, diagnostic only
agc_df$pix <- make_pix(agc_df)

cat(sprintf("[AGC] %d rows loaded (target %s). The paper used 3,676.\n",
            nrow(agc_df), agc_response))
print(count(agc_df, loc))


# -----------------------------------------------------------------------------
# STEP 4: Final hyperparameters of the proxy models (no new tuning)
# -----------------------------------------------------------------------------
read_best <- function(stage) {
  path <- file.path(result_dir, paste0("grid_rf_best_", stage, ".csv"))
  if (!file.exists(path)) stop("Hyperparameter file not found: ", path)
  b <- read_csv(path, show_col_types = FALSE)
  hp <- list(ntree           = b$ntree[1],
             mtry            = b$mtry[1],
             nodesize        = b$nodesize[1],
             sample_fraction = b$sample_fraction[1])
  if (test_mode) hp$ntree <- 50
  hp
}

hp <- list(PA     = read_best("PA"),
           MORPH3 = read_best("MORPH3"),
           PCT    = read_best("PCT"),
           CINDEX = read_best("CINDEX"),
           AGB    = read_best("AGB"))

cat("\nProxy model hyperparameters (fixed):\n")
print(bind_rows(hp, .id = "stage"))
cat("AGC model hyperparameters (fixed, Table S3):\n")
print(as_tibble(agc_hp))


# -----------------------------------------------------------------------------
# STEP 5: Helper functions
# -----------------------------------------------------------------------------

# 5a. Fit one proxy model with fixed hyperparameters. For classification with
#     one class only, a constant probability is returned instead of a forest.
fit_stage <- function(train, response, predictors, hp, seed, stage) {
  if (nrow(train) < min_samples) {
    stop(sprintf("[%s] only %d training rows in this fold", stage, nrow(train)))
  }
  y <- train[[response]]
  if (is.factor(y)) {
    y <- droplevels(y)
    if (nlevels(y) < 2) {
      cat(sprintf("   [%s] one class only (%s); constant probability used\n",
                  stage, levels(y)))
      return(list(type = "const", level = levels(y)))
    }
  }
  set.seed(seed)
  model <- randomForest(
    x        = as.data.frame(train[, predictors]),
    y        = y,
    ntree    = hp$ntree,
    mtry     = min(hp$mtry, length(predictors)),
    nodesize = hp$nodesize,
    sampsize = max(1, floor(hp$sample_fraction * nrow(train))),
    replace  = FALSE
  )
  list(type = "rf", model = model)
}

# 5b. Predict with a fitted proxy model.
#     classes = NULL -> regression (numeric vector)
#     classes = c()  -> classification (probability matrix)
predict_stage <- function(fit, newdata, predictors, classes = NULL) {
  X  <- as.data.frame(newdata[, predictors])
  ok <- complete.cases(X)

  if (is.null(classes)) {
    out <- rep(NA_real_, nrow(X))
    if (any(ok)) out[ok] <- predict(fit$model, newdata = X[ok, , drop = FALSE])
    return(out)
  }

  out <- matrix(NA_real_, nrow(X), length(classes), dimnames = list(NULL, classes))
  if (!any(ok)) return(out)
  out[ok, ] <- 0
  if (fit$type == "const") {
    out[ok, fit$level] <- 1
  } else {
    p <- predict(fit$model, newdata = X[ok, , drop = FALSE], type = "prob")
    out[ok, colnames(p)] <- p
  }
  out
}

# 5c. Rebuild all proxy predictors for the AGC rows.
#     keep_fn(df) returns the rows that each proxy model may use for training
#     in this fold or iteration. target_df = the AGC rows that receive the
#     new proxy values. The AGB model is trained on the point estimate.
build_fold_proxies <- function(keep_fn, seed, target_df = agc_df) {
  out <- target_df

  # (1) Presence probability
  pa_train <- keep_fn(pa_df)
  f_pa <- fit_stage(pa_train, "PA", base_vars, hp$PA, seed + 1, "PA")
  out$PA_prob <- predict_stage(f_pa, out, base_vars, classes = c("0", "1"))[, "1"]

  # (2) Morphology probabilities
  mo_train <- keep_fn(morph_df)
  f_mo <- fit_stage(mo_train, "morph3", base_vars, hp$MORPH3, seed + 2, "MORPH3")
  p_mo <- predict_stage(f_mo, out, base_vars, classes = morph_classes)
  for (k in seq_along(morph_classes)) out[[morph_pred[k]]] <- p_mo[, morph_classes[k]]

  # (3) Percent cover. The P_* inputs of the SPC training rows come from this
  #     fold's morphology model, as in the workflow (R/workflow/09_cover_model.R).
  spc_train <- keep_fn(spc_df)
  p_mo_spc  <- predict_stage(f_mo, spc_train, base_vars, classes = morph_classes)
  for (k in seq_along(morph_classes)) spc_train[[morph_pred[k]]] <- p_mo_spc[, morph_classes[k]]
  f_spc <- fit_stage(spc_train, "tSPC", spc_vars, hp$PCT, seed + 3, "SPC")
  out$tSPC_pred <- predict_stage(f_spc, out, spc_vars)

  # (4) Carbon index, clipped to 0-1 as in R/workflow/12_carbon_index_apply.R
  ci_train <- keep_fn(ci_df)
  f_ci <- fit_stage(ci_train, "carbon_index", base_vars, hp$CINDEX, seed + 4, "CINDEX")
  out$carbon_index <- pmin(pmax(predict_stage(f_ci, out, base_vars), 0), 1)

  # (5) Above-ground biomass (point-estimate target)
  agb_train <- keep_fn(agb_df)
  f_agb <- fit_stage(agb_train, agb_response, base_vars, hp$AGB, seed + 5, "AGB")
  out$tAGB_pred <- predict_stage(f_agb, out, base_vars)

  cat(sprintf("   Proxy training rows | PA %d | MORPH3 %d | SPC %d | CINDEX %d | AGB %d\n",
              nrow(pa_train), nrow(mo_train), nrow(spc_train),
              nrow(ci_train), nrow(agb_train)))

  out %>% filter(if_all(all_of(agc_predictors), ~ !is.na(.x)))
}

# 5d. Random forest with fixed settings (every AGC model)
fit_rf_fixed <- function(train, predictors, response, hp_set, seed) {
  set.seed(seed)
  randomForest(
    x        = as.data.frame(train[, predictors]),
    y        = train[[response]],
    ntree    = hp_set$ntree,
    mtry     = min(hp_set$mtry, length(predictors)),
    nodesize = hp_set$nodesize,
    sampsize = max(1, floor(hp_set$sample_fraction * nrow(train))),
    replace  = FALSE
  )
}

# 5e. Accuracy metrics
#     R2   = squared correlation (as caret::R2); NA for a constant
#     MEC  = model efficiency, 1 - SSE/SST (negative when worse than the mean)
#     Bias = mean prediction minus mean observation; Bias_pct relative to the mean
calc_metrics <- function(obs, pred) {
  ok   <- !is.na(obs) & !is.na(pred)
  obs  <- obs[ok]
  pred <- pred[ok]
  sse  <- sum((obs - pred)^2)
  sst  <- sum((obs - mean(obs))^2)
  r2   <- if (length(unique(pred)) > 1 && length(unique(obs)) > 1) cor(pred, obs)^2 else NA_real_
  tibble(RMSE     = sqrt(mean((obs - pred)^2)),
         MAE      = mean(abs(obs - pred)),
         R2       = r2,
         MEC      = if (sst > 0) 1 - sse / sst else NA_real_,
         mean_obs = mean(obs),
         mean_pred = mean(pred),
         Bias_pct = 100 * (mean(pred) - mean(obs)) / mean(obs))
}

# 5f. Score several prediction vectors against the same observations
score_models <- function(obs, preds) {
  map_dfr(names(preds), function(mn) {
    calc_metrics(obs, preds[[mn]]) %>% mutate(Model = mn, .before = 1)
  })
}

# 5g. Pixel-grouped test set, stratified by region. Pixels are sorted first,
#     so the draw does not depend on row order.
draw_test_pix_sorted <- function(df, seed) {
  set.seed(seed)
  df %>%
    distinct(loc, pix) %>%
    arrange(loc, pix) %>%
    group_by(loc) %>%
    filter(pix %in% sample(pix, max(1, round(test_frac * n())))) %>%
    pull(pix)
}

# 5h. Small helpers
mins_since <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))
q_lo <- function(x) quantile(x, 0.025, na.rm = TRUE)
q_hi <- function(x) quantile(x, 0.975, na.rm = TRUE)

model_levels <- c("Proxy workflow", "GSE bands only", "Mean-value baseline")


# -----------------------------------------------------------------------------
# STEP 6: Locations
# -----------------------------------------------------------------------------
all_locs <- sort(unique(c(pa_df$loc, morph_df$loc, spc_df$loc,
                          ci_df$loc, agb_df$loc, agc_df$loc)))

valid_locs <- agc_df %>%
  count(loc) %>%
  filter(n >= min_samples) %>%
  pull(loc) %>%
  sort()

cat("\nLocations in proxy and AGC data:\n"); print(all_locs)
cat("Locations used as AGC test sites:\n");   print(valid_locs)


# =============================================================================
# PART A: LOLO-CV with all proxy models refitted in each fold
#   Baseline = mean AGC of the training rows (all other regions), i.e. one
#   default value for an unsampled region.
# =============================================================================
lolo_df <- tibble()

if (run_lolo) {
  cat("\n=====================================================\n")
  cat(" PART A: Full-pipeline LOLO-CV (point-estimate target)\n")
  cat("=====================================================\n")

  lolo_locs    <- if (test_mode) valid_locs[1] else valid_locs
  lolo_results <- list()

  for (i in seq_along(lolo_locs)) {
    test_loc <- lolo_locs[i]
    ckpt     <- file.path(ckpt_dir, paste0("lolo_", test_loc, ".rds"))

    if (file.exists(ckpt)) {
      cat(sprintf("[LOLO] %s loaded from checkpoint\n", test_loc))
      lolo_results[[test_loc]] <- readRDS(ckpt)
      next
    }

    cat(sprintf("\n[LOLO %d/%d] Held-out location: %s\n", i, length(lolo_locs), test_loc))
    t0 <- Sys.time()

    res <- tryCatch({
      fold_seed  <- seed_base + 100 * i
      train_locs <- setdiff(all_locs, test_loc)

      # All proxy models are trained without the held-out location
      fold_df  <- build_fold_proxies(function(d) filter(d, loc %in% train_locs), fold_seed)
      train_df <- filter(fold_df, loc != test_loc)
      test_df  <- filter(fold_df, loc == test_loc)
      if (nrow(test_df) == 0) stop("no test rows after filtering")

      m_px  <- fit_rf_fixed(train_df, agc_predictors,      agc_response, agc_hp, fold_seed + 10)
      m_gse <- fit_rf_fixed(train_df, gse_only_predictors, agc_response, agc_hp, fold_seed + 20)

      preds <- list(
        "Proxy workflow"      = predict(m_px,  newdata = as.data.frame(test_df[, agc_predictors])),
        "GSE bands only"      = predict(m_gse, newdata = as.data.frame(test_df[, gse_only_predictors])),
        "Mean-value baseline" = rep(mean(train_df[[agc_response]]), nrow(test_df))
      )

      score_models(test_df[[agc_response]], preds) %>%
        mutate(TestLocation = test_loc,
               n_train      = nrow(train_df),
               n_test       = nrow(test_df),
               # Diagnostic: predicted vs field carbon index at the unseen site
               R2_cindex_heldout = suppressWarnings(
                 cor(test_df$carbon_index, test_df$carbon_index_obs, use = "complete.obs")^2),
               .before = 1)
    }, error = function(e) {
      message(sprintf("[LOLO] %s failed: %s", test_loc, conditionMessage(e)))
      NULL
    })

    if (!is.null(res)) {
      saveRDS(res, ckpt)
      lolo_results[[test_loc]] <- res
      cat(sprintf("%s | Train %d | Test %d | %.1f min\n",
                  test_loc, res$n_train[1], res$n_test[1], mins_since(t0)))
      print(res %>% select(Model, RMSE, R2, MEC, Bias_pct) %>%
              mutate(across(where(is.double), ~ round(.x, 3))))
    }
  }

  lolo_df <- bind_rows(lolo_results)
  if (nrow(lolo_df) > 0) {
    write_csv(lolo_df, file.path(result_dir, paste0("agc_lolocv", suffix, ".csv")))
  }
}


# =============================================================================
# PART B: Cross-site transfer with all proxy models refitted per source
#   Each single-region workflow (proxies and AGC) sees only its own region and
#   is applied to every region. The source region itself is a training fit,
#   not a test, and is flagged as such.
#   Baseline = mean AGC of the source region.
# =============================================================================
cs_df <- tibble()

if (run_cross_site) {
  cat("\n=====================================================\n")
  cat(" PART B: Cross-site transfer (point-estimate target)\n")
  cat("=====================================================\n")

  cs_sources <- if (test_mode) valid_locs[1] else valid_locs
  cs_results <- list()

  for (i in seq_along(cs_sources)) {
    src  <- cs_sources[i]
    ckpt <- file.path(ckpt_dir, paste0("crossSite_", src, ".rds"))

    if (file.exists(ckpt)) {
      cat(sprintf("[Cross-site] %s loaded from checkpoint\n", src))
      cs_results[[src]] <- readRDS(ckpt)
      next
    }

    cat(sprintf("\n[Cross-site %d/%d] Model source: %s\n", i, length(cs_sources), src))
    t0 <- Sys.time()

    res <- tryCatch({
      fold_seed <- seed_base + 1000 + 100 * i

      fold_df  <- build_fold_proxies(function(d) filter(d, loc == src), fold_seed)
      train_df <- filter(fold_df, loc == src)
      if (nrow(train_df) < min_samples) stop(sprintf("only %d AGC training rows", nrow(train_df)))

      m_px     <- fit_rf_fixed(train_df, agc_predictors,      agc_response, agc_hp, fold_seed + 10)
      m_gse    <- fit_rf_fixed(train_df, gse_only_predictors, agc_response, agc_hp, fold_seed + 20)
      src_mean <- mean(train_df[[agc_response]])

      map_dfr(sort(unique(fold_df$loc)), function(test_loc) {
        test_df <- filter(fold_df, loc == test_loc)
        preds <- list(
          "Proxy workflow"      = predict(m_px,  newdata = as.data.frame(test_df[, agc_predictors])),
          "GSE bands only"      = predict(m_gse, newdata = as.data.frame(test_df[, gse_only_predictors])),
          "Mean-value baseline" = rep(src_mean, nrow(test_df))
        )
        score_models(test_df[[agc_response]], preds) %>%
          mutate(ModelSource = src,
                 TestSite    = test_loc,
                 Type        = if (test_loc == src) "Training fit" else "Transfer",
                 train_n     = nrow(train_df),
                 test_n      = nrow(test_df),
                 .before = 1)
      })
    }, error = function(e) {
      message(sprintf("[Cross-site] %s failed: %s", src, conditionMessage(e)))
      NULL
    })

    if (!is.null(res)) {
      saveRDS(res, ckpt)
      cs_results[[src]] <- res
      cat(sprintf("%s done | %.1f min\n", src, mins_since(t0)))
      print(res %>% filter(Model == "Proxy workflow") %>%
              select(TestSite, Type, RMSE, R2, Bias_pct) %>%
              mutate(across(where(is.double), ~ round(.x, 2))))
    }
  }

  cs_df <- bind_rows(cs_results)
  if (nrow(cs_df) > 0) {
    write_csv(cs_df, file.path(result_dir, paste0("agc_cross_site", suffix, ".csv")))
  }
}


# =============================================================================
# PART C: Tables and figures for PARTS A and B
# =============================================================================

# C1. LOLO table for the paper (one row per held-out region)
if (nrow(lolo_df) > 0) {
  lolo_table <- lolo_df %>%
    mutate(Model = factor(Model, levels = model_levels)) %>%
    arrange(TestLocation, Model) %>%
    select(TestLocation, n_test, Model, RMSE, MAE, R2, MEC, Bias_pct) %>%
    mutate(across(c(RMSE, MAE), ~ round(.x, 2)),
           across(c(R2, MEC), ~ round(.x, 3)),
           Bias_pct = round(Bias_pct, 1)) %>%
    pivot_wider(names_from = Model, values_from = c(RMSE, MAE, R2, MEC, Bias_pct),
                names_glue = "{.value}_{Model}")

  out_lt <- file.path(result_dir, paste0("agc_lolocv_table", suffix, ".csv"))
  write_csv(lolo_table, out_lt)
  cat(sprintf("\nLOLO table saved: %s\n", out_lt))
  print(lolo_df %>%
          select(TestLocation, Model, RMSE, R2, MEC, Bias_pct) %>%
          mutate(across(where(is.double), ~ round(.x, 3))),
        n = Inf, width = Inf)

  # C2. Figure: RMSE and regional bias per held-out region
  p_lolo <- lolo_df %>%
    select(TestLocation, Model, RMSE, Bias_pct) %>%
    pivot_longer(c(RMSE, Bias_pct), names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(recode(Metric, RMSE = "RMSE (gC m-2)", Bias_pct = "Bias of regional mean (%)"),
                           levels = c("RMSE (gC m-2)", "Bias of regional mean (%)")),
           Model  = factor(Model, levels = model_levels)) %>%
    ggplot(aes(x = TestLocation, y = Value, colour = Model, shape = Model)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = c("Proxy workflow" = "black",
                                   "GSE bands only" = "grey50",
                                   "Mean-value baseline" = "#c0504d")) +
    labs(title = "LOLO-CV: AGC at an unseen region", x = "Held-out region",
         y = NULL, colour = NULL, shape = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  print(p_lolo)
  ggsave(file.path(result_dir, paste0("agc_lolocv", suffix, ".png")),
         p_lolo, width = 8, height = 7, dpi = 300)
}

# C3. Cross-site RMSE matrix (proxy workflow). The diagonal is a training fit.
#     The last row adds the LOLO result: a model trained on the five other
#     regions together.
if (nrow(cs_df) > 0) {
  cs_px <- cs_df %>% filter(Model == "Proxy workflow")

  cs_matrix <- cs_px %>%
    mutate(Cell = ifelse(Type == "Training fit",
                         sprintf("(%.2f)", RMSE), sprintf("%.2f", RMSE))) %>%
    select(ModelSource, TestSite, Cell) %>%
    pivot_wider(names_from = TestSite, values_from = Cell)

  if (nrow(lolo_df) > 0) {
    lolo_row <- lolo_df %>%
      filter(Model == "Proxy workflow") %>%
      transmute(ModelSource = "Other five regions (LOLO)", TestSite = TestLocation,
                Cell = sprintf("%.2f", RMSE)) %>%
      pivot_wider(names_from = TestSite, values_from = Cell)
    cs_matrix <- bind_rows(cs_matrix, lolo_row)
  }

  out_mx <- file.path(result_dir, paste0("agc_cross_site_matrix", suffix, ".csv"))
  write_csv(cs_matrix, out_mx)
  cat(sprintf("\nCross-site matrix saved: %s (diagonal in brackets = training fit)\n", out_mx))
  print(cs_matrix, width = Inf)

  cat("\nTransfer to other regions (off-diagonal), by model:\n")
  print(cs_df %>%
          filter(Type == "Transfer") %>%
          group_by(Model) %>%
          summarise(n_pairs   = n(),
                    RMSE_mean = mean(RMSE), RMSE_min = min(RMSE), RMSE_max = max(RMSE),
                    R2_mean   = if (all(is.na(R2))) NA_real_ else mean(R2, na.rm = TRUE),
                    abs_bias_pct_mean = mean(abs(Bias_pct)),
                    .groups = "drop") %>%
          mutate(across(where(is.double), ~ round(.x, 3))),
        width = Inf)

  # C4. Heatmap (proxy workflow)
  p_cs <- cs_px %>%
    mutate(ModelSource = factor(ModelSource, levels = rev(valid_locs)),
           Label = ifelse(Type == "Training fit",
                          sprintf("(%.1f)", RMSE), sprintf("%.1f", RMSE))) %>%
    ggplot(aes(x = TestSite, y = ModelSource, fill = RMSE)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = Label), size = 3.8) +
    scale_fill_gradient(low = "#8fd19e", high = "#e8706f",
                        name = expression("RMSE (gC m"^-2*")")) +
    labs(title = "Cross-site AGC RMSE (proxy workflow, point-estimate target)",
         subtitle = "Values in brackets are training fits (source = test region)",
         x = "Test region", y = "Model source") +
    theme_minimal(base_size = 12)

  print(p_cs)
  ggsave(file.path(result_dir, paste0("agc_cross_site_heatmap", suffix, ".png")),
         p_cs, width = 9, height = 5.5, dpi = 300)
}


# =============================================================================
# PART D: Pixel hold-out (nested, pixel-grouped Monte Carlo within regions)
#   In each iteration:
#     1. Per region, 30% of pixels are drawn as the test set.
#     2. Every proxy model is refitted without the test pixels and samples.
#     3. Three AGC predictions are scored on the same test pixels: proxy
#        workflow, GSE-only and the regional mean of the training rows.
# =============================================================================
mc_df <- tibble()

if (run_nested_mc) {
  cat("\n=====================================================\n")
  cat(" PART D: Nested, pixel-grouped Monte Carlo (S5)\n")
  cat("=====================================================\n")

  mc_iters   <- if (test_mode) 1 else n_iter_mc
  mc_results <- list()

  for (it in seq_len(mc_iters)) {
    ckpt <- file.path(ckpt_dir, sprintf("nestedMC_iter_%02d.rds", it))

    if (file.exists(ckpt)) {
      cat(sprintf("[Nested MC] iteration %d loaded from checkpoint\n", it))
      mc_results[[it]] <- readRDS(ckpt)
      next
    }

    cat(sprintf("\n[Nested MC %d/%d]\n", it, mc_iters))
    t0 <- Sys.time()

    res <- tryCatch({
      it_seed  <- seed_base + 9000 + 100 * it     # same seed as step S4
      test_pix <- draw_test_pix_sorted(agc_df, it_seed)
      test_ids <- agc_df$gee_id[agc_df$pix %in% test_pix]

      keep_fn  <- function(x) filter(x, !(pix %in% test_pix), !(gee_id %in% test_ids))
      fold_df  <- build_fold_proxies(keep_fn, it_seed)
      train_df <- filter(fold_df, !(pix %in% test_pix))
      test_df  <- filter(fold_df,   pix %in% test_pix)
      if (nrow(test_df) == 0) stop("no test rows after filtering")

      m_px  <- fit_rf_fixed(train_df, agc_predictors,      agc_response, agc_hp, it_seed + 10)
      m_gse <- fit_rf_fixed(train_df, gse_only_predictors, agc_response, agc_hp, it_seed + 20)

      reg_mean <- train_df %>%
        group_by(loc) %>%
        summarise(mu = mean(.data[[agc_response]]), .groups = "drop")

      test_df <- test_df %>% left_join(reg_mean, by = "loc")
      test_df$pred_px  <- predict(m_px,  newdata = as.data.frame(test_df[, agc_predictors]))
      test_df$pred_gse <- predict(m_gse, newdata = as.data.frame(test_df[, gse_only_predictors]))

      score_scope <- function(d, scope) {
        score_models(d[[agc_response]],
                     list("Proxy workflow"      = d$pred_px,
                          "GSE bands only"      = d$pred_gse,
                          "Mean-value baseline" = d$mu)) %>%
          mutate(Scope = scope, n_test = nrow(d), .before = 1)
      }

      bind_rows(
        score_scope(test_df, "All regions"),
        map_dfr(sort(unique(test_df$loc)), function(l) score_scope(filter(test_df, loc == l), l))
      ) %>%
        mutate(Iteration = it, n_train = nrow(train_df))
    }, error = function(e) {
      message(sprintf("[Nested MC] iteration %d failed: %s", it, conditionMessage(e)))
      NULL
    })

    if (!is.null(res)) {
      saveRDS(res, ckpt)
      mc_results[[it]] <- res
      overall <- filter(res, Scope == "All regions")
      cat(sprintf("Iteration %d | Train %d | Test %d | %.1f min\n",
                  it, overall$n_train[1], overall$n_test[1], mins_since(t0)))
      print(overall %>% select(Model, RMSE, R2, MEC, Bias_pct) %>%
              mutate(across(where(is.double), ~ round(.x, 3))))
    }
  }

  mc_df <- bind_rows(mc_results)

  if (nrow(mc_df) > 0) {
    write_csv(mc_df, file.path(result_dir, paste0("agc_nestedMC_iterations", suffix, ".csv")))

    mc_summary <- mc_df %>%
      group_by(Scope, Model) %>%
      summarise(n_iter         = n(),
                n_test         = round(mean(n_test)),
                RMSE_mean      = mean(RMSE), RMSE_lower = q_lo(RMSE), RMSE_upper = q_hi(RMSE),
                MAE_mean       = mean(MAE),
                R2_mean        = if (all(is.na(R2))) NA_real_ else mean(R2, na.rm = TRUE),
                R2_lower       = q_lo(R2), R2_upper = q_hi(R2),
                MEC_mean       = mean(MEC, na.rm = TRUE),
                Bias_pct_mean  = mean(Bias_pct),
                Bias_pct_lower = q_lo(Bias_pct), Bias_pct_upper = q_hi(Bias_pct),
                .groups = "drop") %>%
      mutate(Model = factor(Model, levels = model_levels)) %>%
      arrange(Scope != "All regions", Scope, Model) %>%
      mutate(Model = as.character(Model))

    out_sum <- file.path(result_dir, paste0("agc_nestedMC_summary", suffix, ".csv"))
    write_csv(mc_summary, out_sum)
    cat(sprintf("\nNested Monte Carlo summary saved: %s\n", out_sum))
    print(mc_summary %>%
            select(Scope, Model, n_test, RMSE_mean, R2_mean, R2_lower, R2_upper,
                   MEC_mean, Bias_pct_mean) %>%
            mutate(across(where(is.double), ~ round(.x, 3))),
          n = Inf, width = Inf)

    # Paired comparison on the same test pixels: proxy workflow minus GSE-only
    paired <- mc_df %>%
      filter(Scope == "All regions", Model %in% c("Proxy workflow", "GSE bands only")) %>%
      select(Iteration, Model, RMSE, R2) %>%
      pivot_wider(names_from = Model, values_from = c(RMSE, R2))
    if (all(c("R2_Proxy workflow", "R2_GSE bands only") %in% names(paired))) {
      cat(sprintf("\nProxy minus GSE-only (same test pixels): dR2 = %.3f (proxy better in %d of %d), dRMSE = %.2f\n",
                  mean(paired$`R2_Proxy workflow` - paired$`R2_GSE bands only`),
                  sum(paired$`R2_Proxy workflow` > paired$`R2_GSE bands only`), nrow(paired),
                  mean(paired$`RMSE_Proxy workflow` - paired$`RMSE_GSE bands only`)))
    }

    # Figure: RMSE per region across iterations
    p_mc <- mc_df %>%
      mutate(Scope = factor(Scope, levels = c("All regions", sort(valid_locs))),
             Model = factor(Model, levels = model_levels)) %>%
      ggplot(aes(x = Scope, y = RMSE, fill = Model)) +
      geom_boxplot(width = 0.7, outlier.size = 1, position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = c("Proxy workflow" = "grey40",
                                   "GSE bands only" = "grey80",
                                   "Mean-value baseline" = "#e6b8b7")) +
      labs(title = "Nested, pixel-grouped Monte Carlo (within-region accuracy)",
           x = NULL, y = expression("RMSE (gC m"^-2*")"), fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

    print(p_mc)
    ggsave(file.path(result_dir, paste0("agc_nestedMC_boxplot", suffix, ".png")),
           p_mc, width = 9, height = 5.5, dpi = 300)
  }
}


# =============================================================================
# PART E: Monte Carlo of the percent cover and AGB models
#   AGB : target AGB_pred (point estimate)
#   SPC : samples with seagrass present only (tSPC > 0)
#   Random rows = single samples held out. Pixel-grouped = all samples of a pixel kept
#   together; for SPC the P_* inputs then come from a morphology model that
#   did not see the test pixels.
# =============================================================================
proxy_df <- tibble()

if (run_proxy_mc) {
  cat("\n=====================================================\n")
  cat(" PART E: Proxy stages (AGB point target, SPC presence only)\n")
  cat("=====================================================\n")

  proxy_stages <- list(
    AGB         = list(df = agb_df, response = agb_response, predictors = base_vars, hp = hp$AGB),
    SPC_presence = list(df = filter(spc_df, tSPC > 0), response = "tSPC",
                        predictors = spc_vars, hp = hp$PCT)
  )

  px_iters <- if (test_mode) 1 else n_iter_mc

  for (stage in names(proxy_stages)) {
    st <- proxy_stages[[stage]]
    cat(sprintf("\n[Proxy MC] %s: %d rows\n", stage, nrow(st$df)))

    for (mode in c("Random rows", "Pixel-grouped")) {
      tag <- if (mode == "Random rows") "random" else "pixel"

      for (it in seq_len(px_iters)) {
        ckpt <- file.path(ckpt_dir, sprintf("proxyMC_%s_%s_iter_%02d.rds", stage, tag, it))
        if (file.exists(ckpt)) {
          cat(sprintf("[Proxy MC] %s %s iteration %d loaded from checkpoint\n", stage, tag, it))
          next
        }

        t0      <- Sys.time()
        it_seed <- seed_base + 7000 + 100 * it   # same seeds as PART E2 before

        res <- tryCatch({
          d <- st$df

          if (mode == "Random rows") {
            set.seed(it_seed)
            tr    <- sample(nrow(d), floor((1 - test_frac) * nrow(d)))
            train <- d[tr, ]
            test  <- d[-tr, ]
          } else {
            test_pix <- draw_test_pix_sorted(d, it_seed)
            train    <- filter(d, !(pix %in% test_pix))
            test     <- filter(d,   pix %in% test_pix)

            if (stage == "SPC_presence") {
              test_ids <- test$gee_id
              mo_train <- filter(morph_df, !(pix %in% test_pix), !(gee_id %in% test_ids))
              f_mo <- fit_stage(mo_train, "morph3", base_vars, hp$MORPH3, it_seed + 2, "MORPH3")
              p_tr <- predict_stage(f_mo, train, base_vars, morph_classes)
              p_te <- predict_stage(f_mo, test,  base_vars, morph_classes)
              for (k in seq_along(morph_classes)) {
                train[[morph_pred[k]]] <- p_tr[, morph_classes[k]]
                test[[morph_pred[k]]]  <- p_te[, morph_classes[k]]
              }
            }
          }

          fit  <- fit_stage(train, st$response, st$predictors, st$hp, it_seed + 1, stage)
          pred <- predict_stage(fit, test, st$predictors)
          calc_metrics(test[[st$response]], pred) %>%
            select(RMSE, MAE, R2, MEC) %>%
            pivot_longer(everything(), names_to = "Metric", values_to = "Value") %>%
            mutate(Stage = stage, Mode = mode, Iteration = it,
                   n_train = nrow(train), n_test = nrow(test))
        }, error = function(e) {
          message(sprintf("[Proxy MC] %s %s iteration %d failed: %s",
                          stage, tag, it, conditionMessage(e)))
          NULL
        })

        if (!is.null(res)) {
          saveRDS(res, ckpt)
          cat(sprintf("[Proxy MC] %s %s iteration %d | %s | %.1f min\n", stage, tag, it,
                      paste(sprintf("%s %.3f", res$Metric, res$Value), collapse = " | "),
                      mins_since(t0)))
        }
      }
    }
  }
}

proxy_files <- list.files(ckpt_dir, pattern = "^proxyMC_.*\\.rds$", full.names = TRUE)

if (length(proxy_files) > 0) {
  proxy_df <- map_dfr(proxy_files, readRDS)
  write_csv(proxy_df, file.path(result_dir, paste0("proxy_stage_MC_iterations", suffix, ".csv")))

  reported_s1 <- tribble(
    ~Stage,         ~Metric, ~Reported_TableS1,
    "SPC_presence", "R2",     0.57,
    "SPC_presence", "RMSE",  16.97,
    "AGB",          "R2",     0.50,
    "AGB",          "RMSE",  55.47
  )

  proxy_summary <- proxy_df %>%
    group_by(Stage, Mode, Metric) %>%
    summarise(n_iter = n(),
              n_test = round(mean(n_test)),
              mean   = mean(Value, na.rm = TRUE),
              lower  = q_lo(Value),
              upper  = q_hi(Value),
              .groups = "drop") %>%
    left_join(reported_s1, by = c("Stage", "Metric")) %>%
    arrange(Stage, Metric, desc(Mode))

  out_px <- file.path(result_dir, paste0("proxy_stage_MC_summary", suffix, ".csv"))
  write_csv(proxy_summary, out_px)
  cat(sprintf("\nProxy-stage summary saved: %s\n", out_px))
  print(proxy_summary %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)
}

cat(sprintf("\nFinished. Total time: %.1f min\n", mins_since(t_start)))
if (test_mode) {
  cat("Test run OK? Set test_mode <- FALSE and run the full script.\n")
}
