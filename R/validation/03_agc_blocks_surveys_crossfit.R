# =============================================================================# 03_agc_blocks_surveys_crossfit.R
# Further leak-free validation designs of the AGC models. Run it after
# 02_agc_pixel_lolo_transfer.R: it reuses STEPS 0 to 6 of that script (data and
# helper functions) and its PART D results (pixel hold-out). The part letters
# continue those of 02 (A to E).
#
#   PART F  New survey of the same region: all samples of one region and survey
#           year are held out. Proxy workflow and GSE-only.
#   PART G  Uncertainty of the reference AGC values: the 95% interval of each
#           sample as an SD, compared with the RMSE of the pixel hold-out.
#   PART H  Proxy models trained on all samples (including the test samples)
#           vs refitted with the training samples only.
#   PART I  Spatial blocks of 100, 250, 500 and 1000 m held out inside each
#           region (GSE-only and mean-value baseline), with the distance from
#           each test sample to the nearest training sample.
#   PART J  Distance from the test samples to the nearest training sample in
#           the new-survey, LOLO-CV and single-region transfer designs.
#   PART K  Proxy workflow in the 250 m block design.
#   PART L  Cross-fitted proxy values for the training samples, in the pixel
#           hold-out and 250 m block designs, and the share of the proxy layers
#           in the predictor importance.
#
# Used in the paper: Tables 3, S4, S7, S8, S9 and S12, Figure 4 and Figure S3b.
#
# Outputs (result/, suffix _ext)
#   PART F  agc_survey_transfer_ext.csv
#   PART G  agc_error_floor_ext.csv
#   PART H  agc_proxy_sensitivity_iterations_ext.csv,
#           agc_proxy_sensitivity_summary_ext.csv
#   PART I  agc_block_iterations_ext.csv, agc_block_summary_ext.csv,
#           agc_block_ext.png
#   PART J  agc_design_distances_ext.csv
#   PART K  agc_block_proxy_iterations_ext.csv, agc_block_proxy_summary_ext.csv
#   PART L  agc_crossfit_iterations_ext.csv, agc_crossfit_summary_ext.csv,
#           agc_crossfit_paired_ext.csv, agc_crossfit_proxy_r2_ext.csv
#   checkpoint_ext/   one .rds per finished survey, repeat or inner group
#
# How to run
#   1. Keep ext_test_mode <- TRUE (1 repeat, 50 trees). Outputs get "_TEST".
#   2. Set ext_test_mode <- FALSE and run again. Finished parts are saved as
#      checkpoints, so a stopped run can be started again.
#   Approximate run time of the full run: PART F 20-40 min, PART H 1-2 h,
#   PART I 20-40 min, PART J a few minutes, PART K 1-2 h. PART L refits all
#   proxy models (cf_k + 1) times per repeat, so it takes about (cf_k + 1)
#   times as long as PART K for each design.
# =============================================================================

# Paths: open seagrass-agc-proxy-indonesia.Rproj in RStudio (or start R in the
# repository root), so that here() returns the repository folder.
library(here)

ext_test_mode <- TRUE

run_survey   <- TRUE   # PART F
run_floor    <- TRUE   # PART G
run_proxysen <- TRUE   # PART H
run_blocks   <- TRUE   # PART I
run_dist     <- TRUE   # PART J
run_blockproxy <- TRUE # PART K
run_crossfit   <- TRUE # PART L

# PART L settings
cf_k       <- 5                                   # inner groups for cross-fitting
cf_designs <- c("Block 250 m", "Pixel-grouped")   # run order; drop one to skip it

# Names of the coordinate columns in dataAGC_YYYY.csv (decimal degrees,
# WGS84: xcoord = longitude, ycoord = latitude)
lon_col <- "xcoord"
lat_col <- "ycoord"


# -----------------------------------------------------------------------------
# STEP 0: Reuse STEPS 0 to 6 of 02_agc_pixel_lolo_transfer.R
#   Only the lines before its PART A are run, so its PARTS A to E are not
#   repeated. test_mode is set by ext_test_mode above.
# -----------------------------------------------------------------------------
setup_file <- here("R", "validation", "02_agc_pixel_lolo_transfer.R")
if (!file.exists(setup_file)) stop("Setup script not found: ", setup_file)

setup_lines <- readLines(setup_file, warn = FALSE)
end_setup   <- grep("^# PART A: LOLO-CV", setup_lines)[1] - 2
if (is.na(end_setup)) stop("PART A not found in ", setup_file)

setup_lines <- setup_lines[seq_len(end_setup)]
setup_lines <- sub("^test_mode <- (TRUE|FALSE)",
                   paste0("test_mode <- ", ext_test_mode), setup_lines)
eval(parse(text = setup_lines))

ext_suffix   <- if (test_mode) "_ext_TEST" else "_ext"
ext_ckpt_dir <- file.path(result_dir, paste0("checkpoint", ext_suffix))
dir.create(ext_ckpt_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- function(name) file.path(result_dir, paste0(name, ext_suffix, ".csv"))
fig_path <- function(name) file.path(result_dir, paste0(name, ext_suffix, ".png"))

t_ext <- Sys.time()


# -----------------------------------------------------------------------------
# STEP 1: Survey year and interval bounds of every AGC row
#   The setup keeps neither. gee_id + pix (GSE values of that year) is unique.
# -----------------------------------------------------------------------------
agc_extra <- map_dfr(agc_years, function(y) {
  p <- file.path(result_dir, sprintf("dataAGC_%d.csv", y))
  if (!file.exists(p)) return(NULL)
  d <- read_csv(p, show_col_types = FALSE, guess_max = 100000) %>%
    select(any_of(c("gee_id", "AGC_low", "AGC_up", lon_col, lat_col, gse_vars))) %>%
    mutate(gee_id = as.character(gee_id),
           across(-gee_id, ~ suppressWarnings(as.numeric(.x))))
  if (!all(gse_vars %in% names(d))) return(NULL)
  for (cc in c(lon_col, lat_col)) if (!cc %in% names(d)) d[[cc]] <- NA_real_
  d$pix <- make_pix(d)
  transmute(d, gee_id, pix, year = y, AGC_low, AGC_up,
            lon = .data[[lon_col]], lat = .data[[lat_col]])
}) %>%
  distinct(gee_id, pix, .keep_all = TRUE)

agc_df <- agc_df %>% left_join(agc_extra, by = c("gee_id", "pix"))
if (any(is.na(agc_df$year))) {
  stop(sprintf("%d AGC rows without a survey year. Check gee_id and GSE columns.",
               sum(is.na(agc_df$year))))
}
cat("\nAGC samples per region and survey year:\n")
print(count(agc_df, loc, year), n = Inf)

# Survey year of every point in the dataPA files (used to drop a whole survey
# from the proxy training data in PART F)
pa_year <- map_dfr(pa_files, function(p) {
  y <- as.integer(sub("^dataPA_([0-9]{4})\\.csv$", "\\1", basename(p)))
  read_csv(p, show_col_types = FALSE,
           col_types = cols(.default = col_character()),
           col_select = c(gee_id)) %>%
    transmute(gee_id, year = y)
}) %>%
  filter(!is.na(gee_id)) %>%
  distinct() %>%
  left_join(loc_lookup, by = "gee_id")

# Within-region reference accuracy (PART D of 02_agc_pixel_lolo_transfer.R)
ref_path <- file.path(result_dir, "agc_nestedMC_summary_pointTarget.csv")
ref_within <- if (file.exists(ref_path)) {
  read_csv(ref_path, show_col_types = FALSE) %>%
    filter(Scope != "All regions") %>%
    transmute(Region = Scope, Model, RMSE_within = RMSE_mean, MEC_within = MEC_mean)
} else {
  cat("agc_nestedMC_summary_pointTarget.csv not found; no within-region reference\n")
  tibble(Region = character(), Model = character(),
         RMSE_within = numeric(), MEC_within = numeric())
}


# =============================================================================
# PART F: Transfer between surveys of the same region
#   Test = all AGC samples of one region in one survey year (n >= min_survey).
#   Training scope:
#     Same region, other years : AGC and proxy models use only the other
#                                surveys of that region
#     All other data           : all regions, without the test survey
#   In both scopes, every proxy model is refitted without any point of the
#   test survey (dataPA year and region) and without the test pixels.
# =============================================================================
min_survey <- 20

if (run_survey) {
  cat("\n=====================================================\n")
  cat(" PART F: Transfer between surveys of the same region\n")
  cat("=====================================================\n")

  survey_tab <- agc_df %>% count(loc, year, name = "n")
  multi_loc  <- survey_tab %>% filter(n >= min_survey) %>% count(loc) %>%
    filter(n >= 2) %>% pull(loc)
  surveys <- survey_tab %>% filter(loc %in% multi_loc, n >= min_survey) %>%
    arrange(loc, year)
  if (test_mode) surveys <- head(surveys, 1)
  cat("Test surveys:\n"); print(surveys)

  sv_results <- list()

  for (si in seq_len(nrow(surveys))) {
    r <- surveys$loc[si]; yb <- surveys$year[si]
    key  <- sprintf("%s_%d", r, yb)
    ckpt <- file.path(ext_ckpt_dir, paste0("survey_", key, ".rds"))
    if (file.exists(ckpt)) {
      cat(sprintf("[Survey] %s loaded from checkpoint\n", key))
      sv_results[[key]] <- readRDS(ckpt)
      next
    }

    cat(sprintf("\n[Survey %d/%d] Test: %s %d\n", si, nrow(surveys), r, yb))
    t0 <- Sys.time()

    survey_ids <- pa_year %>% filter(loc == r, year == yb) %>% pull(gee_id)
    test_pix   <- agc_df %>% filter(loc == r, year == yb) %>% pull(pix)
    test_ids   <- agc_df %>% filter(loc == r, year == yb) %>% pull(gee_id)
    drop_test  <- function(d) filter(d, !(gee_id %in% c(survey_ids, test_ids)),
                                     !(pix %in% test_pix))

    scopes <- list(
      "Same region, other years" = function(d) filter(drop_test(d), loc == r),
      "All other data"           = drop_test
    )

    res <- map_dfr(names(scopes), function(sc) {
      tryCatch({
        s <- seed_base + 40000 + 100 * si + if (sc == "All other data") 50 else 0
        fold_df  <- build_fold_proxies(scopes[[sc]], s)
        test_df  <- filter(fold_df, loc == r, year == yb)
        train_df <- filter(fold_df, !(loc == r & year == yb))
        if (sc == "Same region, other years") train_df <- filter(train_df, loc == r)
        if (nrow(train_df) < min_samples) stop(sprintf("only %d AGC training rows", nrow(train_df)))

        m_px  <- fit_rf_fixed(train_df, agc_predictors,      agc_response, agc_hp, s + 10)
        m_gse <- fit_rf_fixed(train_df, gse_only_predictors, agc_response, agc_hp, s + 20)

        preds <- list(
          "Proxy workflow"      = predict(m_px,  newdata = as.data.frame(test_df[, agc_predictors])),
          "GSE bands only"      = predict(m_gse, newdata = as.data.frame(test_df[, gse_only_predictors])),
          "Mean-value baseline" = rep(mean(train_df[[agc_response]]), nrow(test_df))
        )
        score_models(test_df[[agc_response]], preds) %>%
          mutate(Region = r, TestYear = yb, Scope = sc,
                 TrainYears = paste(sort(unique(train_df$year[train_df$loc == r])), collapse = "+"),
                 n_train = nrow(train_df), n_test = nrow(test_df), .before = 1)
      }, error = function(e) {
        message(sprintf("[Survey] %s, %s failed: %s", key, sc, conditionMessage(e)))
        NULL
      })
    })

    if (nrow(res) > 0) {
      saveRDS(res, ckpt)
      sv_results[[key]] <- res
      cat(sprintf("%s done | %.1f min\n", key, mins_since(t0)))
      print(res %>% select(Scope, Model, RMSE, R2, MEC, Bias_pct) %>%
              mutate(across(where(is.double), ~ round(.x, 3))), width = Inf)
    }
  }

  sv_df <- bind_rows(sv_results)
  if (nrow(sv_df) > 0) {
    sv_df <- sv_df %>% left_join(ref_within, by = c("Region", "Model"))
    write_csv(sv_df, out_path("agc_survey_transfer"))
    cat("Survey transfer table saved\n")
  }
}


# =============================================================================
# PART G: Error floor from the uncertainty of the field AGC estimates
#   target_sd = (AGC_up - AGC_low) / (2 x 1.96), an approximate SD of each
#   field estimate (the intervals are asymmetric, so this is a rough value).
#   If a model predicted the true AGC exactly, its RMSE against the point
#   estimates would still be about sqrt(mean(target_sd^2)).
#   NOTE: the intervals come from field variance (standard error of shoot
#   density) plus model variance of the cover-to-biomass conversion. They leave
#   out the residual error of the AGB model and the OCC uncertainty, so this
#   floor is a LOWER bound.
# =============================================================================
if (run_floor) {
  cat("\n=====================================================\n")
  cat(" PART G: Error floor from field AGC uncertainty\n")
  cat("=====================================================\n")

  floor_rows <- agc_df %>%
    mutate(target_sd = (AGC_up - AGC_low) / (2 * 1.96),
           half_low  = AGC_pred - AGC_low,
           half_up   = AGC_up - AGC_pred)

  summarise_floor <- function(d) {
    d %>% summarise(n = n(),
                    mean_AGC       = mean(AGC_pred),
                    rmse_floor     = sqrt(mean(target_sd^2)),
                    mean_half_low  = mean(half_low),
                    mean_half_up   = mean(half_up),
                    .groups = "drop")
  }

  floor_tab <- bind_rows(
    floor_rows %>% group_by(Region = loc) %>% summarise_floor(),
    floor_rows %>% summarise_floor() %>% mutate(Region = "All regions", .before = 1)
  )

  if (file.exists(ref_path)) {
    ref_all <- read_csv(ref_path, show_col_types = FALSE) %>%
      filter(Model == "GSE bands only") %>%
      transmute(Region = Scope, RMSE_within_GSE = RMSE_mean)
    floor_tab <- floor_tab %>%
      left_join(ref_all, by = "Region") %>%
      mutate(share_of_MSE_CI = rmse_floor^2 / RMSE_within_GSE^2)
  }

  write_csv(floor_tab, out_path("agc_error_floor"))
  cat("Error floor table saved\n")
  print(floor_tab %>% mutate(across(where(is.double), ~ round(.x, 2))), width = Inf)
}


# =============================================================================
# PART H: Proxy models trained on all samples vs refitted (Table S4)
#   Cells (split x proxy fit):
#     Pixel-grouped x refitted   : PART D of 02_agc_pixel_lolo_transfer.R (read)
#     Pixel-grouped x all samples: proxies trained on all rows, built once
#     Random rows   x all samples: same proxies, 30% of rows per region
#     Random rows   x refitted   : proxies refitted without the test rows
#   Same seeds as PART D, so the pixel-grouped cells use the same test pixels.
# =============================================================================
draw_test_rows <- function(df, seed) {
  set.seed(seed)
  df %>%
    distinct(loc, gee_id) %>%
    arrange(loc, gee_id) %>%
    group_by(loc) %>%
    filter(gee_id %in% sample(gee_id, max(1, round(test_frac * n())))) %>%
    pull(gee_id)
}

if (run_proxysen) {
  cat("\n=====================================================\n")
  cat(" PART H: Proxy sensitivity (point-estimate target)\n")
  cat("=====================================================\n")

  ks_iters   <- if (test_mode) 1 else n_iter_mc
  ks_results <- list()

  insample_file <- file.path(ext_ckpt_dir, "proxies_all_samples.rds")
  if (file.exists(insample_file)) {
    agc_insample <- readRDS(insample_file)
  } else {
    cat("Building proxy layers from all samples (once)\n")
    agc_insample <- build_fold_proxies(function(d) d, seed_base + 30000)
    saveRDS(agc_insample, insample_file)
  }

  for (it in seq_len(ks_iters)) {
    ckpt <- file.path(ext_ckpt_dir, sprintf("proxysen_iter_%02d.rds", it))
    if (file.exists(ckpt)) {
      cat(sprintf("[ProxySen] iteration %d loaded from checkpoint\n", it))
      ks_results[[it]] <- readRDS(ckpt)
      next
    }
    cat(sprintf("\n[ProxySen %d/%d]\n", it, ks_iters))
    t0 <- Sys.time()

    res <- tryCatch({
      it_seed  <- seed_base + 9000 + 100 * it          # same as PART D
      test_pix <- draw_test_pix_sorted(agc_df, it_seed)
      test_row <- draw_test_rows(agc_df, it_seed)

      score_cell <- function(fold_df, is_test, split, proxies, seed) {
        train_df <- fold_df[!is_test(fold_df), ]
        test_df  <- fold_df[ is_test(fold_df), ]
        m <- fit_rf_fixed(train_df, agc_predictors, agc_response, agc_hp, seed)
        calc_metrics(test_df[[agc_response]],
                     predict(m, newdata = as.data.frame(test_df[, agc_predictors]))) %>%
          mutate(Split = split, Proxies = proxies, Iteration = it,
                 n_train = nrow(train_df), n_test = nrow(test_df), .before = 1)
      }
      is_test_pix <- function(d) d$pix %in% test_pix
      is_test_row <- function(d) d$gee_id %in% test_row

      fold_rows <- build_fold_proxies(function(d) filter(d, !(gee_id %in% test_row)), it_seed + 500)

      bind_rows(
        score_cell(agc_insample, is_test_pix, "Pixel-grouped", "Trained on all samples", it_seed + 10),
        score_cell(agc_insample, is_test_row, "Random rows",   "Trained on all samples", it_seed + 10),
        score_cell(fold_rows,    is_test_row, "Random rows",   "Refitted",               it_seed + 10)
      )
    }, error = function(e) {
      message(sprintf("[ProxySen] iteration %d failed: %s", it, conditionMessage(e)))
      NULL
    })

    if (!is.null(res)) {
      saveRDS(res, ckpt)
      ks_results[[it]] <- res
      cat(sprintf("Iteration %d | %.1f min\n", it, mins_since(t0)))
      print(res %>% select(Split, Proxies, RMSE, R2, MEC) %>%
              mutate(across(where(is.double), ~ round(.x, 3))))
    }
  }

  ks_df <- bind_rows(ks_results)

  # Add the pixel-grouped x refitted cell from PART D (proxy workflow, all regions)
  partD_path <- file.path(result_dir, "agc_nestedMC_iterations_pointTarget.csv")
  if (file.exists(partD_path)) {
    partD <- read_csv(partD_path, show_col_types = FALSE) %>%
      filter(Scope == "All regions", Model == "Proxy workflow",
             Iteration %in% unique(ks_df$Iteration)) %>%
      transmute(Split = "Pixel-grouped", Proxies = "Refitted", Iteration,
                n_train, n_test, RMSE, MAE, R2, MEC, mean_obs, mean_pred, Bias_pct)
    ks_df <- bind_rows(ks_df, partD)
  }

  if (nrow(ks_df) > 0) {
    write_csv(ks_df, out_path("agc_proxy_sensitivity_iterations"))
    ks_summary <- ks_df %>%
      group_by(Split, Proxies) %>%
      summarise(n_iter = n(),
                RMSE_mean = mean(RMSE), RMSE_lower = q_lo(RMSE), RMSE_upper = q_hi(RMSE),
                R2_mean = mean(R2), R2_lower = q_lo(R2), R2_upper = q_hi(R2),
                MEC_mean = mean(MEC), .groups = "drop")
    write_csv(ks_summary, out_path("agc_proxy_sensitivity_summary"))
    cat("Proxy sensitivity summary saved\n")
    print(ks_summary %>% mutate(across(where(is.double), ~ round(.x, 3))), width = Inf)
  }
}

# =============================================================================
# PART I: Spatial block hold-out inside each region (GSE-only)
#   Coordinates are turned into metres around the centre of each region.
#   Blocks of b x b metres are made per region (all survey years together),
#   and 30% of the blocks of every region are held out, as in PART D.
#   nn_m = distance from a test sample to the nearest training sample of the
#   same region. The pixel-grouped design of PART D is added as a reference
#   (its accuracy from agc_nestedMC_summary_pointTarget.csv, its distances
#   computed here with the same test pixels).
# =============================================================================
block_sizes <- c(100, 250, 500, 1000)   # block edge length (m)
n_iter_blk  <- if (test_mode) 1 else 10

nn_dist <- function(test, train) {
  d <- rep(NA_real_, nrow(test))
  for (l in unique(test$loc)) {
    i_te <- which(test$loc == l)
    tr   <- train[train$loc == l, c("x_m", "y_m")]
    if (nrow(tr) == 0) next
    d[i_te] <- vapply(i_te, function(i)
      sqrt(min((tr$x_m - test$x_m[i])^2 + (tr$y_m - test$y_m[i])^2)), numeric(1))
  }
  d
}

if (run_blocks) {
  cat("\n=====================================================\n")
  cat(" PART I: Spatial block hold-out (GSE-only)\n")
  cat("=====================================================\n")

  if (all(is.na(agc_df$lon)) || all(is.na(agc_df$lat))) {
    cat(sprintf("No coordinates in columns '%s'/'%s'; PART I skipped\n", lon_col, lat_col))
  } else {
    xy_df <- agc_df %>%
      filter(!is.na(lon), !is.na(lat)) %>%
      group_by(loc) %>%
      mutate(x_m = (lon - mean(lon)) * 111320 * cos(mean(lat) * pi / 180),
             y_m = (lat - mean(lat)) * 110540) %>%
      ungroup()
    cat(sprintf("%d of %d AGC rows have coordinates\n", nrow(xy_df), nrow(agc_df)))
    if (any(abs(xy_df$lon) > 180) || any(abs(xy_df$lat) > 90)) {
      stop("Coordinates are not in decimal degrees. Check lon_col and lat_col.")
    }

    # I1. Distances in the pixel-grouped design (no new models)
    pix_nn <- map_dfr(seq_len(n_iter_blk), function(it) {
      it_seed  <- seed_base + 9000 + 100 * it          # same as PART D
      test_pix <- draw_test_pix_sorted(agc_df, it_seed)
      te <- filter(xy_df, pix %in% test_pix)
      te$nn_m <- nn_dist(te, filter(xy_df, !(pix %in% test_pix)))
      bind_rows(tibble(Scope = "All regions", median_nn_m = median(te$nn_m, na.rm = TRUE)),
                te %>% group_by(Scope = loc) %>%
                  summarise(median_nn_m = median(nn_m, na.rm = TRUE), .groups = "drop")) %>%
        mutate(Iteration = it)
    })

    # I2. Block hold-out
    blk_results <- list()
    for (b in block_sizes) {
      for (it in seq_len(n_iter_blk)) {
        key  <- sprintf("block_%04dm_iter_%02d", b, it)
        ckpt <- file.path(ext_ckpt_dir, paste0(key, ".rds"))
        if (file.exists(ckpt)) {
          blk_results[[key]] <- readRDS(ckpt)
          next
        }
        t0 <- Sys.time()
        res <- tryCatch({
          s <- seed_base + 50000 + b + it
          d <- xy_df %>% mutate(block = paste(loc, floor(x_m / b), floor(y_m / b), sep = "_"))
          set.seed(s)
          test_blocks <- d %>%
            distinct(loc, block) %>% arrange(loc, block) %>%
            group_by(loc) %>%
            filter(block %in% sample(block, max(1, round(test_frac * n())))) %>%
            pull(block)
          test_df  <- filter(d,   block %in% test_blocks)
          train_df <- filter(d, !(block %in% test_blocks))

          m <- fit_rf_fixed(train_df, gse_only_predictors, agc_response, agc_hp, s + 20)
          reg_mean <- train_df %>% group_by(loc) %>%
            summarise(mu = mean(.data[[agc_response]]), .groups = "drop")
          test_df <- test_df %>% left_join(reg_mean, by = "loc")
          test_df$pred <- predict(m, newdata = as.data.frame(test_df[, gse_only_predictors]))
          test_df$nn_m <- nn_dist(test_df, train_df)

          sc <- function(x, scope) {
            score_models(x[[agc_response]],
                         list("GSE bands only" = x$pred, "Mean-value baseline" = x$mu)) %>%
              mutate(Scope = scope, n_test = nrow(x),
                     n_blocks = n_distinct(x$block),
                     median_nn_m = median(x$nn_m, na.rm = TRUE), .before = 1)
          }
          bind_rows(sc(test_df, "All regions"),
                    map_dfr(sort(unique(test_df$loc)), function(l) sc(filter(test_df, loc == l), l))) %>%
            mutate(Block_m = b, Iteration = it, n_train = nrow(train_df), .before = 1)
        }, error = function(e) {
          message(sprintf("[Blocks] %s failed: %s", key, conditionMessage(e)))
          NULL
        })
        if (!is.null(res)) {
          saveRDS(res, ckpt)
          blk_results[[key]] <- res
          ov <- filter(res, Scope == "All regions", Model == "GSE bands only")
          cat(sprintf("%4d m, iteration %d | test %d | median nn %.0f m | R2 %.2f | MEC %.2f | %.1f min\n",
                      b, it, ov$n_test, ov$median_nn_m, ov$R2, ov$MEC, mins_since(t0)))
        }
      }
    }

    blk_df <- bind_rows(blk_results)
    if (nrow(blk_df) > 0) {
      write_csv(blk_df, out_path("agc_block_iterations"))

      blk_summary <- blk_df %>%
        group_by(Block_m, Scope, Model) %>%
        summarise(n_iter = n(), n_test = round(mean(n_test)),
                  median_nn_m = mean(median_nn_m),
                  RMSE_mean = mean(RMSE), RMSE_lower = q_lo(RMSE), RMSE_upper = q_hi(RMSE),
                  R2_mean = if (all(is.na(R2))) NA_real_ else mean(R2, na.rm = TRUE),
                  MEC_mean = mean(MEC), MEC_lower = q_lo(MEC), MEC_upper = q_hi(MEC),
                  abs_bias_pct = mean(abs(Bias_pct)), .groups = "drop") %>%
        mutate(Design = paste0("Block ", Block_m, " m"), .before = 1)

      # Pixel-grouped reference (PART D accuracy + distances from I1)
      if (file.exists(ref_path)) {
        pix_ref <- read_csv(ref_path, show_col_types = FALSE) %>%
          filter(Model %in% c("GSE bands only", "Mean-value baseline")) %>%
          transmute(Scope, Model, n_iter, n_test, RMSE_mean, RMSE_lower, RMSE_upper,
                    R2_mean, MEC_mean, abs_bias_pct = abs(Bias_pct_mean)) %>%
          left_join(pix_nn %>% group_by(Scope) %>%
                      summarise(median_nn_m = mean(median_nn_m), .groups = "drop"),
                    by = "Scope") %>%
          mutate(Design = "Pixel-grouped (PART D)", Block_m = 10, .before = 1)
        blk_summary <- bind_rows(pix_ref, blk_summary)
      }

      write_csv(blk_summary, out_path("agc_block_summary"))
      cat("Block hold-out summary saved\n")
      print(blk_summary %>% filter(Scope == "All regions") %>%
              select(Design, Model, n_test, median_nn_m, RMSE_mean, R2_mean, MEC_mean) %>%
              mutate(across(where(is.double), ~ round(.x, 2))), n = Inf, width = Inf)

      p_blk <- blk_summary %>%
        mutate(Scope = factor(Scope, levels = c("All regions", sort(setdiff(unique(Scope), "All regions"))))) %>%
        ggplot(aes(x = Block_m, y = MEC_mean, colour = Model)) +
        geom_hline(yintercept = 0, colour = "grey70") +
        geom_line() + geom_point(size = 1.8) +
        facet_wrap(~ Scope, scales = "free_y") +
        scale_x_log10(breaks = c(10, block_sizes), labels = c("pixel", block_sizes)) +
        scale_colour_manual(values = c("GSE bands only" = "#0072B2", "Mean-value baseline" = "#D55E00")) +
        labs(x = "Size of held-out unit (m)", y = "MEC", colour = NULL) +
        theme_minimal(base_size = 11) +
        theme(legend.position = "bottom")
      ggsave(fig_path("agc_block"), p_blk, width = 9, height = 6, dpi = 300, bg = "white")
    }
  }
}


# =============================================================================
# Helper for PARTS J, K and L: coordinates in metres around each region centre
# =============================================================================
make_xy_df <- function(df) {
  if (all(is.na(df$lon)) || all(is.na(df$lat))) return(NULL)
  out <- df %>%
    filter(!is.na(lon), !is.na(lat)) %>%
    group_by(loc) %>%
    mutate(x_m = (lon - mean(lon)) * 111320 * cos(mean(lat) * pi / 180),
           y_m = (lat - mean(lat)) * 110540) %>%
    ungroup()
  if (any(abs(out$lon) > 180) || any(abs(out$lat) > 90)) {
    stop("Coordinates are not in decimal degrees. Check lon_col and lat_col.")
  }
  out
}


# =============================================================================
# PART J: Distance from test samples to the nearest training sample (no models)
#   Great-circle distance (m) on the coordinates.
#   Survey transfer : test survey vs the other surveys of the same region
#   LOLO-CV         : held-out region vs all other regions
#   Transfer        : test region vs the source region
# =============================================================================
nn_haversine <- function(lon1, lat1, lon2, lat2) {
  to_rad <- pi / 180
  la2 <- lat2 * to_rad; lo2 <- lon2 * to_rad; cl2 <- cos(la2)
  vapply(seq_along(lon1), function(i) {
    la1 <- lat1[i] * to_rad; lo1 <- lon1[i] * to_rad
    a <- sin((la2 - la1) / 2)^2 + cos(la1) * cl2 * sin((lo2 - lo1) / 2)^2
    2 * 6371000 * asin(sqrt(min(a)))
  }, numeric(1))
}

if (run_dist) {
  cat("\n=====================================================\n")
  cat(" PART J: Distances to the nearest training sample\n")
  cat("=====================================================\n")

  xy_d <- make_xy_df(agc_df)
  if (is.null(xy_d)) {
    cat("No coordinates; PART J skipped\n")
  } else {
    dist_row <- function(design, test, source, te, tr) {
      d <- nn_haversine(te$lon, te$lat, tr$lon, tr$lat)
      tibble(Design = design, Test = test, Source = source, n_test = nrow(te),
             median_nn_m = median(d), p10_nn_m = quantile(d, 0.1), p90_nn_m = quantile(d, 0.9))
    }

    sv_tab <- xy_d %>% count(loc, year)
    sv_multi <- sv_tab %>% filter(n >= min_survey) %>% count(loc) %>% filter(n >= 2) %>% pull(loc)
    d_survey <- sv_tab %>% filter(loc %in% sv_multi, n >= min_survey) %>%
      pmap_dfr(function(loc, year, n) {
        te <- filter(xy_d, .data$loc == !!loc, .data$year == !!year)
        tr <- filter(xy_d, .data$loc == !!loc, .data$year != !!year)
        dist_row("Survey transfer", sprintf("%s %d", loc, year), "Other surveys, same region", te, tr)
      })

    d_lolo <- map_dfr(valid_locs, function(l) {
      dist_row("Unseen regions (LOLO-CV)", l, "Other five regions",
               filter(xy_d, loc == l), filter(xy_d, loc != l))
    })

    d_transfer <- expand_grid(Source = valid_locs, Test = valid_locs) %>%
      filter(Source != Test) %>%
      pmap_dfr(function(Source, Test) {
        dist_row("Single-region transfer", Test, Source,
                 filter(xy_d, loc == Test), filter(xy_d, loc == Source))
      })

    dist_tab <- bind_rows(d_survey, d_lolo, d_transfer)
    write_csv(dist_tab, out_path("agc_design_distances"))
    cat("Design distances saved\n")
    print(dist_tab %>% group_by(Design) %>%
            summarise(n_cases = n(), median_nn_km = median(median_nn_m) / 1000, .groups = "drop"))
  }
}


# =============================================================================
# PART K: Proxy workflow in the 250 m block design
#   Test blocks and seeds as in PART I (b = 250 m). Every proxy model is
#   refitted without the points that lie in a test block. Point coordinates
#   come from the dataPA files (gee_id); points without coordinates stay in
#   training unless they are test samples or test pixels.
# =============================================================================
if (run_blockproxy) {
  cat("\n=====================================================\n")
  cat(" PART K: Proxy workflow in the 250 m block design\n")
  cat("=====================================================\n")

  b_q  <- 250
  xy_q <- make_xy_df(agc_df)
  if (is.null(xy_q)) {
    cat("No coordinates; PART K skipped\n")
  } else {
    centres_q <- xy_q %>% group_by(loc) %>%
      summarise(lon0 = mean(lon), lat0 = mean(lat), .groups = "drop")
    pa_xy <- map_dfr(pa_files, function(p) {
      d <- read_csv(p, show_col_types = FALSE, col_types = cols(.default = col_character()))
      if (!all(c("gee_id", lon_col, lat_col) %in% names(d))) return(NULL)
      transmute(d, gee_id,
                lon_p = suppressWarnings(as.numeric(.data[[lon_col]])),
                lat_p = suppressWarnings(as.numeric(.data[[lat_col]])))
    }) %>% filter(!is.na(gee_id), !is.na(lon_p), !is.na(lat_p)) %>%
      distinct(gee_id, .keep_all = TRUE)
    cat(sprintf("%d points with coordinates in the dataPA files\n", nrow(pa_xy)))

    block_of <- function(d) {
      d %>% select(gee_id, loc) %>%
        left_join(pa_xy, by = "gee_id") %>%
        left_join(centres_q, by = "loc") %>%
        transmute(block = ifelse(is.na(lon_p) | is.na(lon0), NA_character_,
                                 paste(loc,
                                       floor((lon_p - lon0) * 111320 * cos(lat0 * pi / 180) / b_q),
                                       floor((lat_p - lat0) * 110540 / b_q), sep = "_"))) %>%
        pull(block)
    }

    q_iters   <- if (test_mode) 1 else n_iter_blk
    q_results <- list()
    for (it in seq_len(q_iters)) {
      ckpt <- file.path(ext_ckpt_dir, sprintf("blockproxy_%04dm_iter_%02d.rds", b_q, it))
      if (file.exists(ckpt)) {
        q_results[[it]] <- readRDS(ckpt)
        next
      }
      t0 <- Sys.time()
      res <- tryCatch({
        s <- seed_base + 50000 + b_q + it                  # same as PART I
        d <- xy_q %>% mutate(block = paste(loc, floor(x_m / b_q), floor(y_m / b_q), sep = "_"))
        set.seed(s)
        test_blocks <- d %>%
          distinct(loc, block) %>% arrange(loc, block) %>%
          group_by(loc) %>%
          filter(block %in% sample(block, max(1, round(test_frac * n())))) %>%
          pull(block)
        test_ids <- d$gee_id[d$block %in% test_blocks]
        test_pix <- d$pix[d$block %in% test_blocks]

        keep_fn <- function(x) {
          blk <- block_of(x)
          x[!(blk %in% test_blocks) & !(x$gee_id %in% test_ids) & !(x$pix %in% test_pix), ]
        }
        fold_df <- build_fold_proxies(keep_fn, s) %>%
          inner_join(select(d, gee_id, pix, block), by = c("gee_id", "pix"))
        train_df <- filter(fold_df, !(block %in% test_blocks))
        test_df  <- filter(fold_df,   block %in% test_blocks)

        m_px  <- fit_rf_fixed(train_df, agc_predictors,      agc_response, agc_hp, s + 10)
        m_gse <- fit_rf_fixed(train_df, gse_only_predictors, agc_response, agc_hp, s + 20)
        reg_mean <- train_df %>% group_by(loc) %>%
          summarise(mu = mean(.data[[agc_response]]), .groups = "drop")
        test_df <- test_df %>% left_join(reg_mean, by = "loc")
        test_df$pred_px  <- predict(m_px,  newdata = as.data.frame(test_df[, agc_predictors]))
        test_df$pred_gse <- predict(m_gse, newdata = as.data.frame(test_df[, gse_only_predictors]))

        sc <- function(x, scope) {
          score_models(x[[agc_response]],
                       list("Proxy workflow"      = x$pred_px,
                            "GSE bands only"      = x$pred_gse,
                            "Mean-value baseline" = x$mu)) %>%
            mutate(Scope = scope, n_test = nrow(x), .before = 1)
        }
        bind_rows(sc(test_df, "All regions"),
                  map_dfr(sort(unique(test_df$loc)), function(l) sc(filter(test_df, loc == l), l))) %>%
          mutate(Block_m = b_q, Iteration = it, n_train = nrow(train_df), .before = 1)
      }, error = function(e) {
        message(sprintf("[BlockProxy] iteration %d failed: %s", it, conditionMessage(e)))
        NULL
      })
      if (!is.null(res)) {
        saveRDS(res, ckpt)
        q_results[[it]] <- res
        ov <- filter(res, Scope == "All regions")
        cat(sprintf("Block proxy iteration %d | MEC proxy %.2f | GSE %.2f | %.1f min\n", it,
                    ov$MEC[ov$Model == "Proxy workflow"], ov$MEC[ov$Model == "GSE bands only"],
                    mins_since(t0)))
      }
    }

    q_df <- bind_rows(q_results)
    if (nrow(q_df) > 0) {
      write_csv(q_df, out_path("agc_block_proxy_iterations"))
      paired <- q_df %>% filter(Scope == "All regions") %>%
        select(Iteration, Model, RMSE) %>%
        pivot_wider(names_from = Model, values_from = RMSE) %>%
        mutate(diff = `Proxy workflow` - `GSE bands only`)
      q_summary <- q_df %>%
        group_by(Scope, Model) %>%
        summarise(n_iter = n(), n_test = round(mean(n_test)),
                  RMSE_mean = mean(RMSE), RMSE_lower = q_lo(RMSE), RMSE_upper = q_hi(RMSE),
                  R2_mean = if (all(is.na(R2))) NA_real_ else mean(R2, na.rm = TRUE),
                  MEC_mean = mean(MEC), MEC_lower = q_lo(MEC), MEC_upper = q_hi(MEC),
                  abs_bias_pct = mean(abs(Bias_pct)), .groups = "drop") %>%
        mutate(Model = factor(Model, levels = model_levels)) %>%
        arrange(Scope != "All regions", Scope, Model) %>%
        mutate(Model = as.character(Model),
               paired_RMSE_diff_mean = ifelse(Scope == "All regions", mean(paired$diff), NA_real_),
               paired_proxy_better_n = ifelse(Scope == "All regions", sum(paired$diff < 0), NA_integer_))
      write_csv(q_summary, out_path("agc_block_proxy_summary"))
      cat("Block proxy summary saved\n")
      print(q_summary %>% filter(Scope == "All regions") %>%
              select(Model, n_test, RMSE_mean, R2_mean, MEC_mean, paired_RMSE_diff_mean, paired_proxy_better_n) %>%
              mutate(across(where(is.double), ~ round(.x, 3))), width = Inf)
    }
  }
}


# =============================================================================
# PART L: Cross-fitted proxies (stacking) in the pixel and 250 m block designs
#   Outer split: the test pixels of PART D (pixel design) or the test blocks of
#   PART K (250 m blocks), with the same seeds.
#   Test samples: proxies from proxy models fitted on all outer training data,
#   exactly as in PARTS D and K.
#   Training samples, two versions:
#     in-sample    : the proxies of PARTS D and K (the proxy models saw these
#                    samples when they were fitted)
#     cross-fitted : the outer training units (pixels or 250 m blocks) are split
#                    per region into cf_k groups. The proxies of group j come
#                    from proxy models fitted without group j (and without the
#                    outer test data).
#   Models scored on the same test samples:
#     Proxy workflow                (in-sample proxies, reproduces PART D or K)
#     Proxy workflow, cross-fitted  (same AGC settings and seed)
#     GSE bands only, Mean-value baseline
#   Diagnostics: R2 between each proxy and the AGC target in the training
#   samples (in-sample and cross-fitted) and in the test samples, and the share
#   of the RF importance (IncNodePurity) held by the proxy layers.
# =============================================================================
model_levels_cf <- c("Proxy workflow", "Proxy workflow, cross-fitted",
                     "GSE bands only", "Mean-value baseline")
proxy_cols <- c(agc_pred_vars, morph_pred)

# Corrected resampled t-test (Nadeau & Bengio, 2003) for paired differences
corrected_t <- function(d, n_train, n_test) {
  d <- d[!is.na(d)]
  J <- length(d)
  if (J < 2 || var(d) == 0) return(c(t = NA_real_, p = NA_real_))
  t <- mean(d) / sqrt((1 / J + n_test / n_train) * var(d))
  c(t = t, p = 2 * pt(-abs(t), df = J - 1))
}

# Share of IncNodePurity held by the proxy layers
proxy_share <- function(m) {
  imp <- randomForest::importance(m)[, "IncNodePurity"]
  sum(imp[names(imp) %in% proxy_cols]) / sum(imp)
}

# R2 between each proxy and the AGC target in one set of rows
proxy_r2 <- function(d, set_name) {
  map_dfr(c("tAGB_pred", "carbon_index", "tSPC_pred", "PA_prob"), function(v) {
    ok <- !is.na(d[[v]]) & !is.na(d[[agc_response]])
    r2 <- if (sum(ok) > 2 && sd(d[[v]][ok]) > 0) cor(d[[v]][ok], d[[agc_response]][ok])^2 else NA_real_
    tibble(Set = set_name, Proxy = v, R2_with_AGC = r2, n = sum(ok))
  })
}

if (run_crossfit) {
  cat("\n=====================================================\n")
  cat(sprintf(" PART L: Cross-fitted proxies (cf_k = %d)\n", cf_k))
  cat("=====================================================\n")

  # L1. Block helpers (as in PART K), rebuilt so that PART L runs on its own
  b_r  <- 250
  xy_r <- make_xy_df(agc_df)
  if (is.null(xy_r) && "Block 250 m" %in% cf_designs) {
    cat("No coordinates; the 250 m block design of PART L is skipped\n")
    cf_designs <- setdiff(cf_designs, "Block 250 m")
  }
  if (!is.null(xy_r)) {
    centres_r <- xy_r %>% group_by(loc) %>%
      summarise(lon0 = mean(lon), lat0 = mean(lat), .groups = "drop")
    pa_xy_r <- map_dfr(pa_files, function(p) {
      d <- read_csv(p, show_col_types = FALSE, col_types = cols(.default = col_character()))
      if (!all(c("gee_id", lon_col, lat_col) %in% names(d))) return(NULL)
      transmute(d, gee_id,
                lon_p = suppressWarnings(as.numeric(.data[[lon_col]])),
                lat_p = suppressWarnings(as.numeric(.data[[lat_col]])))
    }) %>% filter(!is.na(gee_id), !is.na(lon_p), !is.na(lat_p)) %>%
      distinct(gee_id, .keep_all = TRUE)
    block_of_r <- function(d) {
      d %>% select(gee_id, loc) %>%
        left_join(pa_xy_r, by = "gee_id") %>%
        left_join(centres_r, by = "loc") %>%
        transmute(block = ifelse(is.na(lon_p) | is.na(lon0), NA_character_,
                                 paste(loc,
                                       floor((lon_p - lon0) * 111320 * cos(lat0 * pi / 180) / b_r),
                                       floor((lat_p - lat0) * 110540 / b_r), sep = "_"))) %>%
        pull(block)
    }
  }

  design_key <- function(design) {
    if (design == "Pixel-grouped") "pixel" else sprintf("block%04dm", b_r)
  }

  # L2. Outer split of one design and iteration (same seeds as PART D / K).
  #     keep()       : rows each proxy model may use (outer test data removed)
  #     drop_units() : removes one inner group from proxy training data
  #     prep()       : adds the unit used for the inner groups
  #     is_test()    : marks the outer test rows of the AGC data
  outer_def <- function(design, it) {
    if (design == "Pixel-grouped") {
      s <- seed_base + 9000 + 100 * it                     # as PART D
      test_pix <- draw_test_pix_sorted(agc_df, s)
      test_ids <- agc_df$gee_id[agc_df$pix %in% test_pix]
      list(seed       = s,
           keep       = function(x) x[!(x$pix %in% test_pix) & !(x$gee_id %in% test_ids), ],
           drop_units = function(x, units, ids, pixs) x[!(x$pix %in% pixs) & !(x$gee_id %in% ids), ],
           prep       = function(f) mutate(f, unit = pix),
           is_test    = function(f) f$pix %in% test_pix)
    } else {
      s <- seed_base + 50000 + b_r + it                    # as PARTS I and K
      d <- xy_r %>% mutate(block = paste(loc, floor(x_m / b_r), floor(y_m / b_r), sep = "_"))
      set.seed(s)
      test_blocks <- d %>%
        distinct(loc, block) %>% arrange(loc, block) %>%
        group_by(loc) %>%
        filter(block %in% sample(block, max(1, round(test_frac * n())))) %>%
        pull(block)
      test_ids <- d$gee_id[d$block %in% test_blocks]
      test_pix <- d$pix[d$block %in% test_blocks]
      list(seed = s,
           keep = function(x) {
             blk <- block_of_r(x)
             x[!(blk %in% test_blocks) & !(x$gee_id %in% test_ids) & !(x$pix %in% test_pix), ]
           },
           drop_units = function(x, units, ids, pixs) {
             blk <- block_of_r(x)
             x[!(blk %in% units) & !(x$gee_id %in% ids) & !(x$pix %in% pixs), ]
           },
           prep    = function(f) f %>%
             inner_join(select(d, gee_id, pix, block), by = c("gee_id", "pix")) %>%
             mutate(unit = block),
           is_test = function(f) f$block %in% test_blocks)
    }
  }

  # L3. Cross-fitted proxies for the outer training rows
  crossfit_train <- function(od, train, key) {
    set.seed(od$seed + 777)
    grp <- train %>%
      distinct(loc, unit) %>% arrange(loc, unit) %>%
      group_by(loc) %>%
      mutate(grp = rep_len(seq_len(cf_k), n())[sample.int(n())]) %>%
      ungroup()
    train <- left_join(train, grp, by = c("loc", "unit"))

    map_dfr(seq_len(cf_k), function(j) {
      ck <- file.path(ext_ckpt_dir, sprintf("%s_inner_%d.rds", key, j))
      if (file.exists(ck)) return(readRDS(ck))
      g_rows <- filter(train, grp == j)
      if (nrow(g_rows) == 0) return(NULL)
      units <- unique(g_rows$unit)
      ids   <- unique(g_rows$gee_id)
      pixs  <- unique(g_rows$pix)
      keep_j <- function(x) od$drop_units(od$keep(x), units, ids, pixs)
      cat(sprintf("   inner group %d/%d: %d AGC rows\n", j, cf_k, nrow(g_rows)))
      out <- build_fold_proxies(keep_j, od$seed + 1000 * j,
                                target_df = semi_join(agc_df, g_rows, by = c("gee_id", "pix"))) %>%
        inner_join(select(g_rows, gee_id, pix, unit, grp), by = c("gee_id", "pix"))
      saveRDS(out, ck)
      out
    })
  }

  # L4. Run
  r_iters   <- if (test_mode) 1 else n_iter_mc
  r_results <- list()
  for (design in cf_designs) {
    for (it in seq_len(r_iters)) {
      key  <- sprintf("crossfit_%s_k%d_iter_%02d", design_key(design), cf_k, it)
      ckpt <- file.path(ext_ckpt_dir, paste0(key, ".rds"))
      if (file.exists(ckpt)) {
        r_results[[key]] <- readRDS(ckpt)
        next
      }
      cat(sprintf("\n[CrossFit] %s, iteration %d/%d\n", design, it, r_iters))
      t0 <- Sys.time()

      res <- tryCatch({
        od <- outer_def(design, it)
        s  <- od$seed

        ck_outer <- file.path(ext_ckpt_dir, paste0(key, "_outer.rds"))
        fold <- if (file.exists(ck_outer)) readRDS(ck_outer) else {
          f <- od$prep(build_fold_proxies(od$keep, s))
          saveRDS(f, ck_outer)
          f
        }
        train <- fold[!od$is_test(fold), ]
        test  <- fold[ od$is_test(fold), ]
        if (nrow(test) == 0) stop("no test rows after filtering")

        cf_train <- crossfit_train(od, train, key)

        m_std <- fit_rf_fixed(train,    agc_predictors,      agc_response, agc_hp, s + 10)
        m_cf  <- fit_rf_fixed(cf_train, agc_predictors,      agc_response, agc_hp, s + 10)
        m_gse <- fit_rf_fixed(train,    gse_only_predictors, agc_response, agc_hp, s + 20)

        reg_mean <- train %>% group_by(loc) %>%
          summarise(mu = mean(.data[[agc_response]]), .groups = "drop")
        te <- test %>% left_join(reg_mean, by = "loc")
        te$pred_std <- predict(m_std, newdata = as.data.frame(te[, agc_predictors]))
        te$pred_cf  <- predict(m_cf,  newdata = as.data.frame(te[, agc_predictors]))
        te$pred_gse <- predict(m_gse, newdata = as.data.frame(te[, gse_only_predictors]))

        share_std <- proxy_share(m_std)
        share_cf  <- proxy_share(m_cf)

        sc <- function(x, scope) {
          score_models(x[[agc_response]],
                       list("Proxy workflow"               = x$pred_std,
                            "Proxy workflow, cross-fitted" = x$pred_cf,
                            "GSE bands only"               = x$pred_gse,
                            "Mean-value baseline"          = x$mu)) %>%
            mutate(Scope = scope, n_test = nrow(x), .before = 1)
        }
        scores <- bind_rows(sc(te, "All regions"),
                            map_dfr(sort(unique(te$loc)), function(l) sc(filter(te, loc == l), l))) %>%
          mutate(share_proxy_imp = case_when(Model == "Proxy workflow" ~ share_std,
                                             Model == "Proxy workflow, cross-fitted" ~ share_cf,
                                             TRUE ~ NA_real_)) %>%
          mutate(Design = design, Iteration = it,
                 n_train = nrow(train), n_train_cf = nrow(cf_train), .before = 1)

        diag <- bind_rows(proxy_r2(train,    "Training, in-sample proxies"),
                          proxy_r2(cf_train, "Training, cross-fitted proxies"),
                          proxy_r2(test,     "Test")) %>%
          mutate(Design = design, Iteration = it, .before = 1)

        list(scores = scores, diag = diag)
      }, error = function(e) {
        message(sprintf("[CrossFit] %s failed: %s", key, conditionMessage(e)))
        NULL
      })

      if (!is.null(res)) {
        saveRDS(res, ckpt)
        r_results[[key]] <- res
        # inner and outer checkpoints are no longer needed
        unlink(file.path(ext_ckpt_dir, c(paste0(key, "_outer.rds"),
                                         sprintf("%s_inner_%d.rds", key, seq_len(cf_k)))))
        ov <- filter(res$scores, Scope == "All regions")
        cat(sprintf("%s, iteration %d | MEC proxy %.2f | cross-fitted %.2f | GSE %.2f | proxy importance %.0f%% -> %.0f%% | %.1f min\n",
                    design, it,
                    ov$MEC[ov$Model == "Proxy workflow"],
                    ov$MEC[ov$Model == "Proxy workflow, cross-fitted"],
                    ov$MEC[ov$Model == "GSE bands only"],
                    100 * ov$share_proxy_imp[ov$Model == "Proxy workflow"],
                    100 * ov$share_proxy_imp[ov$Model == "Proxy workflow, cross-fitted"],
                    mins_since(t0)))
      }
    }
  }

  # L5. Outputs
  r_scores <- bind_rows(map(r_results, "scores"))
  r_diag   <- bind_rows(map(r_results, "diag"))

  if (nrow(r_scores) > 0) {
    write_csv(r_scores, out_path("agc_crossfit_iterations"))

    r_summary <- r_scores %>%
      group_by(Design, Scope, Model) %>%
      summarise(n_iter = n(), n_test = round(mean(n_test)),
                RMSE_mean = mean(RMSE), RMSE_lower = q_lo(RMSE), RMSE_upper = q_hi(RMSE),
                R2_mean = if (all(is.na(R2))) NA_real_ else mean(R2, na.rm = TRUE),
                MEC_mean = mean(MEC, na.rm = TRUE), MEC_lower = q_lo(MEC), MEC_upper = q_hi(MEC),
                abs_bias_pct = mean(abs(Bias_pct)),
                share_proxy_imp = if (all(is.na(share_proxy_imp))) NA_real_ else mean(share_proxy_imp, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(Model = factor(Model, levels = model_levels_cf)) %>%
      arrange(Design, Scope != "All regions", Scope, Model) %>%
      mutate(Model = as.character(Model))
    write_csv(r_summary, out_path("agc_crossfit_summary"))

    # Paired differences on the same test samples (all regions), A minus B
    pairs_cf <- list(c("Proxy workflow, cross-fitted", "GSE bands only"),
                     c("Proxy workflow, cross-fitted", "Proxy workflow"),
                     c("Proxy workflow",               "GSE bands only"))
    r_paired <- map_dfr(unique(r_scores$Design), function(dz) {
      w <- filter(r_scores, Design == dz, Scope == "All regions")
      map_dfr(pairs_cf, function(pr) {
        a <- filter(w, Model == pr[1]) %>% arrange(Iteration)
        b <- filter(w, Model == pr[2]) %>% arrange(Iteration)
        common <- intersect(a$Iteration, b$Iteration)
        a <- a[a$Iteration %in% common, ]
        b <- b[b$Iteration %in% common, ]
        dR <- a$RMSE - b$RMSE
        ct <- corrected_t(dR, mean(a$n_train), mean(a$n_test))
        tibble(Design = dz, Model_A = pr[1], Model_B = pr[2], n_iter = length(common),
               dRMSE_mean = mean(dR), dMEC_mean = mean(a$MEC - b$MEC),
               A_better_n = sum(dR < 0),
               t_corrected = unname(ct["t"]), p_corrected = unname(ct["p"]))
      })
    })
    write_csv(r_paired, out_path("agc_crossfit_paired"))

    r_diag_sum <- r_diag %>%
      group_by(Design, Set, Proxy) %>%
      summarise(n_iter = n(), n = round(mean(n)),
                R2_with_AGC = mean(R2_with_AGC, na.rm = TRUE), .groups = "drop") %>%
      arrange(Design, Proxy, Set)
    write_csv(r_diag_sum, out_path("agc_crossfit_proxy_r2"))

    cat("Cross-fit results saved\n")
    print(r_summary %>% filter(Scope == "All regions") %>%
            select(Design, Model, n_iter, RMSE_mean, R2_mean, MEC_mean, share_proxy_imp) %>%
            mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)
    print(r_paired %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)
    print(r_diag_sum %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)
  }
}


cat(sprintf("\nExtensions finished. Total time: %.1f min\n", mins_since(t_ext)))
if (test_mode) cat("Test run OK? Set ext_test_mode <- FALSE and run again.\n")
