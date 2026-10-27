# =============================================================================
# 01_proxy_models_and_pixel_sharing.R
# Accuracy of the presence, leaf morphology and carbon index models with two
# kinds of test sets (Table S1 of the paper), and the share of AGC quadrats that
# share a 10 m pixel with another quadrat (Table S5).
#
#   PART A  Pixel sharing. A pixel is a unique set of 64 GSE values in one year.
#   PART B  Monte Carlo with ten repeats and 30% of the samples held out:
#             random rows     single samples are held out
#             pixel-grouped   all samples of a pixel are held out together,
#                             stratified by region
#           Percent cover and AGB are tested in 02_agc_pixel_lolo_transfer.R
#           (PART E), with the final targets.
#   PART C  Summary table and figure
#
# The models use the settings of the grid search (grid_rf_best_*.csv), with no
# new tuning. Every repeat is saved as a checkpoint, so a stopped run can be
# started again and finished repeats are skipped.
#
# Inputs (result/)
#   training_data_for_GEE_PA.csv, training_data_for_GEE_MORPH3_probs.csv,
#   training_data_for_GEE_CINDEX.csv, training_data_for_GEE_AGC.csv,
#   grid_rf_best_{PA,MORPH3,CINDEX}.csv, dataPA_<year>.csv (location of each point)
# Outputs (result/)
#   agc_pixel_sharing_fullPipeline.csv           Table S5
#   proxy_stage_MC_iterations_fullPipeline.csv   one row per repeat and metric
#   proxy_stage_MC_summary_fullPipeline.csv      Table S1 (presence, morphology,
#                                                carbon index)
#   proxy_stage_MC_fullPipeline.png
#
# How to run: keep test_mode <- TRUE for a quick check (outputs get "_TEST"),
# then set test_mode <- FALSE for the full run.
# =============================================================================


# -----------------------------------------------------------------------------
# STEP 0: Libraries, paths and run settings
# -----------------------------------------------------------------------------
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(randomForest)
library(Metrics)
library(caret)
library(ggplot2)

# Paths: open seagrass-agc-proxy-indonesia.Rproj in RStudio (or start R in the
# repository root), so that here() returns the repository folder.
library(here)

result_dir <- here("result")   # input tables and outputs (not distributed, see README)

# TRUE  = quick check (1 repeat, 50 trees)
# FALSE = full run
test_mode <- TRUE

min_samples  <- 10       # minimum rows for a model
seed_base    <- 42
n_iter_mc    <- 10       # Monte Carlo repeats
test_frac    <- 0.3      # share of pixels (or rows) used for testing
pa_threshold <- 0.6      # PA probability threshold used for model evaluation

suffix   <- if (test_mode) "_fullPipeline_TEST" else "_fullPipeline"
ckpt_dir <- file.path(result_dir, paste0("checkpoint", suffix))
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

t_start <- Sys.time()
cat(sprintf("Proxy model validation started (%s mode)\n",
            if (test_mode) "TEST" else "FULL"))


# -----------------------------------------------------------------------------
# STEP 1: Predictor definitions
# -----------------------------------------------------------------------------
env_vars  <- c("depth", "distToLand")
gse_vars  <- paste0("GSE_A", sprintf("%02d", 0:63))
base_vars <- c(env_vars, gse_vars)            # PA, MORPH3 and CINDEX models

morph_classes <- c("mixed_short_plus_mono_short", "mixed_long", "mono_Ea")


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
print(count(loc_lookup, loc))


# -----------------------------------------------------------------------------
# STEP 3: Load the training rows of each model
# -----------------------------------------------------------------------------
# Pixel key: samples in the same 10 m pixel and year have identical GSE values.
# Values are rounded to 5 decimals to avoid tiny differences from CSV export.
make_pix <- function(df) {
  m <- round(as.matrix(df[, gse_vars]), 5)
  do.call(paste, c(as.data.frame(m), sep = "|"))
}

load_stage <- function(file, response, predictors, stage) {
  path <- file.path(result_dir, file)
  if (!file.exists(path)) stop(sprintf("[%s] File not found: %s", stage, path))

  df <- read_csv(path, show_col_types = FALSE, guess_max = 100000)

  needed  <- c("gee_id", response, predictors)
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) {
    stop(sprintf("[%s] Missing columns: %s", stage, paste(missing, collapse = ", ")))
  }

  # Keep only the needed columns, then attach the location of each point
  df <- df %>%
    select(all_of(needed)) %>%
    mutate(gee_id = as.character(gee_id)) %>%
    left_join(loc_lookup, by = "gee_id")

  n_noloc <- sum(is.na(df$loc))
  if (n_noloc > 0) {
    warning(sprintf("[%s] %d rows without a location were removed.", stage, n_noloc))
  }

  df <- df %>%
    filter(!is.na(loc)) %>%
    filter(if_all(all_of(c(response, predictors)), ~ !is.na(.x)))

  df$pix <- make_pix(df)

  cat(sprintf("[%s] %d training rows loaded\n", stage, nrow(df)))
  df
}

# 3a. Presence-absence (classification, levels "0" and "1")
pa_df <- load_stage("training_data_for_GEE_PA.csv", "PA", base_vars, "PA") %>%
  mutate(PA = factor(as.character(PA), levels = c("0", "1")))

# 3b. Leaf morphology (classification, three classes)
morph_df <- load_stage("training_data_for_GEE_MORPH3_probs.csv", "morph3",
                       base_vars, "MORPH3")
unknown_class <- setdiff(unique(as.character(morph_df$morph3)), morph_classes)
if (length(unknown_class) > 0) {
  stop("[MORPH3] Unexpected class names: ", paste(unknown_class, collapse = ", "))
}
morph_df$morph3 <- factor(as.character(morph_df$morph3), levels = morph_classes)

# 3c. Carbon index
ci_df <- load_stage("training_data_for_GEE_CINDEX.csv", "carbon_index",
                    base_vars, "CINDEX")

# 3d. AGC rows (for the pixel sharing in PART A). Only the target, location
#     and GSE columns are needed.
agc_path <- file.path(result_dir, "training_data_for_GEE_AGC.csv")
if (!file.exists(agc_path)) stop("[AGC] File not found: ", agc_path)

agc_df <- read_csv(agc_path, show_col_types = FALSE, guess_max = 100000)

missing_agc <- setdiff(c("gee_id", "loc", "AGC_pred", base_vars), names(agc_df))
if (length(missing_agc) > 0) {
  stop("[AGC] Missing columns: ", paste(missing_agc, collapse = ", "))
}

agc_df <- agc_df %>%
  filter(!is.na(loc), !is.na(AGC_pred), AGC_pred > 0) %>%
  mutate(gee_id = as.character(gee_id)) %>%
  select(gee_id, loc, AGC_pred, all_of(base_vars)) %>%
  filter(if_all(all_of(base_vars), ~ !is.na(.x)))

agc_df$pix <- make_pix(agc_df)

cat(sprintf("[AGC] %d rows loaded\n", nrow(agc_df)))
print(count(agc_df, loc))


# -----------------------------------------------------------------------------
# STEP 4: Settings of the models (grid search, no new tuning)
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
           CINDEX = read_best("CINDEX"))

cat("\nModel settings (fixed):\n")
print(bind_rows(hp, .id = "stage"))


# -----------------------------------------------------------------------------
# STEP 5: Helper functions
# -----------------------------------------------------------------------------

# 5a. Fit one proxy model with fixed hyperparameters.
#     For classification with only one class in the fold, a constant
#     probability is returned instead of a Random Forest.
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
#     classes = NULL  -> regression, returns a numeric vector
#     classes = c(..) -> classification, returns a probability matrix
predict_stage <- function(fit, newdata, predictors, classes = NULL) {
  X  <- as.data.frame(newdata[, predictors])
  ok <- complete.cases(X)

  if (is.null(classes)) {
    out <- rep(NA_real_, nrow(X))
    if (any(ok)) out[ok] <- predict(fit$model, newdata = X[ok, , drop = FALSE])
    return(out)
  }

  out <- matrix(NA_real_, nrow(X), length(classes),
                dimnames = list(NULL, classes))
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

# 5c. Accuracy metrics (R2 = squared correlation, as in caret::R2)
calc_metrics <- function(obs, pred) {
  tibble(RMSE = Metrics::rmse(obs, pred),
         MAE  = Metrics::mae(obs, pred),
         R2   = caret::R2(pred, obs))
}

# 5d. Elapsed time in minutes
mins_since <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

# 5e. Pixel-grouped test set, stratified by region. Pixels are sorted first,
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

# 5f. 2.5th and 97.5th percentiles
q_lo <- function(x) quantile(x, 0.025, na.rm = TRUE)
q_hi <- function(x) quantile(x, 0.975, na.rm = TRUE)


# =============================================================================
# PART A: AGC quadrats that share a 10 m pixel (Table S5)
# =============================================================================
pix_sizes <- count(agc_df, loc, pix, name = "n_quadrats")

pix_summary <- bind_rows(
  pix_sizes %>% mutate(loc = "All regions"),
  pix_sizes
) %>%
  group_by(Region = loc) %>%
  # summarise() works line by line, so the total is computed last
  summarise(
    n_pixels               = n(),
    max_quadrats_per_pixel = max(n_quadrats),
    pct_quadrats_shared    = round(100 * sum(n_quadrats[n_quadrats > 1]) / sum(n_quadrats), 1),
    n_quadrats             = sum(n_quadrats),
    .groups = "drop"
  ) %>%
  relocate(n_quadrats, .after = Region) %>%
  arrange(Region != "All regions", Region)

write_csv(pix_summary, file.path(result_dir, paste0("agc_pixel_sharing", suffix, ".csv")))
cat("\nAGC quadrats per pixel:\n")
print(pix_summary)


# =============================================================================
# PART B: Random-row vs pixel-grouped Monte Carlo (Table S1)
# =============================================================================
mc_iters <- if (test_mode) 1 else n_iter_mc

proxy_stages <- list(
  PA     = list(df = pa_df,    response = "PA",            predictors = base_vars, hp = hp$PA,     type = "binary"),
  MORPH3 = list(df = morph_df, response = "morph3",        predictors = base_vars, hp = hp$MORPH3, type = "multiclass"),
  CINDEX = list(df = ci_df,    response = "carbon_index",  predictors = base_vars, hp = hp$CINDEX, type = "regression")
)

# Metrics for one stage, returned in long format (Metric, Value)
stage_metrics <- function(type, fit, test, st) {
  if (type == "binary") {
    prob <- predict_stage(fit, test, st$predictors, classes = c("0", "1"))[, "1"]
    pred <- ifelse(prob >= pa_threshold, "1", "0")
    obs  <- as.character(test[[st$response]])
    tp <- sum(pred == "1" & obs == "1"); fp <- sum(pred == "1" & obs == "0")
    fn <- sum(pred == "0" & obs == "1")
    precision <- tp / max(1, tp + fp)
    recall    <- tp / max(1, tp + fn)
    tibble(Metric = c("Accuracy", "Precision", "Recall", "F1"),
           Value  = c(mean(pred == obs), precision, recall,
                      ifelse(precision + recall > 0,
                             2 * precision * recall / (precision + recall), NA)))
  } else if (type == "multiclass") {
    pm   <- predict_stage(fit, test, st$predictors, classes = morph_classes)
    pred <- factor(morph_classes[max.col(pm, ties.method = "first")], levels = morph_classes)
    obs  <- factor(as.character(test[[st$response]]), levels = morph_classes)
    tab  <- table(pred, obs)
    po   <- sum(diag(tab)) / sum(tab)
    pe   <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
    tibble(Metric = c("Accuracy", "Kappa"),
           Value  = c(po, ifelse(pe < 1, (po - pe) / (1 - pe), NA)))
  } else {
    pred <- predict_stage(fit, test, st$predictors)
    m    <- suppressWarnings(calc_metrics(test[[st$response]], pred))
    tibble(Metric = c("RMSE", "MAE", "R2"), Value = c(m$RMSE, m$MAE, m$R2))
  }
}

for (stage in names(proxy_stages)) {
  st <- proxy_stages[[stage]]

  for (mode in c("Random rows", "Pixel-grouped")) {
    cat(sprintf("\n[MC] %s | %s\n", stage, mode))

    for (it in seq_len(mc_iters)) {
      tag  <- if (mode == "Random rows") "random" else "pixel"
      ckpt <- file.path(ckpt_dir, sprintf("proxyMC_%s_%s_iter_%02d.rds", stage, tag, it))
      if (file.exists(ckpt)) {
        cat(sprintf("[MC] %s %s iteration %d loaded from checkpoint\n", stage, tag, it))
        next
      }

      t0      <- Sys.time()
      it_seed <- seed_base + 7000 + 100 * it

      res <- tryCatch({
        d <- st$df

        if (mode == "Random rows") {
          # Random 70/30 split of rows
          set.seed(it_seed)
          tr    <- sample(nrow(d), floor((1 - test_frac) * nrow(d)))
          train <- d[tr, ]
          test  <- d[-tr, ]
        } else {
          test_pix <- draw_test_pix_sorted(d, it_seed)
          train    <- filter(d, !(pix %in% test_pix))
          test     <- filter(d,   pix %in% test_pix)
        }

        fit <- fit_stage(train, st$response, st$predictors, st$hp, it_seed + 1, stage)
        stage_metrics(st$type, fit, test, st) %>%
          mutate(Stage = stage, Mode = mode, Iteration = it,
                 n_train = nrow(train), n_test = nrow(test))
      }, error = function(e) {
        message(sprintf("[MC] %s %s iteration %d failed: %s", stage, tag, it, conditionMessage(e)))
        NULL
      })

      if (!is.null(res)) {
        saveRDS(res, ckpt)
        cat(sprintf("[MC] %s %s iteration %d | %s | %.1f min\n", stage, tag, it,
                    paste(sprintf("%s %.3f", res$Metric, res$Value), collapse = " | "),
                    mins_since(t0)))
      }
    }
  }
}


# =============================================================================
# PART C: Summary table and figure
# =============================================================================
proxy_files <- list.files(ckpt_dir, pattern = "^proxyMC_.*\\.rds$", full.names = TRUE)

if (length(proxy_files) > 0) {
  proxy_df <- map_dfr(proxy_files, readRDS)
  write_csv(proxy_df, file.path(result_dir, paste0("proxy_stage_MC_iterations", suffix, ".csv")))

  proxy_summary <- proxy_df %>%
    group_by(Stage, Mode, Metric) %>%
    summarise(n_iter = n(),
              n_test = round(mean(n_test)),
              mean   = mean(Value, na.rm = TRUE),
              lower  = q_lo(Value),
              upper  = q_hi(Value),
              .groups = "drop") %>%
    mutate(Stage = factor(Stage, levels = c("PA", "MORPH3", "CINDEX"))) %>%
    arrange(Stage, Metric, Mode)

  out_px <- file.path(result_dir, paste0("proxy_stage_MC_summary", suffix, ".csv"))
  write_csv(proxy_summary, out_px)
  cat(sprintf("\nProxy-stage summary saved: %s\n", out_px))
  print(proxy_summary %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)

  # Figure: main metric per stage, random rows vs pixel-grouped
  main_metric <- c(PA = "F1", MORPH3 = "Kappa", CINDEX = "R2")

  p_px <- proxy_summary %>%
    filter(Metric == main_metric[as.character(Stage)]) %>%
    arrange(Stage) %>%
    mutate(Panel = paste0(Stage, " (", Metric, ")"),
           Panel = factor(Panel, levels = unique(Panel)),
           Mode  = factor(Mode, levels = c("Random rows", "Pixel-grouped"))) %>%
    ggplot(aes(x = Mode, y = mean, colour = Mode)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.5) +
    facet_wrap(~ Panel, nrow = 1) +
    scale_colour_manual(values = c("Random rows" = "grey55", "Pixel-grouped" = "black")) +
    coord_cartesian(ylim = c(min(0, proxy_summary$lower[proxy_summary$Metric %in% main_metric], na.rm = TRUE), 1)) +
    labs(title = "Proxy stages: random-row vs pixel-grouped Monte Carlo",
         x = NULL, y = "Test-set score", colour = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", axis.text.x = element_blank())

  print(p_px)
  ggsave(file.path(result_dir, paste0("proxy_stage_MC", suffix, ".png")),
         p_px, width = 10, height = 4.5, dpi = 300)
}


cat(sprintf("\nFinished. Total time: %.1f min\n", mins_since(t_start)))
if (test_mode) {
  cat("Test run OK? Set test_mode <- FALSE and run the full script.\n")
}
