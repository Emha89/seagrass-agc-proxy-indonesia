# Seagrass above-ground carbon from satellite embeddings: proxy workflow and leak-free validation

R, Python and Google Earth Engine (GEE) code for the paper:

Hafizt, M., Salsabila, H.N., Lyons, M., Adi, N.S., Wicaksono, P., Hernawan, U.E.,
Phinn, S., McMahon, K., & Roelfsema, C. **A proxy-based workflow for seagrass
above-ground carbon from satellite embeddings: accuracy declines with distance
from field data.** (Citation details will be added after publication.)

The repository holds code only. No field data, model outputs or GEE assets are
included (see [Data availability](#data-availability)).

## Overview

- **Proxy workflow.** Five Random Forest models predict seagrass presence, leaf
  morphology class probabilities, percent cover, above-ground biomass (AGB) and
  a species-weighted carbon index from the 64 bands of the Google Satellite
  Embeddings (GSE), water depth and distance to land. A final model predicts
  above-ground carbon (AGC) from the GSE bands and these proxy layers. The AGB
  and AGC models are trained on the field point estimates.
- **Leak-free validation.** In every repeat, fold and pair, all five proxy
  models are refitted with the training samples only, and all quadrats of one
  10 m pixel stay on the same side of a split. The designs are the pixel
  hold-out, spatial blocks of 100 to 1000 m, a new survey of the same region,
  unseen regions (leave-one-location-out) and single-region transfer. Each is
  compared with a GSE-only model and a mean-value baseline.
- **Accuracy and distance.** The distance from each test sample to the nearest
  training sample is recorded for every design, and the distance from each
  mapped pixel to the nearest AGC sample is computed in GEE.
- **Annual maps.** The models are refitted in GEE with the same settings and
  applied to the annual GSE layers of 2017 to 2024.

## Repository structure

```
seagrass-agc-proxy-indonesia/
├── R/
│   ├── functions/        helper functions sourced by the workflow scripts
│   ├── workflow/         01-18  data preparation and models, in run order
│   ├── validation/       01-03  leak-free validation designs
│   └── figures_tables/   figures and supplementary tables of the paper
├── gee/                  01-12  Google Earth Engine scripts, in run order
├── seagrass-agc-proxy-indonesia.Rproj
├── LICENSE
└── README.md
```

Folders that are created when the code runs and are not part of the
repository:

```
data/      your input files (see below)
result/    all outputs of the R scripts
```

## Requirements

- R with the packages `here`, `dplyr`, `readr`, `purrr`, `tidyr`, `tibble`,
  `stringr`, `randomForest`, `ranger`, `caret`, `Metrics`, `ggplot2`, `ggpubr`,
  `sf`, `terra`, `patchwork`, `ragg`, `scales` and `svglite`
- Python 3 with `pandas` and `matplotlib`
- A Google Earth Engine account

Open `seagrass-agc-proxy-indonesia.Rproj` in RStudio, or start R in the
repository root. All R scripts find the `data/` and `result/` folders with
`here()`, so no paths need to be changed.

## Input data

Put these files in `data/`:

| File | Content |
|---|---|
| `field_points.csv` | One row per field point: `gee_id`, `OBJECTID`, `year`, `loc` (region), `xcoord`, `ycoord` (decimal degrees), `PA` (0/1), `tSPC` (percent cover), `compositio` (key to the species table), and the point estimate and 95% bounds of AGB and AGC (`AGB_pred`, `AGB_low`, `AGB_up`, `AGC_pred`, `AGC_low`, `AGC_up`) |
| `gse_points_<year>.csv` | GSE bands (`GSE_A00` to `GSE_A63`), `depth` and `distToLand` at the field points, exported by `gee/01_extract_gse_at_field_points.js` |
| `study_regions.shp` | Outline of the study regions, with the region name in `loc` |
| `species_composition.csv` | Percent cover of each species (`Ea_SPC`, `Th_SPC`, ...) and the meadow type (`sg_morpho`) for each `compositio` key |

In GEE, set `ASSET_ROOT` at the top of every script to your own asset folder,
and upload there: `field_points`, `study_regions`, `bathymetry` (Allen Coral
Atlas bathymetry) and `wave_2024` (wave and elevation layer, sampled by
`gee/01` but not used in the models). Two more point tables are uploaded during
the run (see below). The distance-to-land layer
(`projects/coral_atlas/global_datasets/osm_distToLand_indo`) is a shared Allen
Coral Atlas asset that needs its own access.

## How to run

### 1. GSE at the field points (GEE)

Run `gee/01_extract_gse_at_field_points.js` and copy the exported
`gse_points_<year>.csv` files to `data/`.

### 2. Data preparation and models (R)

Run the scripts in `R/workflow/` in number order. Two scripts run twice, because
the presence tables need the outputs of later steps:

| Step | Script | Main outputs (in `result/`) |
|---|---|---|
| 01 | `01_prepare_point_data.R` | `dataPA_<year>.csv` |
| 02 | `02_presence_model.R` | presence model, `grid_rf_best_PA.csv` |
| 03 | `03_presence_apply.R` | `predicted_PA_<year>.csv` |
| 04 | `04_presence_thresholds.R` | accuracy and best F1 threshold per year |
| 05 | `05_species_composition.R` | `dataSP_<year>.csv` |
| 06 | `06_carbon_index_field.R` | field carbon index, `dataCARBON_<year>.csv` |
| 07 | `07_morphology_model.R` | morphology model, `training_MORPH3_points_ALL.csv` |
| 08 | `08_morphology_apply.R` | `predicted_MORPH3_probs_<year>.csv` |
| | **Run 01 and 03 again** | adds the carbon index and the morphology probabilities |
| 09 | `09_cover_model.R` | percent cover model |
| 10 | `10_cover_apply.R` | `predicted_tSPC_<year>.csv` |
| 11 | `11_carbon_index_model.R` | carbon index model |
| 12 | `12_carbon_index_apply.R` | `predicted_indexC_<year>.csv` |
| 13 | `13_agb_data.R` | `dataAGB_<year>.csv` (also sourced by 14 and 15) |
| 14 | `14_agb_tuning.R` | `grid_rf_best_AGB.csv` |
| 15 | `15_agb_model.R` | final AGB model, `predicted_tAGB_<year>.csv` |
| 16 | `16_agc_data.R` | `dataAGC_<year>.csv` (also sourced by 17 and 18) |
| 17 | `17_agc_tuning.R` | `grid_rf_best_AGC.csv` |
| 18 | `18_agc_model.R` | final AGC model with the predicted carbon index |

The Random Forest settings of every model come from a grid search over five
random 70/30 splits (Table S3 of the paper).

### 3. Annual maps (GEE)

Upload two point tables to `ASSET_ROOT`: `training_points` (the field points
with their labels after step 08 of the R workflow) and
`training_points_morphology` (`result/training_MORPH3_points_ALL.csv`). Then run:

| Script | Output assets |
|---|---|
| `02_training_points.js` | `training_gse_<year>` |
| `03_training_points_morphology.js` | `training_gse_morphology_<year>` |
| `04_presence_model.js` | `presence_prob_<year>`, `presence_mask_<year>` |
| `05_seagrass_extent.js` | `seagrass_extent_2017_2024` (presence probability of at least 0.8 in any year) |
| `06_morphology_model.js` | `morphology_prob_<year>` |
| `07_cover_model.js` | `cover_<year>` (trained on samples with seagrass present only) |
| `08_agb_model.js` | `agb_<year>` |
| `09_carbon_index_model.js` | `carbon_index_<year>` |
| `10_agc_model.js` | `agc_<year>` |

The GEE models use the same settings as the R models. They are refits, so their
values differ slightly from the R models.

### 4. Validation (R)

Run the scripts in `R/validation/` in number order. Each script starts in a
test mode (`test_mode <- TRUE`, or `ext_test_mode <- TRUE` in script 03) with
one repeat and 50 trees. When the test run works, set it to `FALSE` and run
again. Finished folds and repeats are saved as checkpoints in `result/`, so a
stopped run can be started again. A full run takes several hours per script.

| Script | Designs |
|---|---|
| `01_proxy_models_and_pixel_sharing.R` | presence, morphology and carbon index models with random-row and pixel-grouped test sets; quadrats that share a pixel |
| `02_agc_pixel_lolo_transfer.R` | pixel hold-out, unseen regions (LOLO-CV), single-region transfer; percent cover and AGB models |
| `03_agc_blocks_surveys_crossfit.R` | new survey, uncertainty of the reference AGC, proxy models trained on all samples, spatial blocks, test-to-training distances, proxy workflow in 250 m blocks, cross-fitted proxies |

### 5. Distance to field data, figures and tables

1. Run `gee/11_distance_to_field_samples.js` and copy
   `seagrass_area_by_distance_class.csv` and `seagrass_distance_percentiles.csv`
   to `result/`.
2. Run `gee/12_figure6_layers.js` (set `YEAR`, `REGION` and `CRS` in its
   configuration block) and copy the three `Fig6_*` files to `result/gee/`.
3. Run the scripts in `R/figures_tables/`. Figures are written to
   `result/paper_tables/`. The tables are built with

   ```
   python3 R/figures_tables/supplementary_tables.py result result/supplement
   ```

## Scripts behind each result of the paper

| Result | Scripts |
|---|---|
| Table 3 | `R/validation/02` (PARTS A, B, D), `R/validation/03` (PARTS F, K, L) |
| Figure 3 | `R/figures_tables/figure3_validation_designs.py` (schematic) |
| Figure 4 | `R/validation/02` and `03`, then `R/figures_tables/figure4_figure5_distance.R` |
| Figure 5 | `gee/11`, then `R/figures_tables/figure4_figure5_distance.R` |
| Figure 6 | `gee/12`, then `R/figures_tables/figure6_accuracy_map.R` (the final layout was arranged in QGIS from the same layers) |
| Table S1 | `R/validation/01` (PART B), `R/validation/02` (PART E) |
| Table S2 | species weights in `R/workflow/06_carbon_index_field.R` |
| Table S3 | grid search in `R/workflow/02`, `07`, `09`, `11`, `14` and `17` |
| Table S4 | `R/validation/03` (PART H) |
| Table S5 | `R/validation/01` (PART A) |
| Table S6 | `R/validation/02` (PART D) |
| Table S7 | `R/validation/03` (PART L) |
| Table S8 | `R/validation/03` (PART I) |
| Table S9 | `R/validation/03` (PARTS F, J), `R/validation/02` (PART A) |
| Table S10 | `R/validation/02` (PART B) |
| Table S11 | `gee/11_distance_to_field_samples.js` |
| Table S12 | `R/validation/03` (PART G) |
| Figure S1 | `R/workflow/04`, then `R/figures_tables/figureS1_presence_thresholds.R` |
| Figure S3 | `R/workflow/18` and `R/validation/03` (PART L), then `R/figures_tables/figureS3_agc_importance.R` |

Tables S1 to S10 and S12 are written as CSV files by
`R/figures_tables/supplementary_tables.py`. Figures 1 and 2, Tables 1, 2 and 4
and Figure S2 were made without the code in this repository.

## Data availability

The field data used in the paper are owned by several Indonesian government
institutions and were collected under national survey programmes. Their data
governance rules do not allow public redistribution. Published datasets are
cited in Table 2 of the paper. Other data are available from the corresponding
author on reasonable request and with permission from the data custodians. The
Google Satellite Embeddings are open in the
[GEE Data Catalog](https://developers.google.com/earth-engine/datasets).

The `.gitignore` keeps data files and outputs out of the repository.

## Related repositories

- [seagrass-agc-fieldharmonization-indonesia](https://github.com/Emha89/seagrass-agc-fieldharmonization-indonesia)
- [seagrass-spc-gse-indonesia](https://github.com/Emha89/seagrass-spc-gse-indonesia)
- [seagrass-agc-countrywide-indonesia](https://github.com/Emha89/seagrass-agc-countrywide-indonesia)

## License

MIT License, see [LICENSE](LICENSE). The field data are not covered by this
license.

## Contact

Muhammad Hafizt, School of the Environment, The University of Queensland,
Brisbane, Australia, and National Research and Innovation Agency of Indonesia
(BRIN). m.hafizt@uq.edu.au
