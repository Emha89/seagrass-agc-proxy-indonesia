# Seagrass Above-Ground Carbon Estimation Using a Proxy-Based Framework

R and Google Earth Engine scripts for:

Hafizt, M., Adi, N.S., Lyons, M., Phinn, S., McMahon, K., & Roelfsema, C.
(2026). **A Proxy-Based Strategy for Estimating Seagrass Above-Ground Carbon
from Satellite Observations**.

Companion repository for Study 1 (GSE spectral analysis):
[seagrass-spc-gse-indonesia](https://github.com/Emha89/seagrass-spc-gse-indonesia)

## Overview

Six sequential Random Forest proxy stages predict above-ground carbon (AGC)
from Google Satellite Embeddings (GSE) and field survey data:

`occurrence probability -> leaf morphology probability -> percent cover (SPC)
-> above-ground biomass (AGB) -> species-weighted carbon index -> final AGC`

Models are developed and evaluated in R (Monte Carlo resampling, leave-one-
location-out cross-validation across six regions), then the tuned models are
deployed in Google Earth Engine (GEE) to generate annual AGC maps (2017-2024)
with pixel-level uncertainty.

## Repository structure

```
seagrass-agc-proxy-indonesia/
├── R/
│   ├── 00_shared_functions/    shared plotting / classification / regression helpers
│   ├── 01_occurrence_PA/       stage 1 -- seagrass occurrence probability
│   ├── 02_morphology/          stage 2 -- leaf morphology probability
│   ├── 03_SPC/                 stage 3 -- percent cover
│   ├── 04_AGB/                 stage 4 -- above-ground biomass
│   ├── 05_carbon_index/        stage 5 -- species-weighted carbon index
│   ├── 06_AGC_final/           stage 6 -- final AGC model
│   └── 07_validation/          cross-site and LOLO-CV validation of AGC
├── gee/                         Google Earth Engine deployment scripts
│   ├── 00_extractGSE_trainingPoints.js   samples GSE+terrain at the
│   │                                      ORIGINAL training points --
│   │                                      true source of
│   │                                      GSE_training_<year>_CSV.csv
│   ├── 01_trainingData_prep.js           re-samples GSE at REFINED
│   │                                      (post-R) PA points -- runs
│   │                                      after the R pipeline, not before
│   ├── 01_trainingData_prepMorpho.js     same pattern, for morphology
│   ├── 02_PROB_rf_modelDev.js            PA model, trained + deployed
│   │                                      natively in GEE
│   ├── 02_extent_2017-2025.js            seagrass persistence/extent
│   │                                      maps from PA probability
│   ├── 03_MORPHO_rf_modelDev.js          morphology model, GEE-native
│   ├── 04_COVER_rf_modelDev.js           percent cover model, GEE-native
│   ├── 05_AGB_rf_modelDev.js             AGB model, GEE-native
│   ├── 06_cIndex_rf_modelDev.js          carbon index model, GEE-native
│   ├── 07_AGC_rf_modelDev.js             final AGC model, all proxies
│   │                                      as predictors
│   ├── 07_AGC_rf_modelUncer.js           per-region/year AGC total +
│   │                                      uncertainty (Monte Carlo +
│   │                                      L=30m autocorrelation
│   │                                      correction, active)
│   └── 3rd_paper_Apps.js                 interactive results-preview
│                                          dashboard
└── data/                       field_data_template.xlsx only -- actual data
                                 not included, see Data availability below
```

Only code that the reported Results/Discussion actually rely on is published
here. Exploratory analyses that were tried and not adopted are kept in the
author's own local working copy rather than this repository -- see
"Exploratory work not included" below for what that covers and why.

Each stage folder generally follows a `dataPrep_* -> rf_model_* ->
rf_applyModel_* -> rf_evalModel_*` pattern, except 02_morphology/, which
is a single combined script (data prep, training, prediction, and
evaluation together) plus a separate apply-model and eval script.

## Suggested execution order

The folder numbers (01-06) describe the proxy chain described in the
paper, but do NOT match the order scripts actually need to run in --
several stages depend on outputs from stages that come later in the
numbering. This was worked out step by step while reviewing each script;
the actual required order is:

1. `01_occurrence_PA/dataPrep_func_PA2.R` -- **pass 1** (produces
   `dataPA_<year>.csv` without carbon_index / morphology columns yet)
2. `01_occurrence_PA/rf_model_PA.R` then `rf_applyModel_PA.R` (produces
   `predicted_PA_<year>.csv`)
3. `02_morphology/rf_model_combined_MORPH3.R` then
   `rf_applyModel_MORPH3.R` (needs `predicted_PA_<year>.csv`; produces
   `predicted_MORPH3_probs_<year>.csv`)
4. `05_carbon_index/dataPrep_mergeSpecies.R` (needs `dataPA_<year>.csv`
   from pass 1, plus raw `sg_compos.csv`; produces `dataSP_<year>.csv`)
5. `05_carbon_index/carbon_scaling_join_SP_finalDell.R` (needs
   `dataPA_<year>.csv` from pass 1, plus `dataSP_<year>.csv` from step 4;
   produces `dataCARBON_<year>.csv`)
6. `01_occurrence_PA/dataPrep_func_PA2.R` again -- **pass 2** (now merges
   in carbon_index and the morphology P_* columns)
7. `03_SPC/rf_model_combined_SPC.R` then `rf_applyModel_SPC.R` (training
   needs the pass-2 `dataPA_<year>.csv`; apply needs `predicted_PA_<year>.csv`,
   which already carries morphology P_* via step 3's own merge)
8. `04_AGB/dataPrep_funcGSE_AGB2.R`, `rf_model_AGB2.R`, `rf_applyModel_AGB.R`
9. `05_carbon_index/rf_model_combined_CINDEX.R` then `rf_applyModel_cIndex.R`
10. `06_AGC_final/dataPrep_funcGSE_AGC2.R`, `rf_model_AGC3.R`, `rf_applyModel_AGC3.R`
11. Any `rf_evalModel_*.R` script, and `07_validation/`, can run any time
    after their matching apply-model step above.

## GEE execution order

Like the R-side pipeline, the GEE scripts have a dependency chain that
isn't simply "run 00 through 07 in order" -- the two 07 scripts each
depend on all five upstream model outputs existing first:

1. `00_extractGSE_trainingPoints.js` -- samples GSE+terrain at the
   original training points, exports `GSE_training_<year>_CSV.csv` (the
   file the R pipeline is built on)
2. *(the full R pipeline runs here -- see the R execution order above)*
3. `01_trainingData_prep.js` and `01_trainingData_prepMorpho.js` --
   re-sample GSE at refined/filtered points that have been through the R
   analysis and re-uploaded to GEE, producing `training_embedding_depth3_
   <year>` and `training_morph3_<year>` assets
4. `02_PROB_rf_modelDev.js` -- trains + deploys the PA model, producing
   `RF_probability04022026_<year>`
5. `02_extent_2017-2025.js` -- persistence/extent masks from step 4's
   output
6. `03_MORPHO_rf_modelDev.js` -- needs step 3's `training_morph3_<year>`;
   produces `MORPH3_probs04022026_<year>`
7. `04_COVER_rf_modelDev.js` -- needs step 6's morphology output;
   produces `RF_tSPC04022026_<year>`
8. `05_AGB_rf_modelDev.js` -- independent of steps 6-7 (GSE+depth+
   distToLand only); produces `RF_AGB04022026_<year>`
9. `06_cIndex_rf_modelDev.js` -- independent of steps 6-8; produces
   `RF_cIndex04022026_<year>`
10. `07_AGC_rf_modelDev.js` -- needs all of steps 4, 6, 7, 8, 9;
    produces `RF_AGC04022026_<year>`
11. `07_AGC_rf_modelUncer.js` -- needs the same five upstream outputs as
    step 10, plus step 5's persistence mask; run once per region-year
    combination (`TARGET_LOC` / `TARGET_YEAR` at the top of the script),
    48 runs total (6 regions x 8 years)
12. `3rd_paper_Apps.js` -- interactive viewer; needs all of the above,
    plus the manually-compiled `AGC_TABLE` inside the script itself (the
    accumulated printed output of step 11's 48 runs)

**Asset naming**: every script defines `ASSET_ROOT` at the top and
builds every other path from it (e.g. `ASSET_ROOT + '/output/
RF_probability04022026_' + year`) -- replace `ASSET_ROOT` in each script
with your own GEE project/asset folder and the whole chain resolves
consistently. The Allen Coral Atlas `distToLand` layer is a separate,
externally-hosted shared asset in every script that uses it -- that one
needs its own access request, not just a path change.

**Known naming inconsistency**: `03_MORPHO`/`04_COVER`/`05_AGB`/
`06_cIndex` all reference a persistence mask named `..._2017_2024_mask`,
but `02_extent_2017-2025.js` (kept exactly as originally written)
exports that same mask as `..._2017_2025_mask`. This only affects an
optional map-preview layer in each script, not the actual exported model
rasters (which are unmasked) -- but if setting this up from scratch,
either rename the exported asset or update the four references so they
match. `07_AGC_rf_modelDev.js` and `07_AGC_rf_modelUncer.js` each
reference a different mask (byYearCount and MaxProb respectively) that
does line up with `02_extent_2017-2025.js`'s own naming -- see the
comments in each file for details.

## Data availability

This repository contains **analysis code only**. Field survey data compiled
from third-party sources are subject to the data-sharing policies of the
original providers and are not included here (see the paper's Data
Availability statement).

`data/field_data_template.xlsx` documents the exact column headers each raw
input file needs (GT_Label_Data, GSE_Training_PerYear, Species_Composition
sheets), with a description and example value for every column -- use it
to structure your own data. The sections below summarise the same
information in text form.

To run this code with your own data, place a training file at:

```
data/carbon_scaling_join_SP_finalDell.csv
```

Expected columns include coordinate fields (`xcoord`, `ycoord`), response
variables (e.g. `AGC_pred`), 64 GSE embedding bands (`GSE_A00`-`GSE_A63`),
and proxy predictor columns (`carbon_index`, `tSPC_pred`, `tAGB_pred`,
morphology probability columns).

Also expected under `data/`: yearly `dataSP_<year>.csv` files (raw field
species percent-cover data with columns like `Ea_SPC`, `Th_SPC`, ...,
`sg_morpho`), and `sg_compos.csv` (species composition survey data, keyed
by `compositio`) -- both used by the carbon index stage.

## Reproducing

1. Open `seagrass-agc-proxy-indonesia.Rproj` in RStudio. This sets the
   working directory automatically so every script's `here()` calls resolve
   to the right folder, regardless of where the repository sits on disk.
2. Install required packages: `ranger`, `randomForest`, `caret`, `dplyr`,
   `readr`, `sf`, `ggplot2`, `here`.
3. Run the R scripts following the execution order above (not simply
   folder 01 through 07 in sequence -- see that section for why).
4. GEE scripts in `gee/` run separately in the
   [GEE Code Editor](https://code.earthengine.google.com) -- see the GEE
   execution order above, and each script's own header comments for
   asset-path placeholders that need to point at your own account.

## Exploratory work not included

Three supplementary analyses were carried out locally but are not part of
this repository, because none of them fed into the reported
Results/Discussion:

- **Per-region spatial correlation length (L)**: a variogram-based
  re-estimate of L, fitted separately for each of the 6 regions. Failed
  to converge for 4 of the 5 testable regions (fitted L exceeded the
  region's own physical area).
- **Cluster design effect / ICC (DEFF)**: an alternative effective-sample-
  size estimate treating each region as a survey cluster (k=6). Gave a
  numerically plausible pooled result, but 6 clusters is too few for a
  reliable ICC estimate on its own.
- **Species-type classification**: an earlier attempt to train an RF
  classifier predicting detailed seagrass species type/dominant species
  (`sg_type_co` / `sg_dom`) directly from GSE bands, plus exploratory
  per-species spectral band plots. The paper's adopted approach instead
  uses the 3-class morphology probability (02_morphology/) as the proxy
  for species composition -- confirmed directly in the manuscript text:
  *"the implementation in this study incorporated leaf morphology
  probability classes as a proxy for species composition."* The species
  presence/dominance summary report these scripts also produced is not
  referenced by any other script in the adopted pipeline either.

None of these replaced what's reported in the paper. The first two
informed a two-sentence addition to Study 2's Limitations section (4.4);
the third was superseded by the morphology-classifier approach entirely.
The companion national model in Chapter 5 *did* adopt a DEFF/ICC
replacement for its own uncertainty estimate (33 clusters via DBSCAN,
N_effective=201.3) -- that lives in a separate repository for that
chapter, not here.

## Notes

In dataPrep_func_PA2.R, the xcoord/ycoord that survive into
dataPA_<year>.csv come from the GSE-extraction side, not the GT
field-data side (gt_sf's own coordinates are dropped earlier via
st_as_sf() without remove = FALSE). Left as-is deliberately: the
exact reasoning isn't fully recalled, but this is the version that
already produced the Study 2 / Chapter 4 manuscript results, so it
should not be changed without also re-checking those results.

## Interactive AGC Viewer

A complete time series of annual AGC predictions (2017-2024) across all six
study regions, plus every intermediate proxy layer and pixel-level
uncertainty, is available at:
https://muhammadhafizt.users.earthengine.app/view/seagrassagclocalsindonesia

This app is built from `gee/3rd_paper_Apps.js` in this repository.

## Citation

_(to be added once the paper is accepted)_

## License

This project is licensed under the MIT License -- see the
[LICENSE](LICENSE) file for details. Field data used to produce the
results are not covered by this license -- see Data availability above.

## Contact

Muhammad Hafizt
School of the Environment, The University of Queensland, Brisbane, Australia
National Research and Innovation Agency of Indonesia (BRIN)
m.hafizt@uq.edu.au