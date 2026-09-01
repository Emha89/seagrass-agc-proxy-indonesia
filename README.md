# Seagrass AGC Proxy Pipeline (Indonesia) -- Study 2

Code accompanying: *A Proxy-Based Strategy for Estimating Seagrass Above-Ground
Carbon from Satellite Observations* (Hafizt et al., in preparation).

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
│   ├── 07_validation/          cross-site and LOLO-CV validation of AGC
│   ├── 08_uncertainty/         spatial uncertainty / correlation length (L)
│   └── unsorted_need_confirmation/   role not yet confirmed, see Open items
├── gee/                        Google Earth Engine deployment script(s)
└── data/                       not included -- see Data availability below
```

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

1. `01_occurrence_PA/dataPrep_func_PA.R` -- **pass 1** (produces
   `dataPA_<year>.csv` without carbon_index / morphology columns yet)
2. `01_occurrence_PA/rf_model_PA.R` then `rf_applyModel_PA.R` (produces
   `predicted_PA_<year>.csv`)
3. `02_morphology/rf_model_combined_MORPH3.R` then
   `rf_applyModel_MORPH3.R` (needs `predicted_PA_<year>.csv`; produces
   `predicted_MORPH3_probs_<year>.csv`)
4. `05_carbon_index/carbon_scaling_join_SP_finalDell.R` (needs
   `dataPA_<year>.csv` from pass 1, plus raw `dataSP_<year>.csv`;
   produces `dataCARBON_<year>.csv`)
5. `01_occurrence_PA/dataPrep_func_PA.R` again -- **pass 2** (now merges
   in carbon_index and the morphology P_* columns)
6. `03_SPC/rf_model_combined_SPC.R` then `rf_applyModel_SPC.R` (training
   needs the pass-2 `dataPA_<year>.csv`; apply needs `predicted_PA_<year>.csv`,
   which already carries morphology P_* via step 3's own merge)
7. `04_AGB/dataPrep_funcGSE_AGB2.R`, `rf_model_AGB2.R`, `rf_applyModel_AGB.R`
8. `05_carbon_index/rf_model_combined_CINDEX.R` then `rf_applyModel_cIndex.R`
9. `06_AGC_final/dataPrep_funcGSE_AGC2.R`, `rf_model_AGC3.R`, `rf_applyModel_AGC3.R`
10. Any `rf_evalModel_*.R` script can run any time after its matching
    apply-model step above.

`08_uncertainty/recompute_L_study2.R` and `07_validation/` can run once
the AGC stage (step 9) is complete.

## Data availability

This repository contains **analysis code only**. Field survey data compiled
from third-party sources are subject to the data-sharing policies of the
original providers and are not included here (see the paper's Data
Availability statement).

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
`sg_morpho`), used by the carbon index stage.

## Reproducing

1. Open `seagrass-agc-proxy-indonesia.Rproj` in RStudio. This sets the
   working directory automatically so every script's `here()` calls resolve
   to the right folder, regardless of where the repository sits on disk.
2. Install required packages: `ranger`, `randomForest`, `caret`, `dplyr`,
   `readr`, `sf`, `gstat`, `blockCV`, `ggplot2`, `here`.
3. Run the R scripts following the execution order above (not simply
   folder 01 through 08 in sequence -- see that section for why).
4. GEE scripts in `gee/` run separately in the
   [GEE Code Editor](https://code.earthengine.google.com) -- see comments in
   each script for GEE asset paths that need to point at your own account.

## Open items

- [x] ~~AGC data prep: `dataPrep_funcGSE_AGC2.R` vs `_AGC3.R`~~ -- resolved:
      `AGC2` confirmed (`rf_model_AGC3.R` sources it by name)
- [x] ~~SPC model: `rf_model_SPC.R` vs `rf_model_SPC2.R`~~ -- resolved: the
      actual script is `rf_model_combined_SPC.R` (predictors are
      morphology probabilities, not carbon_index -- an earlier draft
      using carbon_index was superseded)
- [ ] Confirm the final version for the remaining scripts with more than
      one candidate on disk:
  - Occurrence data prep: `dataPrep_func_PA.R` vs `dataPrep_func_PA2.R`
    (only the first was provided; PA2 not yet reviewed)
  - Spatial uncertainty: `recompute_L_study2.R` (included) vs `recompute_L_v2.R`
    vs an older `recompute L.R`
- [ ] Confirm whether `dataPrep_SPECIES.R`, `rf_model_SP.R`,
      `rf_applyModel_SP.R` belong in the pipeline or are superseded
- [ ] Add the GEE deployment script(s)
- [ ] Add the paper citation once accepted
- [ ] Choose a licence (MIT is a common default for research code; the data
      itself is governed separately by the third-party providers' terms)
- [ ] Build a field-data Excel template: one sheet showing the expected
      column headers (matched to what the scripts actually reference --
      gee_id, xcoord/ycoord, year, loc, PA, tSPC, AGB_pred, AGC_pred,
      species cover columns like `Ea_SPC`, sg_morpho, etc.), so anyone
      reproducing the pipeline with their own field data knows the
      required structure without needing the original data itself. All
      R scripts have now been reviewed, so this can be compiled from the
      column names documented across R/.
- Note: in dataPrep_func_PA.R, the xcoord/ycoord that survive into
      dataPA_<year>.csv come from the GSE-extraction side, not the GT
      field-data side (gt_sf's own coordinates are dropped earlier via
      st_as_sf() without remove = FALSE). Left as-is deliberately: the
      exact reasoning isn't fully recalled, but this is the version that
      already produced the Study 2 / Chapter 4 manuscript results, so it
      should not be changed without also re-checking those results.

## Citation

_(to be added once the paper is accepted)_