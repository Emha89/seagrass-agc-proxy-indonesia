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
rf_applyModel_* -> rf_evalModel_*` pattern.

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

## Reproducing

1. Open `seagrass-agc-proxy-indonesia.Rproj` in RStudio. This sets the
   working directory automatically so every script's `here()` calls resolve
   to the right folder, regardless of where the repository sits on disk.
2. Install required packages: `ranger`, `randomForest`, `caret`, `dplyr`,
   `readr`, `sf`, `gstat`, `blockCV`, `ggplot2`, `here`.
3. Run the R scripts in numbered folder order (01 through 08).
4. GEE scripts in `gee/` run separately in the
   [GEE Code Editor](https://code.earthengine.google.com) -- see comments in
   each script for GEE asset paths that need to point at your own account.

## Open items

- [ ] Confirm the final version for scripts with more than one candidate on disk:
  - Occurrence data prep: `dataPrep_func_PA.R` vs `dataPrep_func_PA2.R`
  - SPC model: `rf_model_SPC.R` vs `rf_model_SPC2.R`
  - AGC data prep: `dataPrep_funcGSE_AGC2.R` vs `dataPrep_funcGSE_AGC3.R`
  - Spatial uncertainty: `recompute_L_study2.R` (included) vs `recompute_L_v2.R`
    vs an older `recompute L.R`
- [ ] Confirm whether `dataPrep_SPECIES.R`, `rf_model_SP.R`,
      `rf_applyModel_SP.R` belong in the pipeline or are superseded
- [ ] Add the GEE deployment script(s)
- [ ] Add the paper citation once accepted
- [ ] Choose a licence (MIT is a common default for research code; the data
      itself is governed separately by the third-party providers' terms)

## Citation

_(to be added once the paper is accepted)_
