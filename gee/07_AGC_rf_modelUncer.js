/************************************************************
 * 07_AGC_rf_modelUncer.js
 * AGC UNCERTAINTY SCRIPT -- pixel-level AGC prediction with
 * uncertainty propagation, followed by area-level aggregation and
 * spatial autocorrelation correction (L=30m, active).
 *
 * CONFIRMED VERSION: Study 2 / Chapter 4 retains L=30m in its reported
 * Results/Discussion/Equation 12, so Section 15's autocorrelation
 * correction below is active, not disabled -- this is the version whose
 * output the paper actually reports.
 *
 * OPERATIONAL NOTE: this script computes results for ONE region and
 * ONE year per run (TARGET_LOC / TARGET_YEAR below) -- the AOI itself
 * is filtered to that single region. It is meant to be run repeatedly,
 * once per region-year combination (6 regions x 8 years = 48 runs), not
 * once for the whole study area. The 3rd_paper_Apps.js dashboard's
 * hardcoded AGC_TABLE values match this script's printed output format
 * exactly (total, CI lower, CI upper, error margins) -- that table is
 * this script's accumulated output across all 48 runs.
 *
 * PURPOSE:
 *   Pixel-level AGC prediction with uncertainty propagation,
 *   followed by area-level aggregation and spatial
 *   autocorrelation correction.
 *
 * KEY PRINCIPLE:
 *   - Training data uncertainty (SD_dataTraining) is derived
 *     from field CI and treated as a GLOBAL scalar.
 *   - Pixel-level uncertainty combines:
 *       SD_total = sqrt(SD_model^2 + SD_dataTraining^2)
 *   - SD_model is estimated via repeated random subsampling
 *     of the training set (5 iterations, 80% per iteration,
 *     without replacement) with varying RF seed per iteration,
 *     implemented within GEE for pixel-scale spatial deployment.
 *     This approximates the Monte Carlo resampling framework
 *     applied in R (100 iterations, 70/30 split) for model
 *     evaluation, adapted to GEE computational constraints.
 *     NOTE: SD_model here represents model variance (prediction
 *     spread across iterations), not prediction error against
 *     a holdout set. Higher iteration counts are feasible
 *     when computational resources permit.
 *   - Area-level uncertainty is corrected for spatial
 *     autocorrelation (L = 30 m, three Sentinel-2 pixels).
 *
 * FINAL OUTPUT FOR PAPER:
 *   - totalCarbon_ton
 *   - totalCI_lower_adj, totalCI_upper_adj
 *
 * PIPELINE POSITION: reads training_embedding_depth3_<year> plus the
 * five upstream model outputs (PA_prob, Morph3, SPC, AGB, cIndex),
 * matching the same asset names as 07_AGC_rf_modelDev.js.
 *
 * ASSET PORTABILITY: this script originally referenced the author's own
 * private GEE assets, all under one project/folder root -- replace
 * ASSET_ROOT below with your own GEE project and asset folder. The
 * Allen Coral Atlas distToLand layer is a separate, externally-hosted
 * shared asset (not the author's own) -- using it requires your own
 * access to that specific Coral Atlas asset, not just a path swap.
 ************************************************************/


// ==========================================================
// 1. CONFIGURATION
// ==========================================================

var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

// Response variable name (used consistently)
var RESPONSE    = 'AGC_pred';

// Spatial resolution (Sentinel-2)
var SCALE       = 10;

// Random seed for reproducibility
var SEED        = 42;

// RF ensemble size (kept for compatibility; now used as MC_ITER)
var N_MODELS    = 5;

// Monte Carlo settings (keep lightweight to avoid memory issues)
var MC_ITER     = N_MODELS; // alias to preserve script structure
var MC_FRACTION = 0.80;     // fraction of training used per MC iteration

// Reporting unit -- change these and re-run for each region/year
var TARGET_LOC  = 'Ayau';
var TARGET_YEAR = 2024;

// Training data quality control
var MAX_CI_WIDTH = 100;  // gC/m2


// ==========================================================
// 2. DATA PATHS
// ==========================================================

var TRAIN_PREFIX   = ASSET_ROOT + '/output/training_embedding_depth3_';
var PA_PROB_PREFIX = ASSET_ROOT + '/output/RF_probability04022026_';
var SPC_PREFIX     = ASSET_ROOT + '/output/RF_tSPC04022026_';
var AGB_PREFIX     = ASSET_ROOT + '/output/RF_AGB04022026_';
var CINDEX_PREFIX  = ASSET_ROOT + '/output/RF_cIndex04022026_';

// Morphology probability rasters
var MORPH_PREFIX   = ASSET_ROOT + '/output/MORPH3_probs04022026_';
var MORPH_BANDS = [
  'P_mixed_long',
  'P_mixed_short_plus_mono_short',
  'P_mono_Ea'
];

// Seagrass persistence mask (analysis domain). An alternative
// (byYearCount-based) mask was considered -- see the commented-out line
// below, kept for reference -- but the MaxProb-based mask is the one
// actually used.
//var PERSISTENCE_MASK = ee.Image(
//  ASSET_ROOT + '/output/seagrass_persistence_mask_byYearCount_2017_2025'
//);

var PERSISTENCE_MASK = ee.Image(
  ASSET_ROOT + '/output/Seagrass_MaxProb_0_8_2017_2024_mask'
);

// Export destination
var EXPORT_PREFIX  = ASSET_ROOT + '/output/RF_AGC_UNC_';


// ==========================================================
// 3. RANDOM FOREST PARAMETERS
// ==========================================================

var RF_PARAMS = {
  numberOfTrees: 900,
  variablesPerSplit: 5,
  minLeafPopulation: 9,
  bagFraction: 0.6,
  seed: SEED
};


// ==========================================================
// 4. AOI DEFINITION
// ==========================================================

var AOI = ee.FeatureCollection(
  ASSET_ROOT + '/InputArea/R2_case_study'
).filter(ee.Filter.eq('loc', TARGET_LOC)).geometry();

Map.centerObject(AOI, 12);


// ==========================================================
// 5. FEATURE SETUP
// ==========================================================

// Google Satellite Embedding bands
var GSE_BANDS = ee.List.sequence(0, 63).map(function(i) {
  return ee.String('A').cat(ee.Number(i).format('%02d'));
});

// Full predictor list used by AGC RF
var FEATURE_BANDS = GSE_BANDS
  .cat(['depth','distToLand'])
  .add('PA_prob')
  .cat(MORPH_BANDS)
  .add('tSPC_pred')
  .add('tAGB_pred')
  .add('carbon_index');

print('Predictor bands:', FEATURE_BANDS);


// ==========================================================
// 6. HELPER FUNCTIONS
// ==========================================================

// Generic image loader
function loadImage(prefix, year, name) {
  return ee.Image(prefix + year).rename(name);
}

// Carbon index loader
function loadCarbonIndex(year) {
  return ee.Image(CINDEX_PREFIX + year).rename('carbon_index');
}

// Morphology probability loader
function loadMorphProbs(year) {
  return ee.Image(MORPH_PREFIX + year).select(MORPH_BANDS);
}

// Safe parsing of numeric AGC fields
function cleanAGC(f) {
  function safeParse(v) {
    return ee.Algorithms.If(
      ee.Algorithms.IsEqual(v, 'NA'),
      null,
      ee.Algorithms.If(v === null, null, ee.Number.parse(v))
    );
  }
  return f
    .set(RESPONSE,     safeParse(f.get(RESPONSE)))
    .set('AGC_up',     safeParse(f.get('AGC_up')))
    .set('AGC_low',    safeParse(f.get('AGC_low')))
    .set('AGC_CIwidt', safeParse(f.get('AGC_CIwidt')));
}

// Monte Carlo subsampling (without replacement, stable & lightweight)
function mcSampleTraining(fc, iter, fraction) {
  var withRand = fc.randomColumn('mc_rand', ee.Number(SEED).add(iter));
  var n = fc.size();
  var k = ee.Number(n).multiply(fraction).int();
  return withRand.sort('mc_rand').limit(k);
}


// ==========================================================
// 7. BUILD TRAINING DATASET
// ==========================================================

var TRAIN_YEARS = [2018,2019,2020,2021,2022,2023];

var training_raw = ee.FeatureCollection(TRAIN_YEARS.map(function(year){

  var pts = ee.FeatureCollection(TRAIN_PREFIX + year);

  var layers = loadImage(PA_PROB_PREFIX, year, 'PA_prob')
    .addBands(loadMorphProbs(year))
    .addBands(loadImage(SPC_PREFIX, year, 'tSPC_pred'))
    .addBands(loadImage(AGB_PREFIX, year, 'tAGB_pred'))
    .addBands(loadCarbonIndex(year));

  return layers.sampleRegions({
    collection: pts,
    properties: pts.first().propertyNames(),
    scale: SCALE,
    geometries: true
  }).map(cleanAGC);

})).flatten()
  .filter(ee.Filter.notNull([RESPONSE,'AGC_low','AGC_up','AGC_CIwidt']))
  .filter(ee.Filter.gt(RESPONSE,0))
  .filter(ee.Filter.lte(RESPONSE,2000));


// ==========================================================
// 8. TRAINING CI DIAGNOSTICS + FILTERING
// ==========================================================

var ciReducers = ee.Reducer.min()
  .combine(ee.Reducer.max(),null,true)
  .combine(ee.Reducer.mean(),null,true)
  .combine(ee.Reducer.stdDev(),null,true)
  .combine(ee.Reducer.percentile([25,50,75]),null,true)
  .combine(ee.Reducer.count(),null,true);

var ciStats_before = training_raw.reduceColumns({
  reducer: ciReducers,
  selectors: ['AGC_CIwidt']
});

print('CI width BEFORE filter:', ciStats_before);

var training = training_raw.filter(
  ee.Filter.lte('AGC_CIwidt', MAX_CI_WIDTH)
);

print('Training size after CI filter:', training.size());


// ==========================================================
// 9. TRAINING DATA UNCERTAINTY (GLOBAL SCALAR)
// ==========================================================

var trainingWithSD = training.map(function(f){
  var sd = ee.Number(f.get('AGC_up'))
    .subtract(ee.Number(f.get('AGC_low')))
    .divide(2 * 1.96);
  return f.set('SD_data', sd);
});

// IMPORTANT: single scalar value, computed as the median (not the mean)
// for robustness against outliers
var SD_dataTraining = ee.Number(
  trainingWithSD.reduceColumns({
    reducer: ee.Reducer.median(),
    selectors: ['SD_data']
  }).get('median')
);

print('SD_dataTraining (scalar):', SD_dataTraining);

var sdSummary = trainingWithSD.reduceColumns({
  reducer: ee.Reducer.percentile([10, 50, 90])
    .combine(ee.Reducer.mean(), null, true)
    .combine(ee.Reducer.median(), null, true),
  selectors: ['SD_data']
});

print('SD_data distribution (p10/p50/p90/mean/median):', sdSummary);

// ==========================================================
// 10. PREDICTOR STACK (TARGET YEAR)
// ==========================================================

var embImg = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
  .filterDate(ee.Date.fromYMD(TARGET_YEAR,1,1),
              ee.Date.fromYMD(TARGET_YEAR+1,1,1))
  .mosaic();

var depthImg = ee.Image(
  ASSET_ROOT + '/InputArea/raster_INA_ACA_Bathymetry'
).rename('depth').multiply(-1);

// Distance to land. EXTERNAL ASSET: hosted by the Allen Coral Atlas
// team, not the author's own -- requires your own access to this
// specific shared asset, not just a path change.
var distImg = ee.Image(
  'projects/coral_atlas/global_datasets/osm_distToLand_indo'
).rename('distToLand');

var stack = embImg.select(GSE_BANDS)
  .addBands(depthImg)
  .addBands(distImg)
  .addBands(loadImage(PA_PROB_PREFIX, TARGET_YEAR,'PA_prob'))
  .addBands(loadMorphProbs(TARGET_YEAR))
  .addBands(loadImage(SPC_PREFIX, TARGET_YEAR,'tSPC_pred'))
  .addBands(loadImage(AGB_PREFIX, TARGET_YEAR,'tAGB_pred'))
  .addBands(loadCarbonIndex(TARGET_YEAR))
  .updateMask(PERSISTENCE_MASK)
  .clip(AOI);


// ==========================================================
// 11. ENSEMBLE RF PREDICTION (MONTE CARLO)
// ==========================================================

function ensemblePrediction(responseField){
  return ee.ImageCollection(
    ee.List.sequence(1, MC_ITER).map(function(i){
      i = ee.Number(i);

      // Monte Carlo subsample of training data
      var train_i = mcSampleTraining(training, i, MC_FRACTION);

      // Vary RF seed per iteration
      var rf_i = {
        numberOfTrees: RF_PARAMS.numberOfTrees,
        variablesPerSplit: RF_PARAMS.variablesPerSplit,
        minLeafPopulation: RF_PARAMS.minLeafPopulation,
        bagFraction: RF_PARAMS.bagFraction,
        seed: ee.Number(SEED).add(i).int()
      };

      var model = ee.Classifier.smileRandomForest(rf_i)
        .setOutputMode('REGRESSION')
        .train({
          features: train_i,
          classProperty: responseField,
          inputProperties: FEATURE_BANDS
        });

      return stack.classify(model)
        .rename('AGC_pred')
        .set('mc_iter', i);
    })
  );
}

var predEns = ensemblePrediction(RESPONSE);


// ==========================================================
// 12. PIXEL-LEVEL UNCERTAINTY
// ==========================================================

var mean     = predEns.reduce(ee.Reducer.mean()).rename('AGC_mean');
var SD_model = predEns.reduce(ee.Reducer.stdDev()).rename('SD_model_MC');

var SD_total = SD_model.pow(2)
  .add(ee.Image.constant(SD_dataTraining).pow(2))
  .sqrt()
  .rename('AGC_totalUncertainty');

// Mean SD_model over AOI
var mean_SD_model = SD_model.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: AOI,
  scale: SCALE,
  maxPixels: 1e13
});

print('Mean SD_model (gC/m2):', mean_SD_model);


// ==========================================================
// 13. PIXEL-LEVEL CI (FOR MAPPING ONLY)
// ==========================================================

var z = 1.96;

var CI_lower = mean.subtract(SD_total.multiply(z))
  .rename('AGC_CI_lower')
  .updateMask(PERSISTENCE_MASK);

var CI_upper = mean.add(SD_total.multiply(z))
  .rename('AGC_CI_upper')
  .updateMask(PERSISTENCE_MASK);

var CI_width = CI_upper.subtract(CI_lower)
  .rename('AGC_CI_width')
  .updateMask(PERSISTENCE_MASK);


// ==========================================================
// 14. AREA-LEVEL AGGREGATION (NO AUTOCORRELATION YET)
// ==========================================================

var pixelArea = ee.Image.pixelArea();

var totalCarbon_g = ee.Number(
  mean.multiply(pixelArea)
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: AOI,
      scale: SCALE,
      maxPixels: 1e13
    })
    .values()
    .get(0)
);

var totalCarbon_ton = totalCarbon_g.divide(1e6);


var totalVar_g2 = ee.Number(
  SD_total.pow(2)
    .multiply(pixelArea.pow(2))
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: AOI,
      scale: SCALE,
      maxPixels: 1e13
    })
    .values()
    .get(0)
);

var totalSD_ton = totalVar_g2.sqrt().divide(1e6);



// ==========================================================
// 15. SPATIAL AUTOCORRELATION CORRECTION (FINAL)
// ==========================================================

// Correlation length = 3x3 Sentinel-2 pixels
var L = 30;

// Seagrass area
var seagrassArea_m2 = PERSISTENCE_MASK
  .multiply(pixelArea)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: AOI,
    scale: SCALE,
    maxPixels: 1e13
  }).values().get(0);

// Number of pixels
var N_pixel = ee.Number(seagrassArea_m2).divide(SCALE * SCALE);

// Effective independent samples
var A_corr = Math.PI * L * L;
var N_eff  = ee.Number(seagrassArea_m2).divide(A_corr);

// Variance inflation factor
var inflationFactor = N_pixel.divide(N_eff).sqrt();

// Adjusted SD and CI (FINAL FOR PAPER)
var totalSD_ton_adj = totalSD_ton.multiply(inflationFactor);

var totalCI_lower_adj =
  ee.Number(totalCarbon_ton).subtract(totalSD_ton_adj.multiply(1.96));

var totalCI_upper_adj =
  ee.Number(totalCarbon_ton).add(totalSD_ton_adj.multiply(1.96));


// ==========================================================
// 16. FINAL PRINT (USE THESE VALUES)
// ==========================================================

print('===== FINAL RESULTS =====');
print('Total AGC (ton C):', totalCarbon_ton);
print('95% CI lower (adj):', totalCI_lower_adj);
print('95% CI upper (adj):', totalCI_upper_adj);
print('========================================');


// ==========================================================
// 17. EXPORTS
// ==========================================================

Export.image.toAsset({
  image: mean,
  description: 'AGC_mean_' + TARGET_LOC + '_' + TARGET_YEAR,
  assetId: EXPORT_PREFIX + 'mean_' + TARGET_LOC + '_' + TARGET_YEAR,
  region: AOI,
  scale: SCALE,
  maxPixels: 1e13
});

Export.image.toAsset({
  image: SD_total,
  description: 'AGC_total_uncertainty_' + TARGET_LOC + '_' + TARGET_YEAR,
  assetId: EXPORT_PREFIX + 'total_unc_' + TARGET_LOC + '_' + TARGET_YEAR,
  region: AOI,
  scale: SCALE,
  maxPixels: 1e13
});

Export.image.toAsset({
  image: CI_width,
  description: 'AGC_CIwidth_' + TARGET_LOC + '_' + TARGET_YEAR,
  assetId: EXPORT_PREFIX + 'CIwidth_' + TARGET_LOC + '_' + TARGET_YEAR,
  region: AOI,
  scale: SCALE,
  maxPixels: 1e13
});

print('EXPORT QUEUED for', TARGET_LOC, TARGET_YEAR);


//-----------------------------------------------------------
// ==========================================================
// 13B. PIXEL-LEVEL ERROR VARIATION VISUALISATION
// ==========================================================

// ----------------------------------------------------------
// 1. Absolute uncertainty (gC/m2)
// ----------------------------------------------------------

var absError = SD_total.rename('AGC_absError')
  .updateMask(PERSISTENCE_MASK);

// Visual parameters (adjust if needed)
var visAbs = {
  min: 0,
  max: 40,   // adjust based on your SD range
  palette: ['white','yellow','orange','red']
};

Map.addLayer(absError, visAbs, 'Pixel Absolute Uncertainty (SD_total)',  false);


// ----------------------------------------------------------
// 2. Relative uncertainty (Coefficient of Variation)
//    CV = SD_total / mean
// ----------------------------------------------------------

var CV = SD_total
  .divide(mean)
  .rename('AGC_CV')
  .updateMask(PERSISTENCE_MASK);

// Avoid division artefact where mean ~ 0
CV = CV.updateMask(mean.gt(1));

var visCV = {
  min: 0,
  max: 0.6,  // 0-60% relative uncertainty
  palette: ['blue','cyan','yellow','orange','red']
};

Map.addLayer(CV, visCV, 'Pixel Relative Uncertainty (CV)', false);