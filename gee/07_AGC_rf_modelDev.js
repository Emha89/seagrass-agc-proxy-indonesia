/************************************************************
 * 07_AGC_rf_modelDev.js
 * RF regression for AGC (above-ground carbon) -- the final stage,
 * synthesizing all upstream proxy layers as predictors.
 *
 * Predictors:
 *   - GSE embeddings (64 bands)
 *   - Static layers (depth, distToLand)
 *   - Sequential proxy layers: PA_prob -> Morph3 probs -> tSPC -> tAGB
 *     -> carbon_index
 *
 * PIPELINE POSITION: reads training_embedding_depth3_<year> (from
 * 01_trainingData_prep.js) plus FIVE upstream model outputs -- all five
 * asset-prefix connections below are confirmed exact matches to their
 * source scripts' own export names, date stamps included:
 *   RF_probability04022026_  <- 02_PROB_rf_modelDev.js
 *   MORPH3_probs04022026_    <- 03_MORPHO_rf_modelDev.js
 *   RF_tSPC04022026_         <- 04_COVER_rf_modelDev.js
 *   RF_AGB04022026_          <- 05_AGB_rf_modelDev.js
 *   RF_cIndex04022026_       <- 06_cIndex_rf_modelDev.js
 * RF_PARAMS below (ntree=900, mtry=5, nodesize=9, sample_fraction=0.6)
 * are the R-side tuned hyperparameters; mtry=5 sits at the edge of the
 * narrower 2-5 range R used for this stage specifically (a larger, more
 * correlated predictor set than the other stages). FEATURE_BANDS matches
 * R/06_AGC_final/rf_model_AGC3.R's predictor set exactly -- unlike the
 * earlier stages, the final AGC stage deliberately includes every
 * upstream proxy as a predictor.
 *
 * PERSISTENCE_MASK below uses seagrass_persistence_mask_byYearCount_
 * 2017_2025, which matches 02_extent_2017-2025.js's own export exactly
 * (both say "2025") -- unlike the "Seagrass_MaxProb..._2017_2024" mask
 * referenced in 03_MORPHO/04_COVER/05_AGB/06_cIndex, which does NOT
 * match that same script's "_2017_2025" export. Since this script's
 * reference lines up cleanly, it further supports the earlier scripts'
 * "2024" references being the stale/outdated ones, not
 * 02_extent_2017-2025.js's naming.
 *
 * FIX APPLIED: same as the other stage scripts -- the display-only layer
 * used pred.updateMask(PERSISTENCE_MASK) with PERSISTENCE_MASK as a
 * plain string (an asset path), not an ee.Image. Wrapped in ee.Image(...)
 * here too. Only affects the optional map-preview layer, not the actual
 * exported RF_AGC_<year> assets (which use the unmasked pred directly).
 *
 * ASSET PORTABILITY: this script originally referenced the author's own
 * private GEE assets, all under one project/folder root -- replace
 * ASSET_ROOT below with your own GEE project and asset folder. The
 * Allen Coral Atlas distToLand layer is a separate, externally-hosted
 * shared asset (not the author's own) -- using it requires your own
 * access to that specific Coral Atlas asset, not just a path swap.
 ************************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var RESPONSE = 'AGC_pred';
var K_FOLDS  = 5;
var SCALE    = 10;
var SEED     = 42;

var TRAIN_YEARS = [2018, 2019, 2020, 2021, 2022, 2023];
var APPLY_YEARS = ee.List.sequence(2017, 2024);

// Asset prefixes
var TRAIN_PREFIX   = ASSET_ROOT + '/output/training_embedding_depth3_';

var PA_PROB_PREFIX = ASSET_ROOT + '/output/RF_probability04022026_';
var MORPH_PREFIX   = ASSET_ROOT + '/output/MORPH3_probs04022026_';

var SPC_PREFIX     = ASSET_ROOT + '/output/RF_tSPC04022026_';
var AGB_PREFIX     = ASSET_ROOT + '/output/RF_AGB04022026_';
var CINDEX_PREFIX  = ASSET_ROOT + '/output/RF_cIndex04022026_';

var EXPORT_PREFIX  = ASSET_ROOT + '/output/RF_AGC04022026_';

// Persistence mask (see note above on why this one's naming lines up
// cleanly, unlike the earlier stages' MaxProb mask reference)
var PERSISTENCE_MASK =
  ASSET_ROOT + '/output/seagrass_persistence_mask_byYearCount_2017_2025';

// RF Hyperparameters
var RF_PARAMS = {
  numberOfTrees: 900,
  variablesPerSplit: 5,
  minLeafPopulation: 9,
  bagFraction: 0.6,
  seed: SEED
};


// ==========================================================
// LOAD STATIC LAYERS
// ==========================================================
var AOI = ee.FeatureCollection(
  ASSET_ROOT + '/InputArea/R2_case_study'
).geometry();

var depth = ee.Image(
  ASSET_ROOT + '/InputArea/raster_INA_ACA_Bathymetry'
).rename('depth').multiply(-1);

// Distance to land. EXTERNAL ASSET: hosted by the Allen Coral Atlas
// team, not the author's own -- requires your own access to this
// specific shared asset, not just a path change.
var distToLand = ee.Image(
  'projects/coral_atlas/global_datasets/osm_distToLand_indo'
).rename('distToLand');

var STATIC_BANDS = ['depth', 'distToLand'];


// ==========================================================
// EMBEDDING LAYERS
// ==========================================================
var EMB_COLL = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';

var GSE_BANDS = ee.List.sequence(0, 63).map(function(i){
  return ee.String('A').cat(ee.Number(i).format('%02d'));
});


// ==========================================================
// MORPHOLOGY PROBABILITY BANDS (3-class Morph3)
// ==========================================================
var MORPH_BANDS = [
  'P_mixed_long',
  'P_mixed_short_plus_mono_short',
  'P_mono_Ea'
];


// ==========================================================
// FINAL FEATURE BAND LIST
// ==========================================================
var FEATURE_BANDS = GSE_BANDS
  .cat(STATIC_BANDS)
  .add('PA_prob')
  .cat(MORPH_BANDS)
  .add('tSPC_pred')
  .add('tAGB_pred')
  .add('carbon_index');

print('Predictor bands used in AGC model:', FEATURE_BANDS);


// ==========================================================
// LOAD DYNAMIC LAYERS
// ==========================================================
function loadPAprob(year) {
  return ee.Image(PA_PROB_PREFIX + year).rename('PA_prob');
}

function loadMorphProbs(year) {
  return ee.Image(MORPH_PREFIX + year).select(MORPH_BANDS);
}

function loadSPC(year) {
  return ee.Image(SPC_PREFIX + year).rename('tSPC_pred');
}

function loadAGB(year) {
  return ee.Image(AGB_PREFIX + year).rename('tAGB_pred');
}

function loadCarbonIndex(year) {
  return ee.Image(CINDEX_PREFIX + year).rename('carbon_index');
}


// ==========================================================
// CLEAN NA VALUES (AGC label)
// ==========================================================
function cleanAGC(f) {
  var val = f.get(RESPONSE);

  var safeVal = ee.Algorithms.If(
    ee.Algorithms.IsEqual(val, 'NA'),
    -9999,
    val
  );

  return f.set(RESPONSE, ee.Number.parse(safeVal));
}


// ==========================================================
// ENRICH TRAINING DATA WITH PROXY STACKS
// ==========================================================
function enrichTraining(year){

  var pts = ee.FeatureCollection(TRAIN_PREFIX + year);

  var pa     = loadPAprob(year);
  var morph  = loadMorphProbs(year);
  var spc    = loadSPC(year);
  var agb    = loadAGB(year);
  var cindex = loadCarbonIndex(year);

  // Full predictor stack
  var enriched = pa
    .addBands(morph)
    .addBands(spc)
    .addBands(agb)
    .addBands(cindex);

  // Sample at points
  var samples = enriched.sampleRegions({
    collection: pts,
    properties: pts.first().propertyNames(),
    scale: SCALE,
    geometries: true
  });

  samples = samples.map(cleanAGC);

  return samples.filter(ee.Filter.neq(RESPONSE, -9999));
}


// ==========================================================
// BUILD FULL TRAINING SET
// ==========================================================
var training = ee.FeatureCollection(TRAIN_YEARS.map(enrichTraining)).flatten();

training = training
  .filter(ee.Filter.gt(RESPONSE, 0))
  .filter(ee.Filter.lte(RESPONSE, 2000));

print('Training size:', training.size());
print('Preview:', training.limit(5));


// ==========================================================
// FINAL MODEL TRAINING (Regression RF)
// ==========================================================
var rf = ee.Classifier.smileRandomForest(RF_PARAMS)
  .setOutputMode('REGRESSION')
  .train({
    features: training,
    classProperty: RESPONSE,
    inputProperties: FEATURE_BANDS
  });

print('Final RF AGC model trained');


// ==========================================================
// VARIABLE IMPORTANCE
// ==========================================================
var importance = ee.Dictionary(rf.explain().get('importance'));

var importanceFC = ee.FeatureCollection(
  importance.keys().map(function (k) {
    return ee.Feature(null, {
      feature: k,
      importance: importance.get(k)
    });
  })
).sort('importance', false);

print('Feature importance (Top 15):', importanceFC.limit(15));


// ==========================================================
// APPLY MODEL PER YEAR
// ==========================================================
function buildStack(year){

  var start = ee.Date.fromYMD(year, 1, 1);

  var gse = ee.ImageCollection(EMB_COLL)
    .filterDate(start, start.advance(1, 'year'))
    .filterBounds(AOI)
    .mosaic();

  var pa     = loadPAprob(year);
  var morph  = loadMorphProbs(year);
  var spc    = loadSPC(year);
  var agb    = loadAGB(year);
  var cindex = loadCarbonIndex(year);

  return gse
    .select(GSE_BANDS)
    .addBands(depth)
    .addBands(distToLand)
    .addBands(pa)
    .addBands(morph)
    .addBands(spc)
    .addBands(agb)
    .addBands(cindex);
}


// ==========================================================
// EXPORT AGC PREDICTIONS
// ==========================================================
APPLY_YEARS.getInfo().forEach(function(y){

  print('Applying AGC model for year:', y);

  var stack = buildStack(y);

  var pred = stack.classify(rf).rename('AGC_pred');

  var display = pred.updateMask(ee.Image(PERSISTENCE_MASK));

  Map.addLayer(display, {
    min: 0,
    max: 100,
    palette: ['ffffff', 'c6dbef', '6baed6', '08306b']
  }, 'AGC ' + y, false);

  Export.image.toAsset({
    image: pred.clip(AOI),
    description: 'RF_AGC_' + y,
    assetId: EXPORT_PREFIX + y,
    region: AOI,
    scale: SCALE,
    maxPixels: 1e13
  });

  print('Export queued AGC', y);
});

// ==========================================================
// DISPLAY TRAINING POINTS ON MAP
// ==========================================================

// Color palette for AGC
var agc_palette = ['ffffff', 'fee0d2', 'fc9272', 'de2d26'];

// Buffer points based on AGC value
var bufferedPoints = training.map(function (f) {
  var value = ee.Number(f.get(RESPONSE));

  var bufferSize = value
    .multiply(0.02)
    .clamp(1, 5);

  return f.buffer(bufferSize).copyProperties(f);
});

// Rasterize for bubble effect
var bubbleRaster = ee.Image().float().paint({
  featureCollection: bufferedPoints,
  color: RESPONSE
});

var visBubbles = {
  min: 0,
  max: 300,
  palette: agc_palette,
  opacity: 0.7
};

Map.addLayer(bubbleRaster, visBubbles, 'AGC Training Bubbles', false);

// Big red training points
var styledPoints = training.map(function (f) {
  return f.set('style', {
    color: 'ff0000',
    pointSize: 3,
    width: 0
  });
});

Map.addLayer(
  styledPoints.style({styleProperty: 'style'}),
  {},
  'AGC Training Points',
  false
);