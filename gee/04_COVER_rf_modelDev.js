/************************************************************
 * 04_COVER_rf_modelDev.js
 * RF regression for tSPC (seagrass percent cover).
 *
 * Predictors:
 *   - GSE embedding (A00-A63)
 *   - depth + distToLand
 *   - morphology probability bands (Stage-A output)
 * Response: tSPC (0-100)
 * Training years: 2018-2023
 * Apply years:  2017-2024
 * Includes K-fold cross validation
 *
 * PIPELINE POSITION: reads training_embedding_depth3_<year> (from
 * 01_trainingData_prep.js) and MORPH3_probs04022026_<year> (from
 * 03_MORPHO_rf_modelDev.js -- prefix matches that script's export
 * exactly, date stamp included). RF_PARAMS below (ntree=300, mtry=10,
 * nodesize=9, sample_fraction=0.7) are presumed to be the R-side tuned
 * hyperparameters for the SPC model. FEATURE_BANDS (GSE + depth +
 * distToLand + morphology probabilities) matches the R-side SPC
 * predictor set exactly (R/03_SPC/rf_model_combined_SPC.R).
 *
 * NAMING NOTE: PERSISTENCE_MASK below points at
 * Seagrass_MaxProb_0_8_2017_2024_mask -- the same "2024" suffix also
 * expected by 03_MORPHO_rf_modelDev.js. Two independent scripts now
 * expecting "2024" is fairly strong evidence that the "2025" suffix
 * actually written into 02_extent_2017-2025.js's own EXPORT_PREFIX may
 * be the mistake, rather than these two references being wrong -- worth
 * revisiting that script's naming, but not changed here without
 * confirmation.
 *
 * FIX APPLIED: the display-only layer used
 * pred.updateMask(PERSISTENCE_MASK) with PERSISTENCE_MASK as a plain
 * string (an asset path), not an ee.Image -- unlike
 * 03_MORPHO_rf_modelDev.js, which correctly wraps it in ee.Image(...)
 * first. Wrapped it here the same way. This only affects the optional
 * map-preview layer, not the actual exported RF_tSPC_<year> assets
 * (which use the unmasked pred directly).
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

var RESPONSE = 'tSPC';
var K_FOLDS = 5;
var SCALE = 10;
var SEED = 42;

var TRAIN_YEARS = [2018, 2019, 2020, 2021, 2022, 2023];
var APPLY_YEARS = ee.List.sequence(2017, 2024);

var TRAIN_PREFIX = ASSET_ROOT + '/output/training_embedding_depth3_';

// See NAMING NOTE above about the "2024" vs "2025" suffix
var PERSISTENCE_MASK = ASSET_ROOT + '/output/Seagrass_MaxProb_0_8_2017_2024_mask';

var EXPORT_PREFIX = ASSET_ROOT + '/output/RF_tSPC04022026_';

// ==========================================================
// MORPHOLOGY PROBABILITY BANDS (Stage-A Output)
// ==========================================================
var MORPH_BANDS = [
  'P_mixed_long',
  'P_mixed_short_plus_mono_short',
  'P_mono_Ea'
];

// Morphology probability asset prefix (from 03_MORPHO_rf_modelDev.js)
var MORPH_PROB_PREFIX = ASSET_ROOT + '/output/MORPH3_probs04022026_';

// RF hyperparameters
var RF_PARAMS = {
  numberOfTrees: 300,
  variablesPerSplit: 10,
  minLeafPopulation: 9,
  bagFraction: 0.7,
  seed: SEED
};

// ==========================================================
// LOAD STATIC LAYERS (depth, distToLand)
// ==========================================================
var AOI = ee.FeatureCollection(
  ASSET_ROOT + '/InputArea/R2_case_study'
).geometry();

var depth = ee.Image(
  ASSET_ROOT + '/InputArea/raster_INA_ACA_Bathymetry'
)
  .rename('depth')
  .multiply(-1);

// Distance to land. EXTERNAL ASSET: hosted by the Allen Coral Atlas
// team, not the author's own -- requires your own access to this
// specific shared asset, not just a path change.
var distToLand = ee.Image(
  'projects/coral_atlas/global_datasets/osm_distToLand_indo'
).rename('distToLand');

var STATIC_BANDS = ['depth', 'distToLand'];

// ==========================================================
// GSE embedding bands
// ==========================================================
var EMB_COLL = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';

var GSE_BANDS = ee.List.sequence(0, 63).map(function(i){
  return ee.String('A').cat(ee.Number(i).format('%02d'));
});

// Final predictors
var FEATURE_BANDS =
  GSE_BANDS
    .cat(STATIC_BANDS)
    .cat(MORPH_BANDS);

print('Predictor bands:', FEATURE_BANDS);

// ==========================================================
// Load Morphology Probability Image per year
// ==========================================================
function loadMorphProb(year) {
  return ee.Image(MORPH_PROB_PREFIX + year)
    .select(MORPH_BANDS);
}

// ==========================================================
// Add Morphology Probabilities to training points
// ==========================================================
function enrichTraining(year){

  var pts   = ee.FeatureCollection(TRAIN_PREFIX + year);
  var morph = loadMorphProb(year);

  return morph.sampleRegions({
    collection: pts,
    properties: pts.first().propertyNames(),
    scale: SCALE,
    geometries: true
  });
}

// Build training dataset 2018-2023 with morph probs
var training =
  ee.FeatureCollection(TRAIN_YEARS.map(enrichTraining)).flatten();

// Filter valid tSPC values
training = training
  .filter(ee.Filter.gt(RESPONSE, 0))
  .filter(ee.Filter.lte(RESPONSE, 100));

print('Training sample size:', training.size());
print('Training preview:', training.limit(5));

// ==========================================================
// K-FOLD VALIDATION
// ==========================================================
var trainingWithFold =
  training.randomColumn('rand', SEED).map(function(f){
    return f.set('fold',
      ee.Number(f.get('rand')).multiply(K_FOLDS).int()
    );
  });

var foldIndices = ee.List.sequence(0, K_FOLDS - 1);

var cv = ee.FeatureCollection(
  foldIndices.map(function(k){

    k = ee.Number(k);

    var trainFold =
      trainingWithFold.filter(ee.Filter.neq('fold', k));

    var validFold =
      trainingWithFold.filter(ee.Filter.eq('fold', k));

    var model = ee.Classifier.smileRandomForest(RF_PARAMS)
      .setOutputMode('REGRESSION')
      .train({
        features: trainFold,
        classProperty: RESPONSE,
        inputProperties: FEATURE_BANDS
      });

    return validFold.classify(model, 'pred').map(function(f){

      var y    = ee.Number(f.get(RESPONSE));
      var yhat = ee.Number(f.get('pred'));
      var diff = y.subtract(yhat);

      return f.set({
        fold: k,
        residual: diff,
        abs_res: diff.abs(),
        sq_res: diff.pow(2)
      });

    });
  })
).flatten();

var mae = cv.aggregate_mean('abs_res');
var rmse = ee.Number(cv.aggregate_mean('sq_res')).sqrt();

print('CV MAE:', mae);
print('CV RMSE:', rmse);

// R2
var yMean = ee.Number(cv.aggregate_mean(RESPONSE));

var cv2 = cv.map(function(f){
  var y = ee.Number(f.get(RESPONSE));
  var d = y.subtract(yMean);
  return f.set('sq_tot', d.pow(2));
});

var sse = ee.Number(cv2.aggregate_sum('sq_res'));
var sst = ee.Number(cv2.aggregate_sum('sq_tot'));

var r2 = ee.Number(1).subtract(sse.divide(sst));
print('CV R2:', r2);

// ==========================================================
// TRAIN FULL MODEL ON ALL DATA
// ==========================================================
var rf = ee.Classifier.smileRandomForest(RF_PARAMS)
  .setOutputMode('REGRESSION')
  .train({
    features: training,
    classProperty: RESPONSE,
    inputProperties: FEATURE_BANDS
  });

print('Model trained');

// ==========================================================
// FUNCTION TO BUILD STACK & APPLY MODEL PER YEAR
// ==========================================================
function buildStack(year){

  var start = ee.Date.fromYMD(year, 1, 1);

  var img = ee.ImageCollection(EMB_COLL)
    .filterDate(start, start.advance(1, 'year'))
    .filterBounds(AOI)
    .mosaic();

  var morphProb = loadMorphProb(year);

  return img
    .select(GSE_BANDS)
    .addBands(depth)
    .addBands(distToLand)
    .addBands(morphProb);
}

// ==========================================================
// Apply per year 2017-2024
// ==========================================================
APPLY_YEARS.getInfo().forEach(function(y){

  var stack = buildStack(y);

  var pred = stack.classify(rf).rename('tSPC_pred');

  // ===== Display: Apply mask for map visualization only =====
  var display = pred.updateMask(ee.Image(PERSISTENCE_MASK));

  Map.addLayer(display, {
    min: 35,
    max: 65,
    palette: ['ffffff', 'd9f0a3', '78c679', '238443']
  }, 'tSPC ' + y, false);

  Export.image.toAsset({
    image: pred.clip(AOI),
    description: 'RF_tSPC_' + y,
    assetId: EXPORT_PREFIX + y,
    region: AOI,
    scale: SCALE,
    maxPixels: 1e13
  });

  print('Export queued for', y);
});

// ==========================================================
// DISPLAY TRAINING POINTS ON MAP (unchanged)
// ==========================================================

// Define color palette for tSPC (0-100)
var tSPC_palette = ['ffffff', 'fef0d9', 'fdcc8a', 'fc8d59', 'd7301f'];

// Buffer visualization
var bufferedPoints = training.map(function (f) {
  var value = ee.Number(f.get(RESPONSE));
  var bufferSize = value.clamp(1, 5);
  return f.buffer(bufferSize).copyProperties(f);
});

var bubbleRaster = ee.Image().float().paint({
  featureCollection: bufferedPoints,
  color: RESPONSE
});

Map.addLayer(bubbleRaster, {
  min: 0,
  max: 100,
  palette: tSPC_palette,
  opacity: 0.7
}, 'Training Bubbles (size ~ tSPC)', false);

// Vector style
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
  'Training Points (Big Red)',
  false
);