/************************************************************
 * 05_AGB_rf_modelDev.js
 * RF regression for AGB (above-ground biomass).
 * Predictors: GSE (A00-A63) + depth + distToLand ONLY. No proxy layers
 * used (deliberately excludes tSPC_pred / morphology probabilities, to
 * avoid predictor redundancy across proxy stages -- same principle as
 * the R-side pipeline).
 *
 * PIPELINE POSITION: reads training_embedding_depth3_<year> (from
 * 01_trainingData_prep.js). RF_PARAMS below (ntree=300, mtry=8,
 * nodesize=9, sample_fraction=0.6) are the R-side tuned hyperparameters
 * for the AGB model (per the original script's own comment). Predictor
 * set matches R/04_AGB/rf_model_AGB2.R exactly.
 *
 * NAMING NOTE: PERSISTENCE_MASK below points at
 * Seagrass_MaxProb_0_8_2017_2024_mask -- the same "2024" suffix also
 * expected by 03_MORPHO_rf_modelDev.js and 04_COVER_rf_modelDev.js.
 * Three independent scripts now expecting "2024" is strong evidence the
 * "2025" suffix written into 02_extent_2017-2025.js's own EXPORT_PREFIX
 * is the actual mistake -- worth revisiting that script, not changed
 * here without confirmation.
 *
 * FIX APPLIED: same as 04_COVER_rf_modelDev.js -- the display-only layer
 * used pred.updateMask(PERSISTENCE_MASK) with PERSISTENCE_MASK as a
 * plain string (an asset path), not an ee.Image. Wrapped in ee.Image(...)
 * here too. Only affects the optional map-preview layer, not the actual
 * exported RF_AGB_<year> assets (which use the unmasked pred directly).
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

var RESPONSE = 'AGB_pred';
var K_FOLDS = 5;
var SCALE = 10;
var SEED = 42;

var TRAIN_YEARS = [2018, 2019, 2020, 2021, 2022, 2023];
var APPLY_YEARS = ee.List.sequence(2017, 2024);

var TRAIN_PREFIX = ASSET_ROOT + '/output/training_embedding_depth3_';

var EXPORT_PREFIX = ASSET_ROOT + '/output/RF_AGB04022026_';

// See NAMING NOTE above about the "2024" vs "2025" suffix
var PERSISTENCE_MASK = ASSET_ROOT + '/output/Seagrass_MaxProb_0_8_2017_2024_mask';

// RF Hyperparameters (from R tuning)
var RF_PARAMS = {
  numberOfTrees: 300,
  variablesPerSplit: 8,
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
// EMBEDDING LAYERS
// ==========================================================
var EMB_COLL = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';

var GSE_BANDS = ee.List.sequence(0, 63).map(function (i) {
  return ee.String('A').cat(ee.Number(i).format('%02d'));
});

// Final predictor list (base only)
var FEATURE_BANDS = GSE_BANDS.cat(STATIC_BANDS);

print('Predictor bands (AGB):', FEATURE_BANDS);

// ==========================================================
// TRAINING DATA PREP
// ==========================================================

// Convert string "NA" into numeric safe value
function cleanAGB(f) {
  var val = f.get(RESPONSE);

  var safeVal = ee.Algorithms.If(
    ee.Algorithms.IsEqual(val, 'NA'),
    -9999,
    val
  );

  return f.set(RESPONSE, ee.Number.parse(safeVal));
}

// Load training samples per year (already includes predictors)
function enrichTraining(year) {

  var pts = ee.FeatureCollection(TRAIN_PREFIX + year);

  // No proxy enrichment needed
  var samples = pts.map(cleanAGB);

  // Remove invalid rows
  return samples.filter(ee.Filter.neq(RESPONSE, -9999));
}

// Merge all training years
var training = ee.FeatureCollection(
  TRAIN_YEARS.map(enrichTraining)
).flatten();

// Filter valid AGB values
training = training.filter(ee.Filter.gt(RESPONSE, 0));

print('Training size:', training.size());
print('Preview:', training.limit(5));

// ==========================================================
// K-FOLD VALIDATION
// ==========================================================
var trainingWithFold = training.randomColumn('rand', SEED).map(function (f) {
  return f.set('fold', ee.Number(f.get('rand')).multiply(K_FOLDS).int());
});

var foldIndices = ee.List.sequence(0, K_FOLDS - 1);

var cv = ee.FeatureCollection(
  foldIndices.map(function (k) {

    k = ee.Number(k);

    var trainFold = trainingWithFold.filter(ee.Filter.neq('fold', k));
    var validFold = trainingWithFold.filter(ee.Filter.eq('fold', k));

    var model = ee.Classifier.smileRandomForest(RF_PARAMS)
      .setOutputMode('REGRESSION')
      .train({
        features: trainFold,
        classProperty: RESPONSE,
        inputProperties: FEATURE_BANDS
      });

    return validFold.classify(model, 'pred').map(function (f) {

      var y = ee.Number(f.get(RESPONSE));
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

// Metrics
var mae = cv.aggregate_mean('abs_res');
var rmse = ee.Number(cv.aggregate_mean('sq_res')).sqrt();

print('CV MAE:', mae);
print('CV RMSE:', rmse);

// R2
var yMean = ee.Number(cv.aggregate_mean(RESPONSE));
var cv2 = cv.map(function (f) {
  var y = ee.Number(f.get(RESPONSE));
  return f.set('sq_tot', y.subtract(yMean).pow(2));
});

var sse = ee.Number(cv2.aggregate_sum('sq_res'));
var sst = ee.Number(cv2.aggregate_sum('sq_tot'));
var r2 = ee.Number(1).subtract(sse.divide(sst));

print('CV R2:', r2);

// ==========================================================
// FINAL MODEL TRAINING
// ==========================================================
var rf = ee.Classifier.smileRandomForest(RF_PARAMS)
  .setOutputMode('REGRESSION')
  .train({
    features: training,
    classProperty: RESPONSE,
    inputProperties: FEATURE_BANDS
  });

print('Final RF AGB model trained');

// ==========================================================
// APPLY MODEL PER YEAR
// ==========================================================
function buildStack(year) {

  var start = ee.Date.fromYMD(year, 1, 1);

  var img = ee.ImageCollection(EMB_COLL)
    .filterDate(start, start.advance(1, 'year'))
    .filterBounds(AOI)
    .mosaic();

  // Predictor stack = GSE + static only
  return img
    .select(GSE_BANDS)
    .addBands(depth)
    .addBands(distToLand);
}

APPLY_YEARS.getInfo().forEach(function (y) {

  var stack = buildStack(y);

  var pred = stack.classify(rf).rename('AGB_pred');

  // Display only
  var display = pred.updateMask(ee.Image(PERSISTENCE_MASK));

  Map.addLayer(display, {
    min: 50,
    max: 200,
    palette: ['ffffff', 'fee6ce', 'fdae6b', 'e6550d']
  }, 'AGB ' + y, false);

  // Export full raster
  Export.image.toAsset({
    image: pred.clip(AOI),
    description: 'RF_AGB_' + y,
    assetId: EXPORT_PREFIX + y,
    region: AOI,
    scale: SCALE,
    maxPixels: 1e13
  });

  print('Export queued for', y);
});

// ==========================================================
// DISPLAY TRAINING POINTS ON MAP (AGB)
// ==========================================================

var agb_palette = [
  'ffffff', 'f7f7f7', 'd9f0a3',
  'addd8e', '78c679', '31a354', '006837'
];

// Bubble visualization
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
  max: 300,
  palette: agb_palette,
  opacity: 0.7
}, 'AGB Training Bubbles', false);

var styledPoints = training.map(function (f) {
  return f.set('style', {
    color: 'ff0000',
    pointSize: 3,
    width: 0
  });
});

Map.addLayer(
  styledPoints.style({ styleProperty: 'style' }),
  {},
  'AGB Training Points',
  false
);