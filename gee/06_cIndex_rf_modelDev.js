/************************************************************
 * 06_cIndex_rf_modelDev.js
 * RF regression for carbon index (cIndex).
 * Predictors: GSE (A00-A63) + depth + distToLand ONLY
 * Response: carbon_ind (0-1)
 * Training years: 2018-2023
 * Apply years: 2017-2024
 * Includes K-fold cross validation
 *
 * PIPELINE POSITION: reads training_embedding_depth3_<year> (from
 * 01_trainingData_prep.js). RF_PARAMS below (ntree=500, mtry=6,
 * nodesize=9, sample_fraction=0.6) are the R-side tuned hyperparameters
 * for this stage.
 *
 * CONFIRMED: the response property here, carbon_ind, is the R-side
 * carbon index model's own predictions (carbon_index_pred from
 * R/05_carbon_index/rf_applyModel_cIndex.R) at each training point's
 * location -- not a direct field measurement, and not the raw
 * species-weighted formula output from
 * carbon_scaling_join_SP_finalDell.R directly. Since carbon index can't
 * be measured directly from satellite data alone (the species-weighted
 * formula needs species composition survey data, which isn't available
 * everywhere), this GEE script trains a GEE-native model to approximate
 * R's carbon index model using only GSE + depth + distToLand, so it can
 * be deployed across the full raster extent where species data doesn't
 * exist. In effect, R's predictions become this model's training
 * targets -- a similar hyperparameter hand-off pattern to the other
 * stages, but here the training labels themselves (not just the tuned
 * hyperparameters) come from R.
 *
 * NAMING NOTE: PERSISTENCE_MASK below points at
 * Seagrass_MaxProb_0_8_2017_2024_mask -- the same "2024" suffix also
 * expected by 03_MORPHO, 04_COVER, and 05_AGB. Four independent scripts
 * now expecting "2024" is strong evidence the "2025" suffix written into
 * 02_extent_2017-2025.js's own EXPORT_PREFIX is the actual mistake --
 * worth revisiting that script, not changed here without confirmation.
 *
 * FIX APPLIED: same as 04_COVER/05_AGB -- the display-only layer used
 * pred.updateMask(PERSISTENCE_MASK) with PERSISTENCE_MASK as a plain
 * string (an asset path), not an ee.Image. Wrapped in ee.Image(...)
 * here too. Only affects the optional map-preview layer, not the actual
 * exported RF_cIndex_<year> assets (which use the unmasked pred
 * directly).
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

var RESPONSE = 'carbon_ind';
var K_FOLDS = 5;
var SCALE = 10;
var SEED = 42;

var TRAIN_YEARS = [2018, 2019, 2020, 2021, 2022, 2023];
var APPLY_YEARS = ee.List.sequence(2017, 2024);

var TRAIN_PREFIX = ASSET_ROOT + '/output/training_embedding_depth3_';

var EXPORT_PREFIX = ASSET_ROOT + '/output/RF_cIndex04022026_';

// See NAMING NOTE above about the "2024" vs "2025" suffix
var PERSISTENCE_MASK = ASSET_ROOT + '/output/Seagrass_MaxProb_0_8_2017_2024_mask';

// RF hyperparameters
var RF_PARAMS = {
  numberOfTrees: 500,
  variablesPerSplit: 6,
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
// EMBEDDING BANDS
// ==========================================================
var EMB_COLL = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';

var GSE_BANDS = ee.List.sequence(0, 63).map(function (i) {
  return ee.String('A').cat(ee.Number(i).format('%02d'));
});

// Final predictor list (base only)
var FEATURE_BANDS = GSE_BANDS.cat(STATIC_BANDS);

print('Predictor bands (Carbon Index):', FEATURE_BANDS);

// ==========================================================
// BUILD TRAINING DATASET (2018-2023)
// ==========================================================

function loadTraining(year) {
  var fc = ee.FeatureCollection(TRAIN_PREFIX + year);
  return fc;
}

var training = ee.FeatureCollection(
  TRAIN_YEARS.map(loadTraining)
).flatten();

// ==========================================================
// FILTER TRAINING DATA
// ==========================================================

training = training
  .filter(ee.Filter.gt(RESPONSE, 0))
  .filter(ee.Filter.lte(RESPONSE, 1));

print('Filtered training sample size:', training.size());
print('Training preview:', training.limit(5));

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
// TRAIN FULL MODEL ON ALL TRAINING DATA
// ==========================================================

var rf = ee.Classifier.smileRandomForest(RF_PARAMS)
  .setOutputMode('REGRESSION')
  .train({
    features: training,
    classProperty: RESPONSE,
    inputProperties: FEATURE_BANDS
  });

print('Carbon Index model trained');

// ==========================================================
// APPLY MODEL PER YEAR (2017-2024)
// ==========================================================

function buildStack(year) {

  var start = ee.Date.fromYMD(year, 1, 1);

  var img = ee.ImageCollection(EMB_COLL)
    .filterDate(start, start.advance(1, 'year'))
    .filterBounds(AOI)
    .mosaic();

  return img
    .select(GSE_BANDS)
    .addBands(depth)
    .addBands(distToLand);
}

APPLY_YEARS.getInfo().forEach(function (y) {

  var stack = buildStack(y);

  var pred = stack.classify(rf).rename('cIndex_pred');

  // Display mask only
  var display = pred.updateMask(ee.Image(PERSISTENCE_MASK));

  Map.addLayer(display, {
    min: 0,
    max: 0.2,
    palette: ['ffffff', 'fcae91', 'fb6a4a', 'cb181d']
  }, 'cIndex ' + y, false);

  // Export full raster
  Export.image.toAsset({
    image: pred.clip(AOI),
    description: 'RF_cIndex_' + y,
    assetId: EXPORT_PREFIX + y,
    region: AOI,
    scale: SCALE,
    maxPixels: 1e13
  });

  print('Export queued for', y);
});