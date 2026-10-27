/*****************************************************************************
 * 09_carbon_index_model.js
 * Random Forest model of the carbon index (0-1) in GEE (Table S3 of the paper:
 * 500 trees, mtry 6, minimum leaf 9, bag fraction 0.6). Predictors: the 64 GSE
 * bands, depth and distance to land. The response column is 'carbon_ind'
 * (the carbon index, name shortened by the shapefile upload).
 *
 * Inputs  (ASSET_ROOT)  training_gse_<year>, study_regions, bathymetry,
 *                       seagrass_extent_2017_2024 (map preview only)
 * Outputs (ASSET_ROOT)  carbon_index_<year>
 *
 * Set ASSET_ROOT below to your own GEE asset folder. All paths are built
 * from it, except the shared Allen Coral Atlas distance-to-land layer,
 * which needs its own access.
 ****************************************************************************/

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

var TRAIN_PREFIX = ASSET_ROOT + '/training_gse_';

var EXPORT_PREFIX = ASSET_ROOT + '/carbon_index_';

// Mapped seagrass extent (05_seagrass_extent.js), used only for the map preview
var PERSISTENCE_MASK = ASSET_ROOT + '/seagrass_extent_2017_2024';

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
  ASSET_ROOT + '/study_regions'
).geometry();

var depth = ee.Image(
  ASSET_ROOT + '/bathymetry'
)
  .rename('depth')
  .multiply(-1);

// Distance to land (shared Allen Coral Atlas asset, needs its own access)
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