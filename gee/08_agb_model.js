/*****************************************************************************
 * 08_agb_model.js
 * Random Forest model of above-ground biomass in GEE, trained on the field
 * point estimate (AGB_pred) with the settings of the R model (Table S3 of the
 * paper: 300 trees, mtry 8, minimum leaf 9, bag fraction 0.6). Predictors: the
 * 64 GSE bands, depth and distance to land.
 *
 * Inputs  (ASSET_ROOT)  training_gse_<year>, study_regions, bathymetry,
 *                       seagrass_extent_2017_2024 (map preview only)
 * Outputs (ASSET_ROOT)  agb_<year>
 *
 * Set ASSET_ROOT below to your own GEE asset folder. All paths are built
 * from it, except the shared Allen Coral Atlas distance-to-land layer,
 * which needs its own access.
 ****************************************************************************/

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

var TRAIN_PREFIX = ASSET_ROOT + '/training_gse_';

var EXPORT_PREFIX = ASSET_ROOT + '/agb_';

// Mapped seagrass extent (05_seagrass_extent.js), used only for the map preview
var PERSISTENCE_MASK = ASSET_ROOT + '/seagrass_extent_2017_2024';

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