/*****************************************************************************
 * 07_cover_model.js
 * Random Forest model of seagrass percent cover (tSPC) in GEE (Table S3 of the
 * paper: 300 trees, mtry 10, minimum leaf 9, bag fraction 0.7). Predictors: the
 * 64 GSE bands, depth, distance to land and the morphology probabilities.
 * For mapping, the model is trained only on samples with seagrass present
 * (tSPC > 0), because most samples with zero cover are non-seagrass points
 * from the Allen Coral Atlas.
 *
 * Inputs  (ASSET_ROOT)  training_gse_<year>, morphology_prob_<year>,
 *                       study_regions, bathymetry,
 *                       seagrass_extent_2017_2024 (map preview only)
 * Outputs (ASSET_ROOT)  cover_<year>
 *
 * Set ASSET_ROOT below to your own GEE asset folder. All paths are built
 * from it, except the shared Allen Coral Atlas distance-to-land layer,
 * which needs its own access.
 ****************************************************************************/

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

var TRAIN_PREFIX = ASSET_ROOT + '/training_gse_';

// Mapped seagrass extent (05_seagrass_extent.js), used only for the map preview
var PERSISTENCE_MASK = ASSET_ROOT + '/seagrass_extent_2017_2024';

var EXPORT_PREFIX = ASSET_ROOT + '/cover_';

// ==========================================================
// MORPHOLOGY PROBABILITY BANDS (06_morphology_model.js)
// ==========================================================
var MORPH_BANDS = [
  'P_mixed_long',
  'P_mixed_short_plus_mono_short',
  'P_mono_Ea'
];

// Morphology probability asset prefix (from 06_morphology_model.js)
var MORPH_PROB_PREFIX = ASSET_ROOT + '/morphology_prob_';

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