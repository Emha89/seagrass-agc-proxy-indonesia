/*****************************************************************************
 * 10_agc_model.js
 * Random Forest model of above-ground carbon (AGC) in GEE, trained on the field
 * point estimate (AGC_pred) with the settings of Table S3 of the paper (900
 * trees, mtry 5, minimum leaf 9, bag fraction 0.6). Predictors: the 64 GSE
 * bands, depth, distance to land and the proxy layers of 04 and 06 to 09
 * (presence probability, morphology probabilities, percent cover, AGB and the
 * predicted carbon index). The model is applied to every year from 2017 to 2024.
 *
 * Inputs  (ASSET_ROOT)  training_gse_<year>, presence_prob_<year>,
 *                       morphology_prob_<year>, cover_<year>, agb_<year>,
 *                       carbon_index_<year>, study_regions, bathymetry,
 *                       seagrass_extent_2017_2024 (map preview only)
 * Outputs (ASSET_ROOT)  agc_<year>
 *
 * Set ASSET_ROOT below to your own GEE asset folder. All paths are built
 * from it, except the shared Allen Coral Atlas distance-to-land layer,
 * which needs its own access.
 ****************************************************************************/

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
var TRAIN_PREFIX   = ASSET_ROOT + '/training_gse_';

var PA_PROB_PREFIX = ASSET_ROOT + '/presence_prob_';
var MORPH_PREFIX   = ASSET_ROOT + '/morphology_prob_';

var SPC_PREFIX     = ASSET_ROOT + '/cover_';
var AGB_PREFIX     = ASSET_ROOT + '/agb_';
var CINDEX_PREFIX  = ASSET_ROOT + '/carbon_index_';

var EXPORT_PREFIX  = ASSET_ROOT + '/agc_';

// Mapped seagrass extent (05_seagrass_extent.js), used only for the map preview
var PERSISTENCE_MASK =
  ASSET_ROOT + '/seagrass_extent_2017_2024';

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
  ASSET_ROOT + '/study_regions'
).geometry();

var depth = ee.Image(
  ASSET_ROOT + '/bathymetry'
).rename('depth').multiply(-1);

// Distance to land (shared Allen Coral Atlas asset, needs its own access)
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