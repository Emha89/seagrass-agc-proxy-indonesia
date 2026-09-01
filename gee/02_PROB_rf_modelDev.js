/************************************************************
 * 02_PROB_rf_modelDev.js
 * Trains the PA (seagrass occurrence probability) Random Forest model
 * natively in GEE (ee.Classifier.smileRandomForest), runs a 5-fold
 * cross-validation check, then applies the final model to produce
 * probability and binary presence-mask rasters for every year in the
 * apply range.
 *
 * - No .reproject() calls, to avoid pixel overflow
 * - Exports probability and mask per year
 *
 * PIPELINE POSITION: reads the training_embedding_depth3_<year> assets
 * produced by 01_trainingData_prep.js (TRAIN_ASSET_PREFIX below matches
 * that script's EXPORT_BASE exactly). RF_PARAMS below (ntree=700, mtry=9,
 * nodesize=3, sample_fraction=0.8) are presumed to be the best_params
 * output from the R-side grid search tuning (R/01_occurrence_PA/
 * rf_model_PA.R), hardcoded here for GEE deployment rather than re-tuned
 * in GEE -- confirm this matches the actual R tuning result.
 *
 * NOTE: training uses years 2018-2023 only (available_years below); the
 * model is then applied across 2017-2025 (APPLY_YEARS), including years
 * outside the training range.
 *
 * THRESHOLD NOTE (worth checking against what's reported in the paper):
 * three different operating thresholds appear across the pipeline --
 * the R-side final model uses 0.6 (rf_model_PA.R); the K-fold CV
 * evaluation below (which the printed confusion matrix/accuracy/F1
 * are based on) uses PROB_THRESHOLD = 0.7; but the actual exported
 * probability/mask rasters for 2017-2024 use APPLY_THRESHOLD = 0.5 (only
 * 2025, which has no APPLY_THRESHOLD entry, falls back to 0.7). The
 * printed CV accuracy therefore does not describe the same threshold as
 * the maps actually exported.
 *
 * Removed from this cleaned version: two inactive commented-out
 * alternatives from the original script (a 9-band reduced GSE_BANDS
 * subset, and mtry=7) -- the active configuration (all 64 GSE bands,
 * mtry=9) is unchanged.
 *
 * NAMING: asset IDs keep the original run-date stamp (04022026) exactly
 * as used in the source scripts, e.g. RF_probability04022026_<year> --
 * an earlier cleaned draft of this file dropped that stamp, but
 * 02_extent_2017-2025.js references RF_probability04022026_ by that
 * exact name, so it's kept here for that downstream script to resolve
 * correctly.
 *
 * ASSET PORTABILITY: this script originally referenced the author's own
 * private GEE assets. Below, ASSET_ROOT collects the paths that shared
 * one project/folder root into a single placeholder -- replace
 * ASSET_ROOT with your own GEE project and asset folder. The
 * Allen Coral Atlas distToLand layer is a separate, externally-hosted
 * shared asset (not the author's own) -- using it requires your own
 * access to that specific Coral Atlas asset, not just a path swap.
 ************************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var RESPONSE = 'PA';
var K_FOLDS = 5;
var SCALE = 10;
var SEED = 42;

var START_YEAR = 2017;
var END_YEAR = 2025;

var EXPORT_BASE = ASSET_ROOT + '/output';
var TRAIN_ASSET_PREFIX = ASSET_ROOT + '/output/training_embedding_depth3_';

var PROB_THRESHOLD = 0.7;
var APPLY_THRESHOLD = {
  2017: 0.5, 2018: 0.5, 2019: 0.5,
  2020: 0.5, 2021: 0.5, 2022: 0.5,
  2023: 0.5, 2024: 0.5
};

// Static layers toggle
var STATIC_LAYER_CONFIG = {
  'depth': { include: true },
  'slope': { include: false },
  'rugosity': { include: false },
  'distToLand': { include: true },
  'wave': { include: false }
};

// RF hyperparameters (see PIPELINE POSITION note above)
var RF_PARAMS = {
  numberOfTrees: 700, // ntree
  variablesPerSplit: 9, // mtry
  minLeafPopulation: 3, // nodesize
  bagFraction: 0.8, // sample_fraction
  seed: SEED
};

var APPLY_YEARS = ee.List.sequence(START_YEAR, END_YEAR);

// ==========================================================
// DEFINE PREDICTOR BANDS
// ==========================================================
var GSE_BANDS = ee.List.sequence(0, 63).map(function (i) {
  return ee.String('A').cat(ee.Number(i).format('%02d'));
});

var STATIC_BANDS = ee.List([]);
if (STATIC_LAYER_CONFIG['depth'].include) STATIC_BANDS = STATIC_BANDS.add('depth');
if (STATIC_LAYER_CONFIG['slope'].include) STATIC_BANDS = STATIC_BANDS.add('slope');
if (STATIC_LAYER_CONFIG['rugosity'].include) STATIC_BANDS = STATIC_BANDS.add('rugosity');
if (STATIC_LAYER_CONFIG['distToLand'].include) STATIC_BANDS = STATIC_BANDS.add('distToLand');
if (STATIC_LAYER_CONFIG['wave'].include) STATIC_BANDS = STATIC_BANDS.cat(['wElevation', 'wHeight', 'wPeriod']);

var FEATURE_BANDS = GSE_BANDS.cat(STATIC_BANDS);
print('Predictors:', FEATURE_BANDS);

// ==========================================================
// LOAD AOI (ALL LOCATIONS)
// ==========================================================
var AOI_FC = ee.FeatureCollection(ASSET_ROOT + '/InputArea/R2_case_study');
var AOI_GEOM = AOI_FC.geometry();
Map.addLayer(AOI_GEOM, {color: 'red'}, 'AOI (All Locations)');
Map.centerObject(AOI_GEOM, 5);

// ==========================================================
// LOAD TRAINING DATA
// ==========================================================
var training_all = ee.FeatureCollection([]);
var available_years = [2018, 2019, 2020, 2021, 2022, 2023];

available_years.forEach(function (y) {
  var fc = ee.FeatureCollection(TRAIN_ASSET_PREFIX + y);
  training_all = training_all.merge(fc);
});
print('Raw training size:', training_all.size());

// Remove nulls
function notNullFor(list) {
  return ee.Filter.and.apply(null, ee.List(list).map(function (b) {
    return ee.Filter.notNull([b]);
  }));
}
training_all = training_all.filter(notNullFor(FEATURE_BANDS.add(RESPONSE)));
print('Filtered training size:', training_all.size());

// ==========================================================
// K-FOLD VALIDATION
// ==========================================================
var withFold = training_all.randomColumn('rand', SEED).map(function (f) {
  return f.set('fold', ee.Number(f.get('rand')).multiply(K_FOLDS).int());
});
print('Fold distribution:', withFold.aggregate_histogram('fold'));

var foldIndices = ee.List.sequence(0, K_FOLDS - 1);
var cvPredictions = foldIndices.map(function (k) {
  k = ee.Number(k);
  var train = withFold.filter(ee.Filter.neq('fold', k));
  var valid = withFold.filter(ee.Filter.eq('fold', k));

  var rf = ee.Classifier.smileRandomForest(RF_PARAMS).train({
    features: train,
    classProperty: RESPONSE,
    inputProperties: FEATURE_BANDS
  });

  var prob = rf.setOutputMode('PROBABILITY');
  return valid.classify(prob, 'prob').map(function (f) {
    return f.set('pred', ee.Number(f.get('prob')).gte(PROB_THRESHOLD).int());
  });
});

var cvResult = ee.FeatureCollection(cvPredictions).flatten();
var cm = cvResult.errorMatrix(RESPONSE, 'pred');
print('Confusion matrix:', cm);
print('Accuracy:', cm.accuracy());

// Manual metrics
var arr = ee.Array(cm.array());
var TP = arr.get([1, 1]);
var FP = arr.get([0, 1]);
var FN = arr.get([1, 0]);

var precision = TP.divide(TP.add(FP));
var recall = TP.divide(TP.add(FN));
var f1 = precision.multiply(recall).multiply(2).divide(precision.add(recall));

print('Precision:', precision);
print('Recall:', recall);
print('F1 Score:', f1);

// ==========================================================
// FINAL MODEL
// ==========================================================
var finalRF = ee.Classifier.smileRandomForest(RF_PARAMS).train({
  features: training_all,
  classProperty: RESPONSE,
  inputProperties: FEATURE_BANDS
});
print('Final model trained');

// ==========================================================
// VARIABLE IMPORTANCE (Feature Ranking)
// ==========================================================
var importance = finalRF.explain().get('importance');
var importanceDict = ee.Dictionary(importance);
var keys = importanceDict.keys();

// Convert dictionary to feature collection
var importanceFeatures = keys.map(function (k) {
  var key = ee.String(k);
  var val = ee.Number(importanceDict.get(key));
  return ee.Feature(null, {
    feature: key,
    importance: val
  });
});

var sortedFC = ee.FeatureCollection(importanceFeatures).sort('importance', false);
print('Feature Importance (sorted):', sortedFC.limit(10));

// Bar chart for top features
var chart = ui.Chart.feature.byFeature({
  features: sortedFC.limit(66),
  xProperty: 'feature',
  yProperties: ['importance']
})
.setChartType('ColumnChart')
.setOptions({
  title: 'Top 10 Feature Importance',
  hAxis: {title: 'Feature'},
  vAxis: {title: 'Importance'},
  legend: {position: 'none'},
  colors: ['#444']
});
print(chart);

// (Optional) Export variable importance as CSV
// Export.table.toDrive({
//   collection: sortedFC,
//   description: 'export_RF_feature_importance',
//   fileFormat: 'CSV'
// });

// ==========================================================
// STACK BUILDER (No .reproject())
// ==========================================================
function buildStackForYear(y) {
  y = ee.Number(y);
  var emb = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
    .filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y.add(1), 1, 1))
    .mosaic();

  var depth = ee.Image(ASSET_ROOT + '/InputArea/raster_INA_ACA_Bathymetry')
    .multiply(-1).rename('depth');

  var slope = ee.Terrain.slope(depth).rename('slope');

  var rugosity = slope.reduceNeighborhood({
    reducer: ee.Reducer.stdDev(),
    kernel: ee.Kernel.circle(9, 'meters')
  }).unmask(0).rename('rugosity');

  // Distance to land. EXTERNAL ASSET: hosted by the Allen Coral Atlas
  // team, not the author's own -- requires your own access to this
  // specific shared asset, not just a path change.
  var distToLand = ee.Image('projects/coral_atlas/global_datasets/osm_distToLand_indo')
    .rename('distToLand');

  var wave = ee.Image(ASSET_ROOT + '/InputArea/era5_gebco_multiband_2024_450m_indo')
    .select([0, 1, 2], ['wElevation', 'wHeight', 'wPeriod']);

  var stack = emb;
  if (STATIC_LAYER_CONFIG['depth'].include) stack = stack.addBands(depth);
  if (STATIC_LAYER_CONFIG['slope'].include) stack = stack.addBands(slope);
  if (STATIC_LAYER_CONFIG['rugosity'].include) stack = stack.addBands(rugosity);
  if (STATIC_LAYER_CONFIG['distToLand'].include) stack = stack.addBands(distToLand);
  if (STATIC_LAYER_CONFIG['wave'].include) stack = stack.addBands(wave);

  var mask = depth.lt(0).and(depth.gte(-500));
  return stack.updateMask(mask).clip(AOI_GEOM);
}

// ==========================================================
// APPLY MODEL BY YEAR
// ==========================================================
APPLY_YEARS.evaluate(function (years) {
  years.forEach(function (y) {
    var yearNum = parseInt(y);
    var stack = buildStackForYear(yearNum).select(FEATURE_BANDS);
    var prob = stack.classify(finalRF.setOutputMode('PROBABILITY')).rename('prob');
    var thr = APPLY_THRESHOLD[yearNum] || PROB_THRESHOLD;
    var mask = prob.gte(thr).rename('seagrass').selfMask();

    // Add to map (optional)
    Map.addLayer(prob, {min: 0, max: 1, palette: ['#d73027', '#fee08b', '#1a9850']}, 'Prob ' + yearNum, false);
    Map.addLayer(mask, {palette: ['#1a9850']}, 'Mask ' + yearNum, false);

    // Export: probability
    Export.image.toAsset({
      image: prob.clip(AOI_GEOM),
      description: 'export_RF_prob_' + yearNum,
      assetId: EXPORT_BASE + '/RF_probability04022026_' + yearNum,
      region: AOI_GEOM,
      scale: SCALE,
      crs: 'EPSG:4326',
      maxPixels: 1e13
    });

    // Export: binary mask
    Export.image.toAsset({
      image: mask.clip(AOI_GEOM),
      description: 'export_RF_mask_' + yearNum,
      assetId: EXPORT_BASE + '/RF_seagrass_mask_' + yearNum,
      region: AOI_GEOM,
      scale: SCALE,
      crs: 'EPSG:4326',
      maxPixels: 1e13
    });

    print('Export queued for year:', yearNum);
  });

  Map.centerObject(AOI_GEOM, 6);
});