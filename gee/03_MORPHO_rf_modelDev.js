/************************************************************
 * 03_MORPHO_rf_modelDev.js
 * RF MORPH3 (multi-class) model: training, K-fold cross-validation,
 * apply + export per year.
 *
 * Output per year:
 *   - Probability bands per class: P_<class>
 *   - Predicted class (index)
 *
 * Includes an optional interactive viewer: probability gradient map for
 * a selectable class, with an adjustable confidence-threshold mask, plus
 * a UI dropdown + slider.
 *
 * PIPELINE POSITION: reads training_morph3_<year> assets produced by
 * 01_trainingData_prepMorpho.js (TRAIN_ASSET_PREFIX below matches that
 * script's EXPORT_BASE exactly). RF_PARAMS below (ntree=300, mtry=6,
 * nodesize=7, sample_fraction=0.7) are presumed to be the R-side tuned
 * hyperparameters for the morphology model, hardcoded here for GEE
 * deployment. STATIC_BANDS (depth + distToLand only, no slope/rugosity/
 * wave) matches the R-side morphology predictor set exactly
 * (R/02_morphology/rf_model_combined_MORPH3.R).
 *
 * NAMING MISMATCH (display only, does not affect the actual MORPH3
 * exports): PERSISTENCE_MASK below points at
 * Seagrass_MaxProb_0_8_2017_2024_mask, but 02_extent_2017-2025.js (kept
 * exactly as originally written) exports that mask as
 * Seagrass_MaxProb_0_8_2017_2025_mask -- "2024" vs "2025". This mask is
 * used only for the optional interactive map preview (explicitly noted
 * in the original script as "DISPLAY ONLY" / "EXPORTS UNCHANGED, NO
 * persistence mask"), so a mismatch here would affect the preview layer
 * only, not the exported probability/prediction assets. Kept exactly as
 * in the original script since I can't tell which year suffix is meant
 * to be correct without access to the actual assets.
 *
 * MINOR NOTE: the try/catch around each year's training asset load
 * likely doesn't catch "asset not found" errors the way it looks like it
 * should -- ee.FeatureCollection() on a bad path doesn't throw
 * synchronously in GEE's client-side JS, the error only surfaces later
 * when that reference is actually evaluated. Left as-is; not a
 * functional problem, just worth knowing if a year silently fails to
 * load.
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

var RESPONSE = 'morph3';
var K_FOLDS  = 5;
var SCALE    = 10;
var SEED     = 42;

var START_YEAR = 2017;
var END_YEAR   = 2025;

// Persistence mask (DISPLAY ONLY) -- see NAMING MISMATCH note above
var PERSISTENCE_MASK =
  ASSET_ROOT + '/output/Seagrass_MaxProb_0_8_2017_2024_mask';

// Predictors (match R morph model: depth + distToLand + 64 GSE bands)
var USE_SLOPE_RUGOSITY = false;  // keep FALSE to mirror R (depth+distToLand only)
var USE_WAVE           = false;  // keep FALSE (excluded in R)

// Export
var EXPORT_BASE = ASSET_ROOT + '/output';

// Use the SAME training asset prefix produced by 01_trainingData_prepMorpho.js
var TRAIN_ASSET_PREFIX = ASSET_ROOT + '/output/training_morph3_';

// Years you actually have training points
var available_years = [2018, 2019, 2020, 2021, 2022, 2023];

// RF params (see PIPELINE POSITION note above)
var RF_PARAMS = {
  numberOfTrees: 300,
  variablesPerSplit: 6,     // mtry analogue
  minLeafPopulation: 7,     // nodesize analogue (approx)
  bagFraction: 0.7,
  seed: SEED
};

// ==========================================================
// VISUALIZATION CONTROLS
// ==========================================================
// Default class to visualize on map (must match band name after rename: 'P_<class>')
var DEFAULT_CLASS_TO_VIEW = 'P_mono_Ea';  // <-- change freely
var DEFAULT_THR = 0.7;                   // threshold for confidence mask
var SHOW_UI = true;                      // set false if you don't want UI controls

// Palette for probability
var PROB_VIS = {min: 0, max: 1, palette: ['red', 'yellow', 'green']};

// ==========================================================
// LOAD AOI
// ==========================================================
var AOI_FC = ee.FeatureCollection(
  ASSET_ROOT + '/InputArea/R2_case_study'
);
var AOI_GEOM = AOI_FC.geometry();
Map.addLayer(AOI_GEOM, {color: 'yellow'}, 'AOI R2 Case Study');

// Load persistence mask as image for display
var persistence = ee.Image(PERSISTENCE_MASK);

// ==========================================================
// DEFINE BANDS
// ==========================================================
var GSE_BANDS = ee.List.sequence(0, 63).map(function (i) {
  return ee.String('A').cat(ee.Number(i).format('%02d'));
});

var STATIC_BANDS = ['depth', 'distToLand'];
if (USE_SLOPE_RUGOSITY) STATIC_BANDS = STATIC_BANDS.concat(['slope', 'rugosity']);
if (USE_WAVE) STATIC_BANDS = STATIC_BANDS.concat(['wElevation', 'wHeight', 'wPeriod']);

var FEATURE_BANDS = GSE_BANDS.cat(STATIC_BANDS);
print('MORPH3 uses bands:', FEATURE_BANDS);

// ==========================================================
// LOAD TRAINING DATASET FROM MULTIPLE YEARS
// ==========================================================
var training_all = ee.FeatureCollection([]);

available_years.forEach(function (y) {
  var assetPath = TRAIN_ASSET_PREFIX + y;
  try {
    var fc = ee.FeatureCollection(assetPath);
    training_all = training_all.merge(fc);
    print('Loaded training data for year:', y, '|', assetPath);
  } catch (e) {
    print('Asset not found:', assetPath);
  }
});

print('Total training samples (raw, incl PA=0):', training_all.size());


// ==========================================================
// FILTER: KEEP ONLY MORPH3-LABELLED POINTS (SEAGRASS ONLY)
// ==========================================================

// Many PA=0 points have morph3 = null -- must exclude for morph training
var totalRaw = training_all.size();

// Keep only rows where morph3 exists
training_all = training_all.filter(ee.Filter.notNull([RESPONSE]));

var totalMorph = training_all.size();

print('Training samples with morph3 label:', totalMorph);
print('Removed non-morph points (PA=0):', totalRaw.subtract(totalMorph));


// ==========================================================
// REMOVE FEATURES WITH NULL PREDICTORS
// ==========================================================
function notNullFor(list) {
  return ee.Filter.and.apply(null, ee.List(list).map(function (b) {
    return ee.Filter.notNull([b]);
  }));
}

// Filter rows with missing predictors
var totalBefore = training_all.size();

training_all = training_all.filter(
  notNullFor(FEATURE_BANDS.add(RESPONSE))
);

var totalAfter = training_all.size();

print('Rows removed due to NA predictors:', totalBefore.subtract(totalAfter));
print('Final MORPH3 training size:', totalAfter);
print("morph3 class distribution:", training_all.aggregate_histogram(RESPONSE));

// ==========================================================
// CLASS LIST + SAFE ENCODING (string -> integer)
// ==========================================================
var classList = ee.List(training_all.aggregate_array(RESPONSE))
  .distinct()
  .sort();

print('morph3 classes (sorted):', classList);

// Build dictionary label -> code (0..n-1)
var idxList = ee.List.sequence(0, classList.size().subtract(1));
var labelToCode = ee.Dictionary.fromLists(classList, idxList);
var codeToLabel = ee.Dictionary.fromLists(idxList, classList); // optional

var RESPONSE_ID = 'morph3_id';

// Add numeric class property for training stability
training_all = training_all.map(function(f){
  var lbl = f.get(RESPONSE);
  var code = labelToCode.get(lbl);
  return f.set(RESPONSE_ID, code);
});

print('Encoded morph3 -> morph3_id (0..n-1)');

// ==========================================================
// K-FOLD CROSS-VALIDATION (TRAINING + EVALUATION)
// ==========================================================
var withFold = training_all.randomColumn('rand', SEED).map(function (f) {
  return f.set('fold', ee.Number(f.get('rand')).multiply(K_FOLDS).int());
});
print('K-Fold distribution:', withFold.aggregate_histogram('fold'));

var foldIndices = ee.List.sequence(0, K_FOLDS - 1);

var cvPredictions = foldIndices.map(function (k) {
  k = ee.Number(k);

  var train = withFold.filter(ee.Filter.neq('fold', k));
  var valid = withFold.filter(ee.Filter.eq('fold', k));

  var rf = ee.Classifier.smileRandomForest(RF_PARAMS).train({
    features: train,
    classProperty: RESPONSE_ID,     // numeric class
    inputProperties: FEATURE_BANDS
  });

  // Predict numeric class code
  var pred = valid.classify(rf, 'pred_id');

  return pred;
});

var cvResult = ee.FeatureCollection(cvPredictions).flatten();

// Confusion matrix (numeric classes)
var cm = cvResult.errorMatrix(RESPONSE_ID, 'pred_id');
print('Confusion matrix (K-fold):', cm);
print('Overall accuracy:', cm.accuracy());
print('Kappa:', cm.kappa());
print('Producer accuracy (recall) per class:', cm.producersAccuracy());
print('Consumer accuracy (precision) per class:', cm.consumersAccuracy());

// ==========================================================
// FINAL MODEL TRAINING (ON FULL DATASET)
// ==========================================================
var finalRF = ee.Classifier.smileRandomForest(RF_PARAMS).train({
  features: training_all,
  classProperty: RESPONSE_ID,       // numeric class
  inputProperties: FEATURE_BANDS
});
print('Final MORPH3 model trained.');

// ==========================================================
// VARIABLE IMPORTANCE
// ==========================================================
var importance = ee.Dictionary(finalRF.explain().get('importance'));
print('Feature Importance (raw):', importance);

var keys = importance.keys();
var sortedFC = ee.FeatureCollection(keys.map(function (k) {
  k = ee.String(k);
  return ee.Feature(null, {feature: k, importance: ee.Number(importance.get(k))});
})).sort('importance', false);

print('Sorted Feature Importance (descending):', sortedFC);

// Top 10 plot
var top10 = sortedFC.limit(10);
var chart = ui.Chart.feature.byFeature({
  features: top10,
  xProperty: 'feature',
  yProperties: ['importance']
})
.setChartType('ColumnChart')
.setOptions({
  title: 'Top 10 Variable Importance (RF MORPH3)',
  hAxis: {title: 'Feature'},
  vAxis: {title: 'Importance'},
  legend: {position: 'none'},
  colors: ['grey']
});
print(chart);

// ==========================================================
// BUILD STACK FUNCTION FOR PREDICTION PER YEAR
// ==========================================================
function buildStackForYear(y) {
  var emb = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
    .filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(ee.Number(y).add(1), 1, 1))
    .mosaic();

  var depth = ee.Image(ASSET_ROOT + '/InputArea/raster_INA_ACA_Bathymetry')
    .multiply(-1).rename('depth');

  // Distance to land. EXTERNAL ASSET: hosted by the Allen Coral Atlas
  // team, not the author's own -- requires your own access to this
  // specific shared asset, not just a path change.
  var distToLand = ee.Image('projects/coral_atlas/global_datasets/osm_distToLand_indo')
    .rename('distToLand');

  var stack = emb
    .addBands(depth)
    .addBands(distToLand);

  if (USE_SLOPE_RUGOSITY) {
    var slope = ee.Terrain.slope(depth).rename('slope');
    var rugosity = slope.reduceNeighborhood({
      reducer: ee.Reducer.stdDev(),
      kernel: ee.Kernel.circle(9, 'meters')
    }).unmask(0).rename('rugosity');
    stack = stack.addBands(slope).addBands(rugosity);
  }

  if (USE_WAVE) {
    var wave = ee.Image(ASSET_ROOT + '/InputArea/era5_gebco_multiband_2024_450m_indo')
      .select([0, 1, 2], ['wElevation', 'wHeight', 'wPeriod']);
    stack = stack.addBands(wave);
  }

  return stack
    .clip(AOI_GEOM)
    .updateMask(depth.abs().lte(500));
}

// ==========================================================
// UI HELPERS (optional)
// ==========================================================
var uiState = {
  selectedBand: DEFAULT_CLASS_TO_VIEW,
  thr: DEFAULT_THR,
  currentYear: null,
  probBands: null
};

function updateMapLayers() {
  if (!uiState.probBands || uiState.currentYear === null) return;

  var band = uiState.selectedBand;
  var thr  = uiState.thr;

  while (Map.layers().length() > 1) {
    Map.layers().remove(Map.layers().get(1));
  }

  // Apply persistence mask for DISPLAY ONLY
  var probDisplay = uiState.probBands.select(band).updateMask(persistence);

  Map.addLayer(
    probDisplay,
    PROB_VIS,
    'Prob ' + band + ' | ' + uiState.currentYear,
    true
  );

  // Threshold mask also uses persistence (display only)
  var mask = probDisplay.gte(thr).selfMask();

  Map.addLayer(
    mask,
    {palette: ['blue']},
    'Mask ' + band + ' >= ' + thr.toFixed(2) + ' | ' + uiState.currentYear,
    false
  );
}

// ==========================================================
// APPLY MODEL TO EACH YEAR AND EXPORT RESULTS
// ==========================================================
var APPLY_YEARS = ee.List.sequence(START_YEAR, END_YEAR);

APPLY_YEARS.evaluate(function (years) {

  var selector, slider, yearSelector;
  if (SHOW_UI) {
    var bandNames = classList.map(function(c){
      return ee.String('P_').cat(ee.String(c));
    }).getInfo();

    selector = ui.Select({
      items: bandNames,
      value: DEFAULT_CLASS_TO_VIEW,
      onChange: function(v){
        uiState.selectedBand = v;
        updateMapLayers();
      }
    });

    slider = ui.Slider({
      min: 0, max: 1, step: 0.01,
      value: DEFAULT_THR,
      onChange: function(v){
        uiState.thr = v;
        updateMapLayers();
      }
    });

    yearSelector = ui.Select({
      items: years.map(function(yy){ return String(yy); }),
      value: String(years[years.length - 1]),
      onChange: function(v){
        var y = parseInt(v, 10);
        uiState.currentYear = y;

        var stack = buildStackForYear(y).select(FEATURE_BANDS);
        var probArr = stack.classify(finalRF.setOutputMode('MULTIPROBABILITY'));
        var probBands = probArr.arrayFlatten([classList]);

        var pNames = classList.map(function(c){
          return ee.String('P_').cat(ee.String(c));
        });
        probBands = probBands.rename(pNames);

        uiState.probBands = probBands;
        updateMapLayers();
      }
    });

    print('MORPH3 Probability Viewer');
    print('Select class probability band:', selector);
    print('Set confidence threshold:', slider);
    print('Preview year:', yearSelector);
  }

  years.forEach(function (y) {
    print('Applying MORPH3 model for year:', y);

    var stack = buildStackForYear(y).select(FEATURE_BANDS);

    var probArr = stack.classify(finalRF.setOutputMode('MULTIPROBABILITY'));

    var probBands = probArr.arrayFlatten([classList]);

    var pNames = classList.map(function(c){
      return ee.String('P_').cat(ee.String(c));
    });
    probBands = probBands.rename(pNames);

    var predIdx = probArr.arrayArgmax().arrayGet([0]).toInt().rename('morph3_pred_id');

    // Initialize UI state with persistence-masked display
    if (uiState.currentYear === null) {
      uiState.currentYear = y;
      uiState.probBands = probBands;
      if (SHOW_UI) updateMapLayers();
    }

    // Exports unchanged (full, no persistence mask)
    Export.image.toAsset({
      image: probBands.clip(AOI_GEOM),
      description: 'export_MORPH3_probs_' + y,
      assetId: EXPORT_BASE + '/MORPH3_probs04022026_' + y,
      region: AOI_GEOM,
      scale: SCALE,
      crs: 'EPSG:4326',
      maxPixels: 1e13
    });

    Export.image.toAsset({
      image: predIdx.clip(AOI_GEOM),
      description: 'export_MORPH3_pred_' + y,
      assetId: EXPORT_BASE + '/MORPH3_pred_' + y,
      region: AOI_GEOM,
      scale: SCALE,
      crs: 'EPSG:4326',
      maxPixels: 1e13
    });

    print('Export queued MORPH3 for year:', y);
  });

  Map.centerObject(AOI_GEOM, 7);
});