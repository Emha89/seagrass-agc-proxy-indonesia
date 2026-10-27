/*****************************************************************************
 * 03_training_points_morphology.js
 * Samples the GSE bands, depth and distance to land at the points with a leaf
 * morphology label (training_MORPH3_points_ALL.csv from
 * R/workflow/07_morphology_model.R, uploaded to GEE), for each survey year.
 *
 * Inputs  (ASSET_ROOT)  training_points_morphology, study_regions, bathymetry
 * Outputs (ASSET_ROOT)  training_gse_morphology_<year>
 *
 * Set ASSET_ROOT below to your own GEE asset folder. All paths are built
 * from it, except the shared Allen Coral Atlas distance-to-land layer,
 * which needs its own access.
 ****************************************************************************/

// ================================
// CONFIGURATION
// ================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var SCALE = 10;

// Morphology training samples (already labelled)
var POINTS = ee.FeatureCollection(
  ASSET_ROOT + '/training_points_morphology'
);

// AOI for masking predictors
var AOI = ee.FeatureCollection(
  ASSET_ROOT + '/study_regions'
).geometry();

// Embedding collection
var EMB_COLL = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';

// Export prefix
var EXPORT_BASE = ASSET_ROOT + '/training_gse_morphology_';


// ================================
// DIAGNOSTICS: RAW SAMPLE COUNTS
// ================================
print("Total Morph3 points:", POINTS.size());
print("Morph3 class histogram:", POINTS.aggregate_histogram("morph3"));
print("Year distribution:", POINTS.aggregate_histogram("year"));


// ================================
// DEFINE PREDICTOR BANDS
// ================================

// GSE 64 bands (native embedding band names, A00-A63)
var GSE_BANDS = ee.List.sequence(0, 63).map(function(i){
  return ee.String('A').cat(ee.Number(i).format('%02d'));
});

// Static layers
// Bathymetry (your own copy, under ASSET_ROOT)
var depth = ee.Image(
  ASSET_ROOT + '/bathymetry'
).multiply(-1).rename("depth");

// Distance to land (shared Allen Coral Atlas asset, needs its own access)
var distToLand = ee.Image(
  'projects/coral_atlas/global_datasets/osm_distToLand_indo'
).rename("distToLand");

// Predictor band list
var FEATURE_BANDS = GSE_BANDS.cat(["depth", "distToLand"]);
print("Predictors for Morph3:", FEATURE_BANDS);


// ================================
// DEPTH MASK (<= 500m)
// ================================
var depthMask = depth.lte(500);


// ================================
// UNIQUE YEARS LOOP
// ================================
var uniqueYears = ee.List(
  POINTS.aggregate_array("year")
).distinct().sort();

print("Unique years in Morph3 dataset:", uniqueYears);


// ================================
// LOOP PER YEAR + EXPORT + DISPLAY
// ================================
uniqueYears.evaluate(function(yearList){

  yearList.forEach(function(year){

    var y = ee.Number(year);
    var yearStr = y.format().getInfo();

    // --------------------------------
    // Filter morph points for this year
    // --------------------------------
    var ptsYear = POINTS.filter(ee.Filter.eq("year", y));

    // --------------------------------
    // Load yearly embedding image
    // --------------------------------
    var emb = ee.ImageCollection(EMB_COLL)
      .filterDate(
        ee.Date.fromYMD(y, 1, 1),
        ee.Date.fromYMD(y.add(1), 1, 1)
      )
      .filterBounds(AOI)
      .mosaic();

    // --------------------------------
    // Build predictor stack
    // --------------------------------
    var stack = emb
      .select(GSE_BANDS)
      .addBands(depth)
      .addBands(distToLand)
      .updateMask(depthMask)
      .clip(AOI);

    // --------------------------------
    // Sample predictors at training points
    // --------------------------------
    var sampled = stack.sampleRegions({
      collection: ptsYear,
      properties: ["gee_id", "morph3", "year"],
      scale: SCALE,
      geometries: true
    });

    // --------------------------------
    // REPORT SUMMARY (BEFORE EXPORT)
    // --------------------------------
    print("======================================");
    print("Year:", yearStr);
    print("Initial morph points:", ptsYear.size());
    print("Sampled with predictors:", sampled.size());
    print("Morph3 histogram:", sampled.aggregate_histogram("morph3"));

    // --------------------------------
    // DISPLAY TRAINING POINTS ON MAP (PER YEAR)
    // --------------------------------
    Map.addLayer(
      sampled,
      {color: "cyan"},
      "Morph3 Training Samples " + yearStr,
      false
    );

    // --------------------------------
    // EXPORT TABLE PER YEAR
    // --------------------------------
    Export.table.toAsset({
      collection: sampled,
      description: "export_training_morph3_" + yearStr,
      assetId: EXPORT_BASE + yearStr
    });

    print("Export queued:", EXPORT_BASE + yearStr);

  });

  // Map center once
  Map.centerObject(AOI, 6);
});