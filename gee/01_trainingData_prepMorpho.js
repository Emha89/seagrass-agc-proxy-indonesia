/****************************************************
 * 01_trainingData_prepMorpho.js
 * Samples the Google Satellite Embedding (GSE) annual composite plus
 * depth and distance-to-land at an already-labelled morphology (morph3)
 * training points asset, for each survey year, and exports one table
 * asset per year.
 *
 * PIPELINE POSITION: companion to 01_trainingData_prep.js -- same
 * pattern (a refined, already-labelled points asset is re-sampled
 * against GSE and exported as a GEE table asset, not a CSV), but for the
 * morphology stage instead of PA. Runs AFTER the R-side morphology
 * pipeline (R/02_morphology/), using points that already carry a morph3
 * label. Feeds a later deployment step not yet reviewed.
 *
 * NOTE: the GSE bands here keep their native names (A00-A63, no "GSE_"
 * prefix), unlike the CSV files used throughout the R pipeline
 * (GSE_A00-GSE_A63). The "GSE_" prefix is added deliberately in
 * 00_extractGSE_trainingPoints.js; the raw embedding image's bands are
 * natively named A00-A63. Not a conflict with the R-side pipeline (this
 * script's output is a separate GEE asset, not one of the reviewed CSV
 * inputs) -- just worth knowing if you're comparing the two.
 *
 * ASSET PORTABILITY: this script originally referenced the author's own
 * private GEE assets. Below, ASSET_ROOT collects the paths that shared
 * one project/folder root into a single placeholder -- replace
 * ASSET_ROOT with your own GEE project and asset folder. The
 * Allen Coral Atlas distToLand layer is a separate, externally-hosted
 * shared asset (not the author's own) -- using it requires your own
 * access to that specific Coral Atlas asset, not just a path swap.
 *
 * The three original "Imports" panel assets (indo_area, ecoreg, osm) are
 * not referenced anywhere in this script's body and have been left out.
 ****************************************************/

// ================================
// CONFIGURATION
// ================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var SCALE = 10;

// Morphology training samples (already labelled)
var POINTS = ee.FeatureCollection(
  ASSET_ROOT + '/TrainingPoint/morphoSample'
);

// AOI for masking predictors
var AOI = ee.FeatureCollection(
  ASSET_ROOT + '/InputArea/R2_case_study'
).geometry();

// Embedding collection
var EMB_COLL = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';

// Export prefix
var EXPORT_BASE = ASSET_ROOT + '/output/training_morph3_';


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
// EXTERNAL ASSET: your own private copy, path under ASSET_ROOT
var depth = ee.Image(
  ASSET_ROOT + '/InputArea/raster_INA_ACA_Bathymetry'
).multiply(-1).rename("depth");

// Distance to land. EXTERNAL ASSET: hosted by the Allen Coral Atlas team,
// not the author's own -- requires your own access to this specific
// shared asset (contact the Coral Atlas project), not just a path change.
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