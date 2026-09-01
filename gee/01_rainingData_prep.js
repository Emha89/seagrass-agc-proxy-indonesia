/****************************************************
 * 01_trainingData_prep.js
 * Samples the Google Satellite Embedding (GSE) annual composite plus
 * static terrain predictors (depth, distance to land, region ID) at a
 * REFINED training points asset, for each survey year, and exports one
 * table asset per year.
 *
 * PIPELINE POSITION: this runs AFTER the R-side pipeline, not before.
 * The points asset (PAsamples, dated) are field points that have gone
 * through the full R analysis (tuning, filtering) and been re-uploaded
 * to GEE -- see 00_extractGSE_trainingPoints.js for the actual earlier
 * step that first extracts GSE + terrain values at the ORIGINAL training
 * points and exports the GSE_training_<year>_CSV.csv files the R
 * pipeline (R/01_occurrence_PA/dataPrep_func_PA2.R and others) is built
 * on. This script's own export (a GEE table asset, not a CSV) appears to
 * feed a later deployment step not yet reviewed.
 *
 * Companion script: 01_trainingData_prepMorpho -- a variant of this
 * script that also joins morphology labels. Not reviewed yet.
 *
 * ASSET PORTABILITY: this script originally referenced the author's own
 * private GEE assets. Below, ASSET_ROOT collects the four paths that
 * shared one project/folder root into a single placeholder -- replace
 * ASSET_ROOT with your own GEE project and asset folder. The
 * Allen Coral Atlas distToLand layer is a separate, externally-hosted
 * shared asset (not the author's own) -- using it requires your own
 * access to that specific Coral Atlas asset, not just a path swap.
 *
 * The two unused "Imports" panel assets from the original script
 * (indo_area, ecoreg, osm -- none referenced anywhere in the script body)
 * have been left out here.
 ****************************************************/

// =====================================================================
// ASSET PATHS -- replace ASSET_ROOT with your own GEE project/folder
// =====================================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var SCALE = 10;

var POINTS = ee.FeatureCollection(
  ASSET_ROOT + '/TrainingPoint/PAsamples'
);

var AOI = ee.FeatureCollection(
  ASSET_ROOT + '/InputArea/R2_case_study'
);

var EMB_COLL = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL';
var EXPORT_BASE = ASSET_ROOT + '/output/training_embedding_depth3_';


// ================================
// REPORT RAW SAMPLE SIZE
// ================================
print("================================================");
print("DATA PREP STARTED (BASE TRAINING SET)");
print("Total POINTS:", POINTS.size());
print("PA Histogram (raw):", POINTS.aggregate_histogram("PA"));
print("================================================");

// Fix PA numeric type
POINTS = POINTS.map(function(f){
  return f.set("PA", ee.Number.parse(f.get("PA")).int());
});

print("PA Histogram (after numeric fix):",
  POINTS.aggregate_histogram("PA")
);


// ================================
// STATIC LAYERS
// ================================
// EXTERNAL ASSET: your own private copy, path under ASSET_ROOT
var raw_depth = ee.Image(
  ASSET_ROOT + '/InputArea/raster_INA_ACA_Bathymetry'
);

var depth = raw_depth.multiply(-1).rename('depth');
var depthMask = depth.lte(500);

// Distance to land. EXTERNAL ASSET: hosted by the Allen Coral Atlas team,
// not the author's own -- requires your own access to this specific
// shared asset (contact the Coral Atlas project), not just a path change.
var distToLand = ee.Image(
  'projects/coral_atlas/global_datasets/osm_distToLand_indo'
).rename('distToLand')
 .updateMask(depthMask);

// AOI loc_id raster
var loc_dict = ee.Dictionary({
  'Ayau': 1, 'Rote': 2, 'Karimunjawa': 3,
  'Bintan': 4, 'Komodo': 5, 'Manado': 6
});

var loc_fc_coded = AOI.map(function (f) {
  return f.set('loc_id', loc_dict.get(f.getString('loc')));
});

var loc_raster = loc_fc_coded.reduceToImage({
  properties: ['loc_id'],
  reducer: ee.Reducer.first()
}).updateMask(depthMask)
  .clip(AOI.geometry())
  .rename('loc_id');

// Terrain stack
var terrainStack = ee.Image([])
  .addBands(depth)
  .addBands(distToLand)
  .addBands(loc_raster);

print("Static bands exported:", terrainStack.bandNames());


// ================================
// MAP DISPLAY STATIC LAYERS
// ================================
Map.centerObject(AOI, 6);

Map.addLayer(depth, {
  min: 0, max: 500,
  palette: ['blue', 'cyan', 'yellow']
}, "Depth (masked)", false);

Map.addLayer(distToLand, {
  min: 0, max: 50000,
  palette: ['green', 'yellow', 'red']
}, "Distance to Land", false);

Map.addLayer(loc_raster.randomVisualizer(), {}, "loc_id raster", false);


// ================================
// YEAR LOOP + EXPORT
// ================================
var uniqueYears = ee.List(
  POINTS.aggregate_array('year_gt')
).distinct().sort();

print("Unique years:", uniqueYears);


// ================================
// EXPORT PER YEAR
// ================================
uniqueYears.evaluate(function(yearList){

  print("================================================");
  print("EXPORT STARTING...");
  print("Total export years:", yearList.length);
  print("================================================");

  yearList.forEach(function(year){

    var y = ee.Number(year);
    var yearStr = y.format().getInfo();

    // Filter points by year
    var pointsYear = POINTS.filter(
      ee.Filter.eq('year_gt', y)
    );

    // Report BEFORE sampling
    print("----------------------------------");
    print("YEAR:", yearStr);
    print("Samples before extraction:", pointsYear.size());
    print("PA histogram:", pointsYear.aggregate_histogram("PA"));

    // Load yearly GSE embedding
    var emb = ee.ImageCollection(EMB_COLL)
      .filterDate(
        ee.Date.fromYMD(y, 1, 1),
        ee.Date.fromYMD(y.add(1), 1, 1)
      )
      .filterBounds(AOI.geometry())
      .mosaic();

    // Combine embedding + static predictors
    var fullImage = emb.addBands(terrainStack);

    // Sample predictors at points
    var sampled = fullImage.sampleRegions({
      collection: pointsYear,
      properties: [
        'PA', 'tSPC', 'gee_id',
        'AGB_pred', 'AGB_low', 'AGB_up',
        'AGC_pred', 'AGC_low', 'AGC_up',
        'AGB_CIwidt', 'AGC_CIwidt',
        'CI_ratio', 'carbon_ind'
      ],
      scale: SCALE,
      geometries: true
    });

    // Report AFTER sampling
    print("Sampled rows exported:", sampled.size());
    print("PA histogram sampled:", sampled.aggregate_histogram("PA"));

    // ================================
    // DISPLAY POINTS ON MAP (Preview)
    // ================================
    Map.addLayer(sampled, {
      color: "yellow"
    }, "Sample Points " + yearStr, false);

    // Export table
    Export.table.toAsset({
      collection: sampled,
      description: 'export_training_RF_' + yearStr,
      assetId: EXPORT_BASE + yearStr
    });

    print("Export queued:", EXPORT_BASE + yearStr);

  });

});