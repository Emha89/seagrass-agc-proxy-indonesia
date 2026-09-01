/************************************************************
 * 02_extent_2017-2025.js
 * Seagrass Persistence Map (2017-2024) -- Year Count Based
 *
 * - Loads yearly RF probability images (2017-2024, from
 *   02_PROB_rf_modelDev.js)
 * - Binarizes each image: 1 if prob >= threshold, else 0
 * - Counts number of years with seagrass presence
 * - Defines persistence based only on year count
 * - Also computes max-probability extent (ever suitable in any year) and
 *   min-probability at biomass sample points (a QA check)
 *
 * PIPELINE POSITION: reads RF_probability04022026_<year> assets produced
 * by 02_PROB_rf_modelDev.js (PROB_ASSET_PREFIX below must match that
 * script's export asset name exactly, including the date stamp).
 *
 * NAMING NOTE: exported names/descriptions in several places use the
 * suffix "2017_2025", but START_YEAR/END_YEAR below are 2017-2024 (8
 * years) -- the actual computation only covers 2017-2024; the "2025" in
 * these labels looks like a leftover from a different version rather
 * than an extra year of data. Not changed here, just flagged.
 *
 * NAMING NOTE 2: biomassAssetList below references
 * training_embedding_depth_<year> (no "3"), which differs from
 * training_embedding_depth3_<year> used as the training source in
 * 01_trainingData_prep.js / 02_PROB_rf_modelDev.js. Confirm whether this
 * is a genuinely different (older?) asset or a naming slip -- kept
 * exactly as in the original script since I can't verify which is
 * intended without access to the assets themselves.
 *
 * ASSET PORTABILITY: this script originally referenced the author's own
 * private GEE assets, all under one project/folder root -- replace
 * ASSET_ROOT below with your own GEE project and asset folder.
 ************************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var START_YEAR = 2017;
var END_YEAR = 2024;

var EXPORT_BASE = ASSET_ROOT + '/output';
var PROB_ASSET_PREFIX = EXPORT_BASE + '/RF_probability04022026_';

var PROB_PRESENCE_THRESHOLD = 0.8; // Threshold per year to consider seagrass present
var MIN_YEARS_PRESENT = 8;         // Minimum years to consider persistent
var EXPORT_DESCRIPTION = 'Seagrass_Persistence_by_YearCount_2017_2025';

var AOI = ee.FeatureCollection(ASSET_ROOT + '/InputArea/R2_case_study').geometry();
Map.addLayer(AOI, {color: 'yellow'}, 'AOI', false);

// Generate list of years
var yearList = [];
for (var y = START_YEAR; y <= END_YEAR; y++) {
  yearList.push(y);
}

// ==========================================================
// OPTIONAL: DISPLAY ALL PROBABILITY MAPS (PER YEAR)
// ==========================================================
yearList.forEach(function (y) {
  var assetId = PROB_ASSET_PREFIX + y;
  var img = ee.Image(assetId);

  Map.addLayer(img, {
    min: 0,
    max: 1,
    palette: ['#d73027', '#fee08b', '#1a9850']
  }, 'Prob Seagrass ' + y, false);  // false = layer off by default
});

// ==========================================================
// LOAD PROBABILITY IMAGES AND BINARY MASKS
// ==========================================================
var presenceImages = yearList.map(function (y) {
  var assetId = PROB_ASSET_PREFIX + y;
  var img = ee.Image(assetId); // Band is 'prob'
  var binarized = img.gte(PROB_PRESENCE_THRESHOLD); // 1 if present, 0 otherwise
  return binarized;
});

var presenceStack = ee.ImageCollection(presenceImages);


// ==========================================================
// COMPUTE PRESENCE COUNT (NUMBER OF YEARS PRESENT)
// ==========================================================
var presenceCount = presenceStack.sum().rename('presence_count');

Map.addLayer(presenceCount, {
  min: 0,
  max: yearList.length,
  palette: ['white', 'yellow', 'green']
}, 'Seagrass Years (prob >= ' + PROB_PRESENCE_THRESHOLD + ')');


// ==========================================================
// DEFINE FINAL PERSISTENCE MASK
// ==========================================================
var persistenceMask = presenceCount.gte(MIN_YEARS_PRESENT)
  .rename('seagrass_persistence')
  .selfMask();

Map.addLayer(persistenceMask, {
  palette: ['#1a9850']
}, 'Persistent Seagrass Mask (by year count only)');


// ==========================================================
// EXPORT: PERSISTENT MASK
// ==========================================================
Export.image.toAsset({
  image: persistenceMask.clip(AOI),
  description: EXPORT_DESCRIPTION + '_mask',
  assetId: EXPORT_BASE + '/seagrass_persistence_mask_byYearCount_2017_2025',
  region: AOI,
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});


// ==========================================================
// EXPORT: PRESENCE COUNT (optional, for visualization)
// ==========================================================
Export.image.toAsset({
  image: presenceCount.clip(AOI),
  description: EXPORT_DESCRIPTION + '_presence_count',
  assetId: EXPORT_BASE + '/seagrass_presence_count_2017_2025',
  region: AOI,
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});


// ==========================================================
// LOG CONFIGURATION
// ==========================================================
print('Year range:', START_YEAR, '-', END_YEAR);
print('Threshold per year for presence:', PROB_PRESENCE_THRESHOLD);
print('Minimum years required for persistence:', MIN_YEARS_PRESENT);

// ==========================================================
// COMPUTE AREA OF PERSISTENT SEAGRASS PER LOCATION (in km2)
// ==========================================================

// Load feature collection with location names
var locations = ee.FeatureCollection(ASSET_ROOT + '/InputArea/R2_case_study');

// Function to compute persistent seagrass area for a given feature
var calculatePersistentArea = function (feature) {
  var locationName = feature.get('loc');

  // Clip the persistence mask to the current feature
  var maskClipped = persistenceMask.clip(feature.geometry());

  // Calculate area in square meters
  var pixelArea = ee.Image.pixelArea().updateMask(maskClipped); // only where seagrass is persistent

  var totalArea = pixelArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: feature.geometry(),
    scale: 10,
    maxPixels: 1e13
  });

  // Convert area to square kilometers (km2)
  var areaKm2 = ee.Number(totalArea.get('area')).divide(1e6);

  return feature.set('persistent_area_km2', areaKm2);
};

// Apply function to all features
var locationsWithArea = locations.map(calculatePersistentArea);

// Print entire feature collection
print('Persistent Seagrass Area (km2) per Location [FeatureCollection]:', locationsWithArea);

// Evaluate and print each location's area to console
locationsWithArea.evaluate(function (fc) {
  print('Detailed Area per Location (km2):');
  fc.features.forEach(function (f) {
    var name = f.properties.loc;
    var area = f.properties.persistent_area_km2;
    print(name + ': ' + area.toFixed(2) + ' km2');
  });
});

// Export as CSV to Google Drive
Export.table.toDrive({
  collection: locationsWithArea,
  description: 'Persistent_Seagrass_Area_Per_Location_2017_2025_km2',
  fileFormat: 'CSV'
});


//----------------------------------------------------------

// ==========================================================
// MAXIMUM PROBABILITY ACROSS YEARS (Per-Pixel)
// Purpose: Map areas that were ever predicted suitable
//          for any seagrass species between 2017-2024
// ==========================================================

// === CONFIGURATION ===
var MAX_PROB_THRESHOLD = 0.8;  // Probability threshold for defining seagrass presence
var THRESHOLD_STR = String(MAX_PROB_THRESHOLD).replace('.', '_');  // e.g., 0.9 -> "0_9"
var EXPORT_PREFIX = 'Seagrass_MaxProb_' + THRESHOLD_STR + '_2017_2025';  // Used in all export names

// === LOAD PROBABILITY IMAGES FOR EACH YEAR ===
var probImages = yearList.map(function (y) {
  var assetId = PROB_ASSET_PREFIX + y;
  return ee.Image(assetId).select(['prob']);  // Ensure consistent band name
});

var probStack = ee.ImageCollection(probImages);

// === REDUCE TO MAXIMUM PROBABILITY VALUE PER PIXEL ACROSS YEARS ===
var maxProbImage = probStack.reduce(ee.Reducer.max()).rename('max_probability');

// === DISPLAY MAXIMUM PROBABILITY IMAGE ON MAP ===
Map.addLayer(maxProbImage, {
  min: 0,
  max: 1,
  palette: ['#d73027', '#fee08b', '#1a9850']
}, 'Maximum Probability (2017-2025)');

// ==========================================================
// DEFINE SEAGRASS MASK (Pixels with max prob >= threshold)
// ==========================================================

var seagrassMask = maxProbImage
  .gte(MAX_PROB_THRESHOLD)
  .rename('seagrass_presence')
  .selfMask();  // Mask non-seagrass areas (false values)

// Add mask to the map
Map.addLayer(seagrassMask, {
  palette: ['#ffff00'],
  opacity: 0.5
}, 'Seagrass Habitat (MaxProb >= ' + MAX_PROB_THRESHOLD + ')');

// ==========================================================
// COMPUTE TOTAL SEAGRASS AREA (MAX PROB >= threshold) PER LOCATION
// ==========================================================

var calculateMaxProbArea = function (feature) {
  var maskClipped = seagrassMask.clip(feature.geometry());

  var pixelArea = ee.Image.pixelArea().updateMask(maskClipped);

  var totalArea = pixelArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: feature.geometry(),
    scale: 10,
    maxPixels: 1e13
  });

  var areaKm2 = ee.Number(totalArea.get('area')).divide(1e6);

  return feature.set('maxprob_area_km2', areaKm2);
};

// Apply the function to all locations
var locationsMaxProbArea = locations.map(calculateMaxProbArea);

// Print result to console
print('Seagrass Area (MaxProb >= ' + MAX_PROB_THRESHOLD + ') per Location:', locationsMaxProbArea);

// Log summary to console (client-side)
locationsMaxProbArea.evaluate(function (fc) {
  print('Area per Location maxThres (km2):');
  fc.features.forEach(function (f) {
    var name = f.properties.loc;
    var area = f.properties.maxprob_area_km2;
    print(name + ': ' + area.toFixed(2) + ' km2');
  });
});

// ==========================================================
// EXPORT: MAXIMUM PROBABILITY IMAGE TO ASSET
// ==========================================================
Export.image.toAsset({
  image: maxProbImage.clip(AOI),
  description: EXPORT_PREFIX + '_maxprob_image',
  assetId: EXPORT_BASE + '/' + EXPORT_PREFIX + '_maxprob_image',
  region: AOI,
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// ==========================================================
// EXPORT: SEAGRASS MASK (MAX PROB >= threshold) TO ASSET
// ==========================================================
Export.image.toAsset({
  image: seagrassMask.clip(AOI),
  description: EXPORT_PREFIX + '_mask',
  assetId: EXPORT_BASE + '/' + EXPORT_PREFIX + '_mask',
  region: AOI,
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// ==========================================================
// EXPORT: TOTAL AREA PER LOCATION TO GOOGLE DRIVE (CSV)
// ==========================================================
Export.table.toDrive({
  collection: locationsMaxProbArea,
  description: EXPORT_PREFIX + '_area_per_location_km2',
  fileFormat: 'CSV'
});


// ==========================================================
// MINIMUM PROBABILITY ACROSS YEARS (Per-Pixel)
// ==========================================================

// Load RF probability images for each year using consistent band name
var probImages2 = yearList.map(function (y) {
  var assetId = PROB_ASSET_PREFIX + y;
  return ee.Image(assetId).select(['prob']);  // Select only the 'prob' band
});

var probStack2 = ee.ImageCollection(probImages2);

// Reduce to minimum value across years per pixel
var minProbImage = probStack2.reduce(ee.Reducer.min()).rename('min_probability');

// Visualize
Map.addLayer(minProbImage, {
  min: 0,
  max: 1,
  palette: ['#d7191c', '#fdae61', '#1a9641']
}, 'Minimum Probability (2017-2025)');


// ==========================================================
// SAMPLE MINIMUM PROBABILITY AT BIOMASS LOCATIONS
// ==========================================================

// List of field sample asset paths (2018-2023). NOTE: uses
// training_embedding_depth_<year> (no "3") -- see NAMING NOTE 2 at the
// top of this file.
var biomassAssetList = [
  ASSET_ROOT + '/output/training_embedding_depth_2018',
  ASSET_ROOT + '/output/training_embedding_depth_2019',
  ASSET_ROOT + '/output/training_embedding_depth_2020',
  ASSET_ROOT + '/output/training_embedding_depth_2021',
  ASSET_ROOT + '/output/training_embedding_depth_2022',
  ASSET_ROOT + '/output/training_embedding_depth_2023'
];

// Load and merge all feature collections
var training_all = ee.FeatureCollection([]);
biomassAssetList.forEach(function (assetId) {
  var fc = ee.FeatureCollection(assetId);
  training_all = training_all.merge(fc);
});

// Keep only field points that have biomass information (AGB_pred column not null)
var biomassPoints = training_all.filter(ee.Filter.notNull(['AGB_pred']));
print('Biomass sample points with AGB_pred:', biomassPoints.size());

// Sample the min probability image at those biomass point locations
var sampledMinProb = minProbImage.sampleRegions({
  collection: biomassPoints,
  scale: 10,
  geometries: true
});
print('Sampled min probability at biomass points (preview):', sampledMinProb.limit(10));

// Aggregate: Find the lowest probability value found among all biomass samples
var minProbValue = sampledMinProb.aggregate_min('min_probability');
print('Overall minimum probability across all biomass points:', minProbValue);

// Filter biomass points below threshold
var belowThreshold = sampledMinProb.filter(ee.Filter.lt('min_probability', PROB_PRESENCE_THRESHOLD));
print('Biomass points below threshold:', belowThreshold.size());

Map.addLayer(belowThreshold, {color: 'red'}, 'Below Threshold Points');