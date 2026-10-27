/*****************************************************************************
 * 05_seagrass_extent.js
 * Mapped seagrass extent of the paper: pixels with a presence probability of
 * at least 0.8 in any year from 2017 to 2024 (maximum of the yearly
 * presence_prob_<year> layers of 04_presence_model.js). The script also prints
 * and exports the mapped area of each region.
 *
 * Inputs  (ASSET_ROOT)  presence_prob_<year>, study_regions
 * Outputs (ASSET_ROOT)  presence_prob_max_2017_2024, seagrass_extent_2017_2024
 *         (Drive)       Seagrass_MaxProb_0_8_2017_2024_area_per_location_km2.csv
 *
 * Set ASSET_ROOT below to your own GEE asset folder.
 ****************************************************************************/

// ==========================================================
// CONFIGURATION
// ==========================================================
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var START_YEAR = 2017;
var END_YEAR = 2024;

var EXPORT_BASE = ASSET_ROOT;
var PROB_ASSET_PREFIX = EXPORT_BASE + '/presence_prob_';

var AOI = ee.FeatureCollection(ASSET_ROOT + '/study_regions').geometry();
Map.addLayer(AOI, {color: 'yellow'}, 'AOI', false);

// Regions (property 'loc'), for the area per region
var locations = ee.FeatureCollection(ASSET_ROOT + '/study_regions');

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
// MAXIMUM PROBABILITY ACROSS YEARS (per pixel)
// Pixels predicted as seagrass (probability >= 0.8) in any year 2017-2024
// ==========================================================

// === CONFIGURATION ===
var MAX_PROB_THRESHOLD = 0.8;  // Probability threshold for defining seagrass presence
var THRESHOLD_STR = String(MAX_PROB_THRESHOLD).replace('.', '_');  // e.g., 0.8 -> "0_8"
var EXPORT_PREFIX = 'Seagrass_MaxProb_' + THRESHOLD_STR + '_2017_2024';  // name of the area table

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
}, 'Maximum Probability (2017-2024)');

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
  assetId: EXPORT_BASE + '/presence_prob_max_2017_2024',
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
  assetId: EXPORT_BASE + '/seagrass_extent_2017_2024',
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
