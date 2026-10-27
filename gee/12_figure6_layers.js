// =============================================================================
// 12_figure6_layers.js
// Layers for Figure 6, for one region and one survey year (Rote 2021 in the paper,
// the survey with the most AGC samples).
//   (a) AGC predicted by the proxy workflow (10_agc_model.js)
//   (b) distance from each mapped pixel to the nearest AGC sample of that year
//   (c) expected accuracy class from this distance:
//         1  <= 10 m
//         2  10-100 m
//         3  100-250 m
//         4  250 m-1 km
//         5  > 1 km         untested within the survey
//       Each class takes the MEC of the pixel or block design whose median
//       distance in this region falls in the class (Table S8).
//       Only samples of the same year are used, because samples of other years
//       correspond to the new-survey design.
//
// Uses the AGC maps of 10_agc_model.js (agc_YYYY).
// Output (Google Drive, folder GEE_exports), with REGION and YEAR in the names:
//   Fig6_<REGION><YEAR>_layers.tif   bands AGC, dist_m, acc_class (10 m, UTM)
//   Fig6_<REGION><YEAR>_S2_RGB.tif   Sentinel-2 true colour of YEAR (20 m)
//   Fig6_<REGION><YEAR>_samples.csv  lon, lat, AGC of the samples of YEAR
// The console prints the mapped area and mean AGC per accuracy class.
// =============================================================================

// -----------------------------------------------------------------------------
// CONFIGURATION  (only this block changes between regions)
// -----------------------------------------------------------------------------
var YEAR       = 2021;
var REGION     = 'Rote';
var CRS        = 'EPSG:32751';   // UTM zone 51S (Rote, about 123 E). Komodo: EPSG:32750
var N_EXPECTED = '1,461';        // AGC samples of Rote 2021 (new-survey test). Komodo 2019: 1,240

var RESPONSE = 'AGC_pred';
var ASSET_ROOT     = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';   // your own GEE asset folder
var OUT            = ASSET_ROOT + '/';
var AGC_ASSET      = OUT + 'agc_' + YEAR;
var TRAIN_PREFIX   = OUT + 'training_gse_';
var PA_PROB_PREFIX = OUT + 'presence_prob_';
var NAME           = 'Fig6_' + REGION + YEAR;

var AOI = ee.FeatureCollection(
  ASSET_ROOT + '/study_regions'
).filter(ee.Filter.eq('loc', REGION)).geometry();

var SCALE  = 10;
var FOLDER = 'GEE_exports';
var SEARCH_RADIUS = 50000;

// -----------------------------------------------------------------------------
// 1. Mapped extent: presence probability >= 0.8 in any year 2017-2024
//    (same extent as 05_seagrass_extent.js and Table S11)
// -----------------------------------------------------------------------------
var EXTENT_YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var extent = ee.ImageCollection(EXTENT_YEARS.map(function (y) {
  return ee.Image(PA_PROB_PREFIX + y).rename('PA_prob');
})).max().gte(0.8).selfMask().clip(AOI);

// -----------------------------------------------------------------------------
// 2. AGC samples of the same year (same cleaning as the AGC model script)
// -----------------------------------------------------------------------------
function cleanAGC(f) {
  var val = f.get(RESPONSE);
  var safeVal = ee.Algorithms.If(ee.Algorithms.IsEqual(val, 'NA'), -9999, val);
  return f.set(RESPONSE, ee.Number.parse(safeVal));
}
var samples = ee.FeatureCollection(TRAIN_PREFIX + YEAR)
  .filterBounds(AOI)
  .map(cleanAGC)
  .filter(ee.Filter.gt(RESPONSE, 0))
  .filter(ee.Filter.lte(RESPONSE, 2000));
print('CHECK AGC samples, ' + REGION + ' ' + YEAR + ' (expected ' + N_EXPECTED + ')', samples.size());

// -----------------------------------------------------------------------------
// 3. Layers
// -----------------------------------------------------------------------------
var agc = ee.Image(AGC_ASSET).rename('AGC').updateMask(extent);

var dist = samples.distance(SEARCH_RADIUS, 5)
  .unmask(SEARCH_RADIUS)
  .rename('dist_m')
  .updateMask(extent);

var accClass = ee.Image(5)
  .where(dist.lte(1000), 4)
  .where(dist.lte(250), 3)
  .where(dist.lte(100), 2)
  .where(dist.lte(10), 1)
  .updateMask(extent)
  .rename('acc_class');

// -----------------------------------------------------------------------------
// 4. Mapped area (ha) and mean AGC per accuracy class (for the text)
// -----------------------------------------------------------------------------
var areaHa = ee.Image.pixelArea().divide(10000).rename('area_ha').updateMask(extent);
print('Mapped area (ha) per accuracy class', areaHa.addBands(accClass).reduceRegion({
  reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'acc_class'}),
  geometry: AOI, scale: SCALE, maxPixels: 1e13, tileScale: 4
}));
print('Mean AGC (gC m-2) per accuracy class', agc.addBands(accClass).reduceRegion({
  reducer: ee.Reducer.mean().group({groupField: 1, groupName: 'acc_class'}),
  geometry: AOI, scale: SCALE, maxPixels: 1e13, tileScale: 4
}));

// -----------------------------------------------------------------------------
// 5. Sentinel-2 true colour of the same year, as background
// -----------------------------------------------------------------------------
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI)
  .filterDate(ee.Date.fromYMD(YEAR, 1, 1), ee.Date.fromYMD(YEAR + 1, 1, 1))
  .map(function (img) {
    var scl = img.select('SCL');
    var ok = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
    return img.updateMask(ok).select(['B4', 'B3', 'B2']);
  })
  .reduce(ee.Reducer.percentile([60]))
  .rename(['B4', 'B3', 'B2']);
var rgb = s2.visualize({min: 0, max: 2500}).clip(AOI);

// -----------------------------------------------------------------------------
// 6. Preview and exports
// -----------------------------------------------------------------------------
Map.centerObject(AOI, 11);
Map.addLayer(rgb, {}, 'Sentinel-2 ' + YEAR, false);
Map.addLayer(agc, {min: 0, max: 80, palette: ['f7fcf5', 'a1d99b', '31a354', '00441b']}, 'AGC ' + YEAR);
Map.addLayer(dist, {min: 0, max: 2000, palette: ['0571b0', 'f7f7f7', 'ca0020']}, 'distance (m)', false);
Map.addLayer(accClass, {min: 1, max: 5, palette: ['08306b', '2563a8', '4a90c8', '8fc1e3', 'd9711c']},
             'expected accuracy class', false);
Map.addLayer(samples, {color: 'black'}, 'AGC samples ' + YEAR, false);

Export.image.toDrive({
  image: agc.toFloat().addBands(dist.toFloat()).addBands(accClass.toFloat()),
  description: NAME + '_layers',
  folder: FOLDER, region: AOI, scale: SCALE, crs: CRS,
  maxPixels: 1e13, fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: rgb,
  description: NAME + '_S2_RGB',
  folder: FOLDER, region: AOI, scale: 20, crs: CRS,
  maxPixels: 1e13, fileFormat: 'GeoTIFF'
});
Export.table.toDrive({
  collection: samples.map(function (f) {
    var xy = f.geometry().coordinates();
    return ee.Feature(null, {lon: xy.get(0), lat: xy.get(1), AGC: f.get(RESPONSE)});
  }),
  description: NAME + '_samples',
  folder: FOLDER, fileFormat: 'CSV', selectors: ['lon', 'lat', 'AGC']
});
