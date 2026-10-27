/*****************************************************************************
 * 01_extract_gse_at_field_points.js
 * Samples the 64 GSE bands, depth and distance to land (plus slope, rugosity
 * and wave variables, which are not used in the models) at every field point,
 * and exports one CSV per year (gse_points_<year>.csv). These files are the
 * inputs of R/workflow/01_prepare_point_data.R.
 *
 * Inputs (ASSET_ROOT)  field_points, bathymetry, wave_2024
 *
 * Set ASSET_ROOT below to your own GEE asset folder. All paths are built
 * from it, except the shared Allen Coral Atlas distance-to-land layer,
 * which needs its own access.
 ****************************************************************************/

// ------------------ [A] INPUTS ------------------
var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var TRAINING_FC_PATH   = ASSET_ROOT + '/field_points';
var DEPTH_ASSET        = ASSET_ROOT + '/bathymetry';
// Distance to land (shared Allen Coral Atlas asset, needs its own access)
var DIST_TO_LAND_ASSET = 'projects/coral_atlas/global_datasets/osm_distToLand_indo';
var WAVE_ASSET_2024    = ASSET_ROOT + '/wave_2024';

// Indonesia bounding-box geometry (plain coordinates, not a private asset)
var indo_area = ee.Geometry.Polygon(
  [[[94.76922917931722, 7.317170170331969],
    [94.76922917931722, -11.309399821791997],
    [141.17547917931722, -11.309399821791997],
    [141.17547917931722, 7.317170170331969]]], null, false);

var INDO_EXTENT_GEOM   = indo_area;
var DEPTH_POSITIVE     = true; // Set true if depth values are positive downward
var YEARS              = ee.List.sequence(2017, 2024); // years to export

// ------------------ [B] LOAD BASE DATA ------------------
var KEY = 'gee_id';

var training_points = ee.FeatureCollection(TRAINING_FC_PATH)
  .map(function(f){ return f.set(KEY, ee.String(f.get(KEY))); });

var depth = ee.Image(DEPTH_ASSET).rename('depth');
if (DEPTH_POSITIVE) {
  depth = depth.multiply(-1); // ensure negative depth values
}
var distToLand = ee.Image(DIST_TO_LAND_ASSET).rename('distToLand');

// slope and rugosity from depth
function getSlopeAndRugosity(depthImg) {
  var slope = ee.Terrain.slope(depthImg).rename("slope");
  var rugosity = slope.reduceNeighborhood({
    reducer: ee.Reducer.stdDev(),
    kernel: ee.Kernel.circle(100, 'meters')
  }).unmask(0).rename("rugosity");
  return slope.addBands(rugosity);
}

// ------------------ [B.2] LOAD WAVE DATA ------------------
// ERA5-GEBCO static multiband raster:
// Band 1 = elevation (m), Band 2 = significant wave height (m), Band 3 = mean wave period (s)

var wave_raw = ee.Image(WAVE_ASSET_2024)
  .rename(['elevation', 'sig_wave_height', 'mean_wave_period']);

// Reproject to 10m resolution to match other layers
var wave_reprojected = wave_raw.reproject({
  crs: 'EPSG:4326',
  scale: 10
});

// ------------------ [C] GSE Builder Function ------------------
function buildGSE(year){
  var ic = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
    .filterDate(ee.Date.fromYMD(year,1,1), ee.Date.fromYMD(year+1,1,1))
    .filterBounds(INDO_EXTENT_GEOM);
  var img = ee.Image(ic.mosaic()).unmask(-9999, false);
  return img.rename(img.bandNames().map(function(b){ return ee.String('GSE_').cat(b); }));
}

// ------------------ [D] LOOP OVER YEARS ------------------
YEARS.getInfo().forEach(function(year) {
  var gseImg = buildGSE(year);
  var slopeRugosity = getSlopeAndRugosity(depth);

  // Build raster stack (all years include static layers)
  var stack = gseImg
    .addBands(depth)
    .addBands(distToLand)
    .addBands(slopeRugosity)
    .addBands(wave_reprojected);

  // Extract pixel values at each training point (preserve all rows)
  var sampled = training_points.map(function(point) {
    var values = stack.reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: point.geometry(),
      scale: 10,
      maxPixels: 1e13
    });

    var lonlat = point.geometry().coordinates();

    return point.set(values)
      .set('lon', lonlat.get(0))
      .set('lat', lonlat.get(1));
  });

  // Count total and valid GSE samples
  var total = sampled.size();
  var valid = sampled.filter(ee.Filter.notNull(['GSE_A01'])).size();

  print('Year:', year);
  print('Total training points:', total);
  print('With valid GSE_A01 values:', valid);

  // Export to CSV
  Export.table.toDrive({
    collection: sampled,
    description: 'gse_points_' + year,
    folder: 'GEE_exports',
    fileNamePrefix: 'gse_points_' + year,
    fileFormat: 'CSV'
  });

  print('Export scheduled for year', year);
});

// ------------------ [E] MAP QA (for last year only) ------------------
var gseLast = buildGSE(2024);
Map.centerObject(INDO_EXTENT_GEOM, 5);

Map.addLayer(
  gseLast.select(['GSE_A01','GSE_A16','GSE_A09']).rename(['R','G','B']),
  {min: -0.3, max: 0.3},
  'GSE RGB 2024', true
);

Map.addLayer(depth, {min: -30, max: 0, palette: ['blue', 'lightblue']}, 'Depth (neg)', false);
Map.addLayer(distToLand, {min: 0, max: 100000, palette: ['green','yellow','red']}, 'Dist to Land', false);
Map.addLayer(training_points, {color: 'FF8800'}, 'Training Points', true);

// ------------------ [F] DEBUG PROJECTIONS ------------------
print('Depth projection:', depth.projection());
print('GSE projection:', gseLast.projection());
print('DistToLand projection:', distToLand.projection());
print('Wave projection:', wave_reprojected.projection());
print('Total training points:', training_points.size());