// =============================================================================
// 11_distance_to_field_samples.js
// How far is each mapped seagrass pixel from the nearest AGC field sample?
//
// Output (two CSV files in Google Drive):
//   seagrass_area_by_distance_class.csv : area (ha) per region and distance class
//   seagrass_distance_percentiles.csv   : 10th/50th/90th percentile, mean and
//                                          pixel count of the distance per region
// Used for Figure 5 (R/figures_tables/figure4_figure5_distance.R) and Table S11.
//
// Assets: the same as the AGC model script (10_agc_model.js).
//   - AGC samples  = training_gse_YYYY (2018-2023), same cleaning
//                    and filters as the AGC model (AGC_pred > 0 and <= 2000)
//   - extent       = presence probability >= 0.8 in any year 2017-2024
//                    (presence_prob_YYYY), as in the Methods of the paper
//   - regions      = assigned from longitude/latitude (the six regions are
//                    400+ km apart), so no region polygon asset is needed
//   - AOI          = study_regions
// =============================================================================

// -----------------------------------------------------------------------------
// CONFIGURATION (same paths as the AGC model script)
// -----------------------------------------------------------------------------
var ASSET_ROOT   = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';   // your own GEE asset folder
var RESPONSE     = 'AGC_pred';
var TRAIN_YEARS  = [2018, 2019, 2020, 2021, 2022, 2023];
var TRAIN_PREFIX = ASSET_ROOT + '/training_gse_';

var EXTENT_YEARS   = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var PA_PROB_PREFIX = ASSET_ROOT + '/presence_prob_';
var PROB_THRESHOLD = 0.8;

var AOI = ee.FeatureCollection(
  ASSET_ROOT + '/study_regions'
).geometry();

var SCALE  = 10;
var FOLDER = 'GEE_exports';

// Distances are computed up to this radius. Pixels farther away get this value
// (they fall in the last class, and p90 may be capped at this value).
var SEARCH_RADIUS = 50000;

// -----------------------------------------------------------------------------
// 1. AGC field samples, all years (same cleaning as the AGC model script)
// -----------------------------------------------------------------------------
function cleanAGC(f) {
  var val = f.get(RESPONSE);
  var safeVal = ee.Algorithms.If(ee.Algorithms.IsEqual(val, 'NA'), -9999, val);
  return f.set(RESPONSE, ee.Number.parse(safeVal));
}

var samples = ee.FeatureCollection(TRAIN_YEARS.map(function (y) {
  return ee.FeatureCollection(TRAIN_PREFIX + y)
    .map(cleanAGC)
    .map(function (f) { return f.set('yr', y); });
})).flatten()
  .filter(ee.Filter.gt(RESPONSE, 0))
  .filter(ee.Filter.lte(RESPONSE, 2000));

// -----------------------------------------------------------------------------
// 2. Regions from longitude/latitude
//    1 Bintan (lon < 107), 2 Karimunjawa (107-115), 3 Komodo (115-121.5),
//    4 Rote (121.5-127, lat < -5), 5 Manado (121.5-127, lat >= -5),
//    6 Ayau (lon >= 127)
// -----------------------------------------------------------------------------
var REGION_NAMES = ee.List(['Bintan', 'Karimunjawa', 'Komodo', 'Rote', 'Manado', 'Ayau']);

function regionCode(lon, lat) {
  lon = ee.Number(lon);
  lat = ee.Number(lat);
  return ee.Number(ee.Algorithms.If(lon.lt(107), 1,
         ee.Algorithms.If(lon.lt(115), 2,
         ee.Algorithms.If(lon.lt(121.5), 3,
         ee.Algorithms.If(lon.lt(127),
           ee.Algorithms.If(lat.lt(-5), 4, 5),
           6)))));
}

samples = samples.map(function (f) {
  var xy = f.geometry().coordinates();
  var code = regionCode(xy.get(0), xy.get(1));
  return f.set('region_code', code, 'region', REGION_NAMES.get(code.subtract(1)));
});

var ll = ee.Image.pixelLonLat();
var lonImg = ll.select('longitude');
var latImg = ll.select('latitude');
var regionImg = ee.Image(6)
  .where(lonImg.lt(127).and(latImg.gte(-5)), 5)
  .where(lonImg.lt(127).and(latImg.lt(-5)), 4)
  .where(lonImg.lt(121.5), 3)
  .where(lonImg.lt(115), 2)
  .where(lonImg.lt(107), 1)
  .rename('region').toInt();

// Checks: compare these numbers with Table 2 (AGC samples per region)
print('Property names of one training point', ee.Feature(samples.first()).propertyNames());
print('AGC samples, all years (R analysis: 3,676)', samples.size());
print('AGC samples per region', samples.aggregate_histogram('region'));
print('AGC samples per year', samples.aggregate_histogram('yr'));

// -----------------------------------------------------------------------------
// 3. Seagrass extent: presence probability >= 0.8 in at least one year
// -----------------------------------------------------------------------------
var probs = ee.ImageCollection(EXTENT_YEARS.map(function (y) {
  return ee.Image(PA_PROB_PREFIX + y).rename('PA_prob');
}));
var extent = probs.max().gte(PROB_THRESHOLD).selfMask().clip(AOI);

// -----------------------------------------------------------------------------
// 4. Distance (m) from every pixel to the nearest AGC sample
// -----------------------------------------------------------------------------
var dist = samples.distance(SEARCH_RADIUS, 5)
  .unmask(SEARCH_RADIUS)
  .rename('dist_m');
var distSeagrass = dist.updateMask(extent);

// Distance classes (upper limits in m). Class k covers (limits[k-1], limits[k]].
var LIMITS = [10, 50, 100, 250, 500, 1000, 2000, 5000, SEARCH_RADIUS];
var LABELS = ['<=10 m', '10-50 m', '50-100 m', '100-250 m', '250-500 m',
              '500-1000 m', '1-2 km', '2-5 km', '>5 km'];
var cls = ee.Image(LIMITS.length - 1);
for (var k = LIMITS.length - 2; k >= 0; k--) {
  cls = cls.where(distSeagrass.lte(LIMITS[k]), k);
}
cls = cls.updateMask(distSeagrass.mask()).rename('cls').toInt();

// -----------------------------------------------------------------------------
// 5. Area (ha) per region and distance class
// -----------------------------------------------------------------------------
var areaHa = ee.Image.pixelArea().divide(10000).rename('area_ha');

var areaStack = areaHa.addBands(cls).addBands(regionImg).updateMask(extent);
var areaRes = areaStack.reduceRegion({
  reducer: ee.Reducer.sum()
    .group({groupField: 1, groupName: 'cls'})
    .group({groupField: 2, groupName: 'region_code'}),
  geometry: AOI, scale: SCALE, maxPixels: 1e13, tileScale: 4
});

var classRows = ee.FeatureCollection(
  ee.List(areaRes.get('groups')).map(function (o) {
    o = ee.Dictionary(o);
    var rc = ee.Number(o.get('region_code'));
    return ee.List(o.get('groups')).map(function (g) {
      g = ee.Dictionary(g);
      var c = ee.Number(g.get('cls'));
      return ee.Feature(null, {
        region: REGION_NAMES.get(rc.subtract(1)),
        region_code: rc,
        class_id: c,
        class_label: ee.List(LABELS).get(c),
        upper_m: ee.List(LIMITS).get(c),
        area_ha: g.get('sum')
      });
    });
  }).flatten()
);

// -----------------------------------------------------------------------------
// 6. Distance percentiles per region and for all regions together
// -----------------------------------------------------------------------------
var pctReducer = ee.Reducer.percentile([10, 50, 90])
  .combine({reducer2: ee.Reducer.mean(), sharedInputs: true})
  .combine({reducer2: ee.Reducer.count(), sharedInputs: true});

var pctRes = distSeagrass.addBands(regionImg).reduceRegion({
  reducer: pctReducer.group({groupField: 1, groupName: 'region_code'}),
  geometry: AOI, scale: SCALE, maxPixels: 1e13, tileScale: 4
});

var pctRows = ee.FeatureCollection(
  ee.List(pctRes.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    var rc = ee.Number(g.get('region_code'));
    return ee.Feature(null, g).set('region', REGION_NAMES.get(rc.subtract(1)));
  })
);

// All regions together: same grouped reducer with one constant group (0),
// so that the column names are the same as in the per-region rows
var pctAllRes = distSeagrass.addBands(ee.Image(0).rename('all').toInt()).reduceRegion({
  reducer: pctReducer.group({groupField: 1, groupName: 'region_code'}),
  geometry: AOI, scale: SCALE, maxPixels: 1e13, tileScale: 4
});
var pctAllRow = ee.Feature(null, ee.Dictionary(ee.List(pctAllRes.get('groups')).get(0)))
  .set('region', 'All');
pctRows = pctRows.merge(ee.FeatureCollection([pctAllRow]));

// -----------------------------------------------------------------------------
// 7. Map and exports. The area and percentile tables are only computed in the
//    export tasks (printing them here may time out).
// -----------------------------------------------------------------------------
Map.centerObject(samples, 5);
Map.addLayer(extent, {palette: ['00ff00']}, 'seagrass extent (prob >= 0.8 any year)', false);
Map.addLayer(distSeagrass, {min: 0, max: 2000, palette: ['0571b0', 'f7f7f7', 'ca0020']},
             'distance to nearest AGC sample (m)');
Map.addLayer(samples, {color: 'black'}, 'AGC samples');

Export.table.toDrive({
  collection: classRows,
  description: 'seagrass_area_by_distance_class',
  folder: FOLDER, fileFormat: 'CSV',
  selectors: ['region', 'region_code', 'class_id', 'class_label', 'upper_m', 'area_ha']
});
Export.table.toDrive({
  collection: pctRows,
  description: 'seagrass_distance_percentiles',
  folder: FOLDER, fileFormat: 'CSV'
});
