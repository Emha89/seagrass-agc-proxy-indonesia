/************************************************************
 * 3rd_paper_Apps.js
 * Interactive GEE App: Seagrass Carbon Proxy Viewer.
 *
 * Lets a reader browse all 8 proxy-chain layers (PA probability,
 * 3-class morphology, percent cover, AGB, carbon index, AGC) for any of
 * the 6 study regions and any year 2017-2024, click a pixel to see all
 * proxy values and an AGC time series at that location, and see the
 * region-level AGC total + 95% CI for the selected region/year.
 *
 * PIPELINE POSITION: reads the same PA_PROB_PREFIX / MORPH_PREFIX /
 * SPC_PREFIX / AGB_PREFIX / CINDEX_PREFIX / AGC_PREFIX assets produced
 * by 02_PROB_rf_modelDev.js through 07_AGC_rf_modelDev.js (all date
 * stamps match exactly). The AGC_TABLE dictionary below (total AGC + CI
 * per region per year) is the accumulated printed output of running
 * 07_AGC_rf_modelUncer.js once per region-year combination (6 regions x
 * 8 years = 48 runs) -- not computed live by this script.
 *
 * The errorInfo strings per layer (F1/accuracy/RMSE/R2 with 95% CI) are
 * likewise pre-computed summary statistics, hardcoded here for display
 * rather than recomputed live.
 *
 * FIX: "Limited within persitence seagrass area." corrected to
 * "persistence".
 *
 * ASSET PORTABILITY: this script originally referenced the author's own
 * private GEE assets, all under one project/folder root -- replace
 * ASSET_ROOT below with your own GEE project and asset folder.
 ************************************************************/

var ASSET_ROOT = 'projects/YOUR-GEE-PROJECT/assets/YOUR-FOLDER';

var PA_PROB_PREFIX = ASSET_ROOT + '/output/RF_probability04022026_';
var MORPH_PREFIX   = ASSET_ROOT + '/output/MORPH3_probs04022026_';
var SPC_PREFIX     = ASSET_ROOT + '/output/RF_tSPC04022026_';
var AGB_PREFIX     = ASSET_ROOT + '/output/RF_AGB04022026_';
var CINDEX_PREFIX  = ASSET_ROOT + '/output/RF_cIndex04022026_';
var AGC_PREFIX     = ASSET_ROOT + '/output/RF_AGC04022026_';

var PERSISTENCE_MASK =
  ASSET_ROOT + '/output/Seagrass_MaxProb_0_8_2017_2024_mask';

var AOI_ASSET = ASSET_ROOT + '/InputArea/R2_case_study';

var SCALE = 10;

var SITES = {
  'Ayau':        [131.04647,  0.36225, 15],
  'Bintan':      [104.56868,  1.23176, 15],
  'Karimunjawa': [110.47921, -5.77212, 15],
  'Komodo':      [119.72583, -8.56826, 15],
  'Manado':      [124.93688,  1.68922, 15],
  'Rote':        [122.81111,-10.77786, 15]
};

var PAL_PA     = ['#ffffff','#fee0d2','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d'];
var PAL_MORPH  = ['#ffffff','#efedf5','#d9d8ea','#bcbddc','#9e9ac8','#807dba','#6a51a3','#4a1486'];
var PAL_SPC    = ['#ffffff','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b'];
var PAL_AGB    = ['#ffffff','#f5e4c3','#e8c98a','#d4a44c','#b8832a','#8c6020','#5c3d10','#3b2407'];
var PAL_CINDEX = ['#ffffff','#fde0ef','#fbb4d4','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'];
var PAL_AGC    = ['#9e0142','#d6404e','#f57547','#fee08b','#cdeab5','#84cdb4','#5b9ec9','#3288bd'];

var LAYERS = [
  {
    label: 'PA Probability', short: 'PA prob',
    getImage: function(y){ return ee.Image(PA_PROB_PREFIX + y).rename('PA_prob'); },
    bandName: 'PA_prob', unit: '0-1', min: 0, max: 1, palette: PAL_PA,
    errorInfo: 'Model F1: 0.99  (95% CI: 0.98-0.99)'
  },
  {
    label: 'Morphology - Mono Enhalus (P_mono_Ea)', short: 'Mono Ea',
    getImage: function(y){ return ee.Image(MORPH_PREFIX + y).select('P_mono_Ea'); },
    bandName: 'P_mono_Ea', unit: '0-1', min: 0, max: 1, palette: PAL_MORPH,
    errorInfo: 'Model Accuracy: 0.87  (95% CI: 0.86-0.88)'
  },
  {
    label: 'Morphology - Mixed Long (P_mixed_long)', short: 'Mixed Long',
    getImage: function(y){ return ee.Image(MORPH_PREFIX + y).select('P_mixed_long'); },
    bandName: 'P_mixed_long', unit: '0-1', min: 0, max: 1, palette: PAL_MORPH,
    errorInfo: 'Model Accuracy: 0.87  (95% CI: 0.86-0.88)'
  },
  {
    label: 'Morphology - Mixed Short (P_mixed_short)', short: 'Mixed Short',
    getImage: function(y){
      return ee.Image(MORPH_PREFIX + y).select('P_mixed_short_plus_mono_short');
    },
    bandName: 'P_mixed_short_plus_mono_short', unit: '0-1', min: 0, max: 1, palette: PAL_MORPH,
    errorInfo: 'Model Accuracy: 0.87  (95% CI: 0.86-0.88)'
  },
  {
    label: 'Seagrass % Cover (SPC)', short: 'SPC',
    getImage: function(y){ return ee.Image(SPC_PREFIX + y).rename('tSPC_pred'); },
    bandName: 'tSPC_pred', unit: '%', min: 0, max: 100, palette: PAL_SPC,
    errorInfo: 'Model RMSE: 16.97%  (95% CI: 16.04-17.80)  R\u00b2: 0.57'
  },
  {
    label: 'Above-Ground Biomass (AGB)', short: 'AGB',
    getImage: function(y){ return ee.Image(AGB_PREFIX + y).rename('tAGB_pred'); },
    bandName: 'tAGB_pred', unit: 'gDW/m2', min: null, max: null, palette: PAL_AGB,
    errorInfo: 'Model RMSE: 55.56 gDW/m\u00b2  (95% CI: 54.80-56.33)  R\u00b2: 0.49'
  },
  {
    label: 'Carbon Index', short: 'C-Index',
    getImage: function(y){ return ee.Image(CINDEX_PREFIX + y).rename('carbon_index'); },
    bandName: 'carbon_index', unit: '0-1', min: 0, max: 1, palette: PAL_CINDEX,
    errorInfo: 'Model RMSE: 0.10  (95% CI: 0.09-0.10)  R\u00b2: 0.67'
  },
  {
    label: 'Above-Ground Carbon (AGC)', short: 'AGC',
    getImage: function(y){ return ee.Image(AGC_PREFIX + y).rename('AGC_pred'); },
    bandName: 'AGC_pred', unit: 'gC/m2', min: null, max: null, palette: PAL_AGC,
    errorInfo: 'Model RMSE: 13.56 gC/m\u00b2  (95% CI: 13.35-13.75)  R\u00b2: 0.54'
  }
];

var LAYER_MAP = {};
for (var _i = 0; _i < LAYERS.length; _i++) {
  LAYER_MAP[LAYERS[_i].label] = LAYERS[_i];
}
var LAYER_LABELS = LAYERS.map(function(l){ return l.label; });

var persistenceMask = ee.Image(PERSISTENCE_MASK);
var aoi = ee.FeatureCollection(AOI_ASSET).geometry();

ui.root.clear();

var mainMap = ui.Map();
mainMap.setOptions('HYBRID');
mainMap.setControlVisibility({ zoomControl: true, mapTypeControl: false });
mainMap.style().set('cursor', 'crosshair');

var controlPanel = ui.Panel({
  style: { width: '310px', padding: '10px 12px',
           backgroundColor: '#16213e' }
});

ui.root.add(ui.SplitPanel({
  firstPanel:  controlPanel,
  secondPanel: mainMap,
  orientation: 'horizontal',
  wipe: false
}));

function mkLabel(text, style) {
  return ui.Label({ value: text, style: style });
}
function divider() {
  return mkLabel('------------------------------',
    { color: '#2e4070', backgroundColor: '#16213e', margin: '5px 0' });
}
function sectionHead(text) {
  return mkLabel(text, {
    fontWeight: 'bold', color: '#a0c4ff',
    backgroundColor: '#16213e', fontSize: '12px', margin: '5px 0 2px 0'
  });
}

controlPanel.add(mkLabel('Seagrass Carbon Proxy Viewer', {
  fontSize: '15px', fontWeight: 'bold',
  color: '#e8f4fd', backgroundColor: '#16213e', margin: '0 0 2px 0'
}));
controlPanel.add(mkLabel(
  'Annual proxy-layer maps (2017-2024) across 6 Indonesian coastal regions.', {
  fontSize: '11px', color: '#8899bb',
  backgroundColor: '#16213e', margin: '0 0 4px 0', whiteSpace: 'wrap'
}));
controlPanel.add(divider());

controlPanel.add(sectionHead('Study Region'));
var siteSelect = ui.Select({
  items: Object.keys(SITES),
  value: 'Ayau',
  style: { width: '282px', margin: '2px 0 8px 0' },
  onChange: function(site) {
    var si = SITES[site];
    mainMap.setCenter(si[0], si[1], si[2]);
    updateLayer();
    updateAGCSummary();
  }
});
controlPanel.add(siteSelect);

controlPanel.add(sectionHead('Year'));
var yearDisp = mkLabel('2024', {
  color: '#ffffff', backgroundColor: '#16213e',
  fontSize: '22px', fontWeight: 'bold', margin: '0 0 2px 0'
});
controlPanel.add(yearDisp);

var yearSlider = ui.Slider({
  min: 2017, max: 2024, value: 2024, step: 1,
  style: { width: '282px', margin: '0 0 8px 0' },
  onChange: function(v) {
    yearDisp.setValue(String(Math.round(v)));
    updateLayer();
    updateAGCSummary();
  }
});
controlPanel.add(yearSlider);

controlPanel.add(sectionHead('Proxy Layer'));
var layerSelect = ui.Select({
  items: LAYER_LABELS,
  value: 'Above-Ground Carbon (AGC)',
  style: { width: '282px', margin: '2px 0 10px 0' },
  onChange: function() { updateLayer(); updateErrorInfo(); }
});
controlPanel.add(layerSelect);

controlPanel.add(divider());

var legendTitle = mkLabel('', {
  fontWeight: 'bold', color: '#a0c4ff',
  backgroundColor: '#16213e', fontSize: '11px',
  margin: '2px 0 4px 0', whiteSpace: 'wrap'
});
var colorBarImg = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select('longitude')
           .visualize({ min: 0, max: 1, palette: PAL_AGC }),
  params: { bbox: [0, 0, 1, 0.1], dimensions: '270x18' },
  style: { stretch: 'horizontal', margin: '0 0 2px 0', maxHeight: '18px' }
});
var legendMinLbl  = mkLabel('0', {
  color: '#aabbcc', backgroundColor: '#16213e', fontSize: '10px'
});
var legendUnitLbl = mkLabel('', {
  color: '#aabbcc', backgroundColor: '#16213e', fontSize: '10px',
  textAlign: 'center', stretch: 'horizontal'
});
var legendMaxLbl  = mkLabel('1', {
  color: '#aabbcc', backgroundColor: '#16213e', fontSize: '10px'
});
controlPanel.add(legendTitle);
controlPanel.add(colorBarImg);
controlPanel.add(ui.Panel({
  widgets: [legendMinLbl, legendUnitLbl, legendMaxLbl],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: { backgroundColor: '#16213e', stretch: 'horizontal' }
}));

function updateLegend(cfg, lo, hi) {
  legendTitle.setValue('Legend: ' + cfg.label + ' (' + cfg.unit + ')');
  colorBarImg.setImage(
    ee.Image.pixelLonLat().select('longitude')
      .multiply(hi - lo).add(lo)
      .visualize({ min: lo, max: hi, palette: cfg.palette })
  );
  colorBarImg.setParams({ bbox: [0, 0, 1, 0.1], dimensions: '270x18' });
  legendMinLbl.setValue(lo.toFixed(2));
  legendMaxLbl.setValue(hi.toFixed(2));
  legendUnitLbl.setValue(cfg.unit);
}

controlPanel.add(divider());

var inspectorStatus = mkLabel('', {
  color: '#ffd580', backgroundColor: '#16213e',
  fontSize: '10px', margin: '0 0 2px 0'
});
controlPanel.add(inspectorStatus);

var errorInfoLabel = mkLabel('', {
  color: '#88aacc', backgroundColor: '#16213e',
  fontSize: '10px', whiteSpace: 'wrap',
  margin: '0 0 2px 0', fontStyle: 'italic'
});
controlPanel.add(errorInfoLabel);

function updateErrorInfo() {
  var cfg = LAYER_MAP[layerSelect.getValue()];
  errorInfoLabel.setValue('Monte Carlo model error: ' + cfg.errorInfo);
}

controlPanel.add(divider());

controlPanel.add(sectionHead('Total AGC -- Regional Estimate'));
controlPanel.add(mkLabel(
  'Limited within persistence seagrass area.',
  { color: '#8899bb', backgroundColor: '#16213e',
    fontSize: '10px', whiteSpace: 'wrap', margin: '0 0 4px 0' }
));

var agcSummaryLabel = mkLabel('', {
  color: '#ffffff', backgroundColor: '#1e2f50',
  fontSize: '11px', whiteSpace: 'wrap',
  margin: '0', padding: '6px',
  border: '1px solid #2e4070'
});
controlPanel.add(agcSummaryLabel);

// Accumulated output of running 07_AGC_rf_modelUncer.js once per
// region-year combination (see header note above) -- not computed live.
var AGC_TABLE = {
  'Ayau': {
    2017: [38.84, 38.04, 39.64, 0.80, 0.80],
    2018: [42.25, 41.45, 43.05, 0.80, 0.80],
    2019: [43.98, 43.18, 44.78, 0.80, 0.80],
    2020: [39.64, 38.84, 40.44, 0.80, 0.80],
    2021: [35.41, 34.61, 36.21, 0.80, 0.80],
    2022: [35.52, 34.72, 36.32, 0.80, 0.80],
    2023: [38.91, 38.10, 39.71, 0.80, 0.80],
    2024: [38.29, 37.49, 39.10, 0.80, 0.80]
  },
  'Bintan': {
    2017: [2169.74, 2163.63, 2175.85, 6.11, 6.11],
    2018: [2464.40, 2458.27, 2470.53, 6.13, 6.13],
    2019: [2625.14, 2619.02, 2631.25, 6.11, 6.11],
    2020: [2291.90, 2285.80, 2298.00, 6.10, 6.10],
    2021: [2220.94, 2214.84, 2227.04, 6.10, 6.10],
    2022: [2046.83, 2040.73, 2052.93, 6.10, 6.10],
    2023: [2357.04, 2350.93, 2363.14, 6.11, 6.11],
    2024: [2374.92, 2368.81, 2381.03, 6.11, 6.11]
  },
  'Karimunjawa': {
    2017: [670.20, 666.67, 673.73, 3.53, 3.53],
    2018: [856.13, 852.58, 859.67, 3.54, 3.54],
    2019: [867.77, 864.23, 871.30, 3.53, 3.53],
    2020: [645.28, 641.75, 648.80, 3.53, 3.53],
    2021: [613.85, 610.32, 617.37, 3.53, 3.53],
    2022: [532.37, 528.84, 535.89, 3.52, 3.52],
    2023: [790.79, 787.26, 794.31, 3.53, 3.53],
    2024: [824.58, 821.05, 828.11, 3.53, 3.53]
  },
  'Komodo': {
    2017: [3594.06, 3586.49, 3601.62, 7.57, 7.57],
    2018: [4198.87, 4191.31, 4206.43, 7.56, 7.56],
    2019: [4094.56, 4087.00, 4102.12, 7.56, 7.56],
    2020: [3521.25, 3513.68, 3528.82, 7.57, 7.57],
    2021: [3414.31, 3406.74, 3421.88, 7.57, 7.57],
    2022: [3361.69, 3354.12, 3369.26, 7.57, 7.57],
    2023: [3934.42, 3926.85, 3941.99, 7.57, 7.57],
    2024: [4096.91, 4089.35, 4104.48, 7.56, 7.56]
  },
  'Manado': {
    2017: [768.03, 763.93, 772.13, 4.10, 4.10],
    2018: [816.86, 812.77, 820.95, 4.09, 4.09],
    2019: [892.98, 888.89, 897.06, 4.09, 4.09],
    2020: [790.59, 786.50, 794.68, 4.09, 4.09],
    2021: [735.26, 731.16, 739.35, 4.09, 4.09],
    2022: [744.90, 740.81, 748.99, 4.09, 4.09],
    2023: [804.66, 800.56, 808.75, 4.09, 4.09],
    2024: [800.49, 796.40, 804.58, 4.09, 4.09]
  },
  'Rote': {
    2017: [2423.08, 2416.85, 2429.31, 6.23, 6.23],
    2018: [2737.89, 2731.66, 2744.11, 6.23, 6.23],
    2019: [2686.20, 2679.97, 2692.43, 6.23, 6.23],
    2020: [2556.37, 2550.14, 2562.60, 6.23, 6.23],
    2021: [2325.59, 2319.36, 2331.81, 6.23, 6.23],
    2022: [2370.38, 2364.15, 2376.61, 6.23, 6.23],
    2023: [2625.08, 2618.85, 2631.30, 6.23, 6.23],
    2024: [2506.89, 2500.66, 2513.11, 6.23, 6.23]
  }
};

function updateAGCSummary() {
  var site = siteSelect.getValue();
  var year = Math.round(yearSlider.getValue());
  var row  = AGC_TABLE[site] ? AGC_TABLE[site][year] : null;
  if (!row) { agcSummaryLabel.setValue('No data available.'); return; }
  agcSummaryLabel.setValue(
    site + '  |  ' + year + '\n' +
    'Total AGC :  ' + row[0].toFixed(2) + ' ton C\n' +
    '95% CI    :  ' + row[1].toFixed(2) + ' \u2013 ' + row[2].toFixed(2) + ' ton C\n' +
    'Error     :  \u2212' + row[3].toFixed(2) + '  /  +' + row[4].toFixed(2) + ' ton C'
  );
}

controlPanel.add(divider());

controlPanel.add(sectionHead('AGC Time-Series at Clicked Pixel'));
controlPanel.add(mkLabel(
  'Click a pixel on the map to plot AGC (2017-2024) at that location.',
  { color: '#8899bb', backgroundColor: '#16213e',
    fontSize: '10px', whiteSpace: 'wrap', margin: '0 0 4px 0' }
));

var chartPanel = ui.Panel({
  style: { backgroundColor: '#16213e', margin: '0', padding: '0' }
});
controlPanel.add(chartPanel);

function buildAGCChart(pt, lat, lon) {
  chartPanel.clear();
  chartPanel.add(mkLabel('Building chart...', {
    color: '#ffd580', backgroundColor: '#16213e', fontSize: '10px'
  }));

  var allYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];

  var agcStack = ee.Image(
    allYears.map(function(y) {
      return ee.Image(AGC_PREFIX + y)
               .updateMask(persistenceMask)
               .rename('AGC_' + y);
    })
  );

  agcStack.reduceRegion({
    reducer:   ee.Reducer.first(),
    geometry:  pt,
    scale:     SCALE,
    maxPixels: 1e6
  }).evaluate(function(result) {
    chartPanel.clear();
    if (!result) {
      chartPanel.add(mkLabel('No AGC data at this pixel.', {
        color: '#ff8888', backgroundColor: '#16213e', fontSize: '10px'
      }));
      return;
    }
    var vals = [];
    var hasData = false;
    for (var yi = 0; yi < allYears.length; yi++) {
      var v = result['AGC_' + allYears[yi]];
      vals.push((v !== null && v !== undefined) ? v : null);
      if (v !== null && v !== undefined) hasData = true;
    }
    if (!hasData) {
      chartPanel.add(mkLabel('Pixel outside seagrass mask or AOI.', {
        color: '#ff8888', backgroundColor: '#16213e', fontSize: '10px'
      }));
      return;
    }
    var features = ee.FeatureCollection(
      allYears.map(function(y) {
        var v = result['AGC_' + y];
        return ee.Feature(null, {
          year: y,
          AGC:  (v !== null && v !== undefined) ? v : -9999
        });
      })
    );
    var chart = ui.Chart.feature.byFeature({
      features: features, xProperty: 'year', yProperties: ['AGC']
    })
    .setChartType('LineChart')
    .setOptions({
      title: 'AGC  (lat ' + lat.toFixed(4) + ', lon ' + lon.toFixed(4) + ')',
      titleTextStyle: { color: '#cce0ff', fontSize: 10, bold: true },
      hAxis: {
        title: 'Year',
        titleTextStyle: { color: '#8899bb', fontSize: 9 },
        textStyle: { color: '#aabbcc', fontSize: 9 },
        gridlines: { color: '#2e4070' }, format: '####'
      },
      vAxis: {
        title: 'AGC (gC/m2)',
        titleTextStyle: { color: '#8899bb', fontSize: 9 },
        textStyle: { color: '#aabbcc', fontSize: 9 },
        gridlines: { color: '#2e4070' }
      },
      series: { 0: { color: '#6baed6', lineWidth: 2, pointSize: 4 } },
      legend: { position: 'none' },
      backgroundColor: '#1a2a40',
      chartArea: { backgroundColor: '#1a2a40', width: '78%', height: '65%' },
      width: 270, height: 190
    });
    chartPanel.add(chart);

    var validVals = vals.filter(function(v){ return v !== null; });
    if (validVals.length > 0) {
      var minV = Math.min.apply(null, validVals);
      var maxV = Math.max.apply(null, validVals);
      chartPanel.add(mkLabel(
        'Min: ' + minV.toFixed(2) + ' gC/m2 (' + allYears[vals.indexOf(minV)] + ')   ' +
        'Max: ' + maxV.toFixed(2) + ' gC/m2 (' + allYears[vals.indexOf(maxV)] + ')',
        { color: '#88aacc', backgroundColor: '#16213e',
          fontSize: '10px', fontStyle: 'italic', margin: '2px 0 0 0' }
      ));
    }
  });
}

controlPanel.add(divider());

controlPanel.add(sectionHead('Pixel Values (click map)'));
controlPanel.add(mkLabel(
  'Click any pixel to see all proxy values at that location.',
  { color: '#8899bb', backgroundColor: '#16213e',
    fontSize: '10px', whiteSpace: 'wrap', margin: '0 0 4px 0' }
));

var coordLabel = mkLabel('', {
  color: '#ffd580', backgroundColor: '#16213e',
  fontSize: '10px', whiteSpace: 'wrap', margin: '0 0 4px 0'
});
controlPanel.add(coordLabel);

var valueRows = [];
var swatchColors = [
  PAL_PA[7], PAL_MORPH[7], PAL_MORPH[7], PAL_MORPH[7],
  PAL_SPC[7], PAL_AGB[7], PAL_CINDEX[7], PAL_AGC[7]
];
for (var ri = 0; ri < LAYERS.length; ri++) {
  var swatch = mkLabel(' ', {
    backgroundColor: swatchColors[ri],
    width: '10px', height: '10px',
    margin: '3px 5px 0 0', border: '1px solid #ffffff'
  });
  var nameL = mkLabel(LAYERS[ri].short + ':', {
    color: '#cce0ff', backgroundColor: '#16213e',
    fontSize: '10px', width: '75px', margin: '2px 4px 0 0'
  });
  var valL = mkLabel('--', {
    color: '#ffffff', backgroundColor: '#16213e',
    fontSize: '10px', fontWeight: 'bold', stretch: 'horizontal'
  });
  var unitL = mkLabel(LAYERS[ri].unit, {
    color: '#8899bb', backgroundColor: '#16213e',
    fontSize: '10px', margin: '2px 0 0 4px'
  });
  controlPanel.add(ui.Panel({
    widgets: [swatch, nameL, valL, unitL],
    layout: ui.Panel.Layout.flow('horizontal'),
    style: { backgroundColor: '#16213e', margin: '1px 0', stretch: 'horizontal' }
  }));
  valueRows.push(valL);
}

controlPanel.add(divider());

controlPanel.add(mkLabel(
  'm.hafizt@uq.edu.au. A Proxy-Based Strategy for Estimating ' +
  'Seagrass Above-Ground Carbon from Satellite Observations. Remote Sensing.',
  { color: '#4a5a70', backgroundColor: '#16213e',
    fontSize: '10px', whiteSpace: 'wrap', margin: '0 0 4px 0' }
));

// ============================================================
// Sentinel-2 helper -- clipped to persistence mask
// ============================================================
function addS2Layer(year) {
  function maskS2(img) {
    var scl = img.select('SCL');
    return img.updateMask(
      scl.neq(1).and(scl.neq(3)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10))
    );
  }
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filterDate(ee.Date.fromYMD(year, 1, 1),
                ee.Date.fromYMD(ee.Number(year).add(1), 1, 1))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
    .map(maskS2)
    .select(['B4', 'B3', 'B2'])
    .reduce(ee.Reducer.percentile([60]))
    .rename(['B4', 'B3', 'B2'])
    .updateMask(persistenceMask)
    .clip(aoi);
  mainMap.addLayer(
    s2,
    { bands: ['B4', 'B3', 'B2'], min: 0, max: 3000, gamma: 1.4 },
    'S2 p60 . ' + year, true
  );
}

// ============================================================
// CORE -- updateLayer()
// ============================================================
function updateLayer() {
  var year = Math.round(yearSlider.getValue());
  var lbl  = layerSelect.getValue();

  yearDisp.setValue(String(year));
  mainMap.layers().reset();

  addS2Layer(year);   // S2 basemap for selected year

  var cfg     = LAYER_MAP[lbl];
  var clipped = cfg.getImage(year)
                  .updateMask(persistenceMask)
                  .clip(aoi);

  if (cfg.min !== null && cfg.max !== null) {
    mainMap.addLayer(
      clipped,
      { min: cfg.min, max: cfg.max, palette: cfg.palette },
      lbl + ' ' + year
    );
    updateLegend(cfg, cfg.min, cfg.max);
  } else {
    clipped.reduceRegion({
      reducer:    ee.Reducer.percentile([2, 98]),
      geometry:   aoi,
      scale:      100,
      maxPixels:  1e9,
      bestEffort: true
    }).evaluate(function(stats) {
      var band = cfg.bandName;
      var lo = (stats && stats[band + '_p2']  !== undefined) ? stats[band + '_p2']  : 0;
      var hi = (stats && stats[band + '_p98'] !== undefined) ? stats[band + '_p98'] : 1;
      mainMap.addLayer(
        clipped,
        { min: lo, max: hi, palette: cfg.palette },
        lbl + ' ' + year
      );
      updateLegend(cfg, lo, hi);
    });
  }
}

var MARKER_LAYER_NAME = 'Clicked Point';

mainMap.onClick(function(coords) {
  var year = Math.round(yearSlider.getValue());
  var pt   = ee.Geometry.Point([coords.lon, coords.lat]);

  coordLabel.setValue(
    'lat ' + coords.lat.toFixed(5) + '   lon ' + coords.lon.toFixed(5) +
    '   |   year: ' + year
  );

  for (var k = 0; k < valueRows.length; k++) {
    valueRows[k].setValue('...');
  }
  inspectorStatus.setValue('Querying all layers...');

  var layers = mainMap.layers();
  for (var li = layers.length() - 1; li >= 0; li--) {
    if (layers.get(li).getName() === MARKER_LAYER_NAME) {
      layers.remove(layers.get(li));
    }
  }
  mainMap.addLayer(
    ee.Image().paint(ee.FeatureCollection([ee.Feature(pt)]), 1, 3),
    { palette: ['#ffff00'] },
    MARKER_LAYER_NAME
  );

  var morphImg = ee.Image(MORPH_PREFIX + year);
  var stack = ee.Image([
    ee.Image(PA_PROB_PREFIX + year).rename('PA_prob'),
    morphImg.select('P_mono_Ea'),
    morphImg.select('P_mixed_long'),
    morphImg.select('P_mixed_short_plus_mono_short'),
    ee.Image(SPC_PREFIX    + year).rename('tSPC_pred'),
    ee.Image(AGB_PREFIX    + year).rename('tAGB_pred'),
    ee.Image(CINDEX_PREFIX + year).rename('carbon_index'),
    ee.Image(AGC_PREFIX    + year).rename('AGC_pred')
  ]);

  var bandNames = [
    'PA_prob', 'P_mono_Ea', 'P_mixed_long',
    'P_mixed_short_plus_mono_short',
    'tSPC_pred', 'tAGB_pred', 'carbon_index', 'AGC_pred'
  ];

  stack.reduceRegion({
    reducer:   ee.Reducer.first(),
    geometry:  pt,
    scale:     SCALE,
    maxPixels: 1e6
  }).evaluate(function(result) {
    if (!result) {
      inspectorStatus.setValue('No data at this location.');
      for (var k = 0; k < valueRows.length; k++) {
        valueRows[k].setValue('--');
      }
      return;
    }

    var allNull = true;
    for (var k = 0; k < bandNames.length; k++) {
      var v = result[bandNames[k]];
      if (v !== null && v !== undefined) {
        allNull = false;
        valueRows[k].setValue(v.toFixed(4));
      } else {
        valueRows[k].setValue('--');
      }
    }
    inspectorStatus.setValue(
      allNull
        ? 'Outside seagrass mask or AOI.'
        : 'Values for year ' + year
    );
  });

  buildAGCChart(pt, coords.lat, coords.lon);
});

(function() {
  var si = SITES['Ayau'];
  mainMap.setCenter(si[0], si[1], si[2]);
  updateLayer();
  updateErrorInfo();
  updateAGCSummary();
})();