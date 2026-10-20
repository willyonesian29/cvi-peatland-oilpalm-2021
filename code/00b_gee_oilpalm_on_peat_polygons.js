// =========================================================================
// 00b · Oil palm on peatland, 2021 (Google Earth Engine Code Editor)
// =========================================================================
// Original script used for the article (cleaned; logic unchanged).
//
// What it does
//   1. Rasterizes the BIG peatland polygons (asset PGaFin, 1:50,000) into a 10 m peat mask.
//   2. Keeps the 2021 oil palm pixels of Descals et al. (2024) that fall on peat,
//      separately for industrial (class 1) and smallholder (class 2) plantations.
//   3. Converts them to polygons (reduceToVectors, 10 m, native resolution) and exports
//      two shapefiles: industrial and smallholder oil palm on peat.
//   No village join is done here: the shapefiles are overlaid with the Ina-Geoportal
//   village polygons in the next step to obtain peat_ip_ha, peat_sp_ha and peat_tp_ha
//   per village (see README, "Spatial processing").
//
// Assets (set to your own paths)
//   PGaFin : BIG national peatland map (1:50,000), polygons
//   Oil palm 2021 extent: awesome-gee-community catalog (public)
// =========================================================================

var peatVector = ee.FeatureCollection('projects/willyonesian29/assets/PGaFin');

// Peat mask: 1 inside peat, masked elsewhere
var peatMaskRaw = ee.Image().byte().paint(peatVector, 1).selfMask();

// -------------------------------------------------------------------------
// STEP A (run once): store the peat mask as an asset. Painting a national
// vector on the fly and chaining more operations on top of it makes the
// computation graph too large ("Description length exceeds maximum",
// reduceToVectors time-outs); exporting it once flattens the graph.
// -------------------------------------------------------------------------
Export.image.toAsset({
  image: peatMaskRaw,
  description: 'peatMask_from_PGaFin',
  assetId: 'projects/willyonesian29/assets/peatMask_baked',
  region: peatVector.geometry(),
  scale: 10,
  maxPixels: 1e13
});

// -------------------------------------------------------------------------
// STEP B: once the Step A task has finished, comment out Step A and use the
// stored mask instead of peatMaskRaw.
// -------------------------------------------------------------------------
// var peatMask = ee.Image('projects/willyonesian29/assets/peatMask_baked');
var peatMask = peatMaskRaw;

// Oil palm 2021 (1 = industrial, 2 = smallholder), restricted to the peat extent
var peatExtent = peatVector.geometry();
var oilPalm2021 = ee.ImageCollection('projects/sat-io/open-datasets/global-oil-palm/GlobalOilPalm_extent_2021')
  .mosaic()
  .clip(peatExtent);

var industrialOilPalmPeat  = oilPalm2021.eq(1).and(peatMask).selfMask();
var smallholderOilPalmPeat = oilPalm2021.eq(2).and(peatMask).selfMask();

// Raster -> polygons at 10 m (keeps polygon area faithful to the source map)
function toPolygons(img, label) {
  return img.reduceToVectors({
    geometry: peatExtent,
    scale: 10,
    geometryType: 'polygon',
    eightConnected: false,
    labelProperty: 'class',
    maxPixels: 1e13,
    tileScale: 16
  }).map(function (f) { return f.set('opclass', label); });
}
var industrialVectors  = toPolygons(industrialOilPalmPeat, 'industrial');
var smallholderVectors = toPolygons(smallholderOilPalmPeat, 'smallholder');

// Do not Map.addLayer() the paint()/and()/selfMask() chains (graph too large)
Map.setCenter(117.3, -0.5, 5);

// Export two shapefiles (no village join, no reducers)
Export.table.toDrive({
  collection: industrialVectors,
  description: 'Industrial_OilPalm_on_Peat_2021_SHP',
  fileFormat: 'SHP'
});
Export.table.toDrive({
  collection: smallholderVectors,
  description: 'Smallholder_OilPalm_on_Peat_2021_SHP',
  fileFormat: 'SHP'
});
