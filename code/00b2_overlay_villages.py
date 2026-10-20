"""
00b2 · Oil palm on peat per village: overlay of the 00b shapefiles with the village polygons.

The article's values were produced with an equivalent overlay in a desktop GIS, with areas
calculated on the WGS 84 ellipsoid. This script reproduces that step with open-source tools, so
the result can be checked against final_merged_dataset.csv.

Inputs (place in this folder):
  Industrial_OilPalm_on_Peat_2021_SHP*.shp   from 00b_gee_oilpalm_on_peat_polygons.js
  Smallholder_OilPalm_on_Peat_2021_SHP*.shp  (Earth Engine may split large exports into several files)
  FinalGeoPortal.gpkg                        Ina-Geoportal village boundaries
  final_merged_dataset.csv                   optional: for the comparison
Output:
  oilpalm_on_peat_by_village.csv   KDEPUM, peat_ip_ha, peat_sp_ha, peat_tp_ha
"""
import glob
import geopandas as gpd
import pandas as pd
from pyproj import Geod
from shapely import force_2d

geod = Geod(ellps="WGS84")
ga = lambda geoms: [abs(geod.geometry_area_perimeter(g)[0]) / 1e4 for g in geoms]   # ha, ellipsoidal

v = gpd.read_file("FinalGeoPortal.gpkg").to_crs(4326)
v["geometry"] = force_2d(v.geometry)
v = v[["KDEPUM", "geometry"]]

out = v[["KDEPUM"]].copy()
for label, col in [("Industrial", "peat_ip_ha"), ("Smallholder", "peat_sp_ha")]:
    files = sorted(glob.glob(f"{label}_OilPalm_on_Peat_2021_SHP*.shp"))
    assert files, f"no {label} shapefile found"
    p = pd.concat([gpd.read_file(f).to_crs(4326) for f in files], ignore_index=True)
    p = gpd.GeoDataFrame(p[["geometry"]], crs=4326)
    print(f"{label}: {len(p):,} polygons from {len(files)} file(s)")
    x = gpd.overlay(v, p, how="intersection", keep_geom_type=True)
    x[col] = ga(x.geometry)
    out = out.merge(x.groupby("KDEPUM", as_index=False)[col].sum(), on="KDEPUM", how="left")
out = out.fillna(0)
out["peat_tp_ha"] = out.peat_ip_ha + out.peat_sp_ha
out.to_csv("oilpalm_on_peat_by_village.csv", index=False)
print(f"wrote oilpalm_on_peat_by_village.csv ({(out.peat_tp_ha > 0).sum()} villages with oil palm on peat)")

try:   # comparison with the values used in the article
    d = pd.read_csv("final_merged_dataset.csv", low_memory=False)
    d = d[d.peat_tp_ha.notna()].drop_duplicates("id_desa")
    d["id"] = d.id_desa.astype(str).str.replace(".", "", regex=False)
    out["id"] = out.KDEPUM.astype(str).str.replace(".", "", regex=False)
    m = d.merge(out, on="id", suffixes=("_used", "_rerun"))
    for c in ["peat_tp_ha", "peat_ip_ha", "peat_sp_ha"]:
        diff = (m[f"{c}_rerun"] - m[f"{c}_used"]).abs()
        print(f"{c}: r = {m[f'{c}_used'].corr(m[f'{c}_rerun']):.4f}; median |diff| = {diff.median():.2f} ha; "
              f"max |diff| = {diff.max():.1f} ha")
except FileNotFoundError:
    pass
