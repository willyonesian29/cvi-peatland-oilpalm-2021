"""
Prepare the village layer for Figure 2 and Figure S1.

Inputs (not redistributed — see README):
  FinalGeoPortal.gpkg                 Ina-Geoportal (BIG) village boundaries of the 901 study villages
  ../data/S2_village_cvi_scores.csv   village-level CVI scores, classes and LISA classes (this repository)
Output:
  villages.gpkg                       901 polygons (2D, WGS 84) with CVI class, LISA class and regency

Villages are joined by Kemendagri village code (KDEPUM without dots = id_desa); all 901 match.
"""
import geopandas as gpd
import pandas as pd
from shapely import force_2d

g = gpd.read_file("FinalGeoPortal.gpkg").to_crs(4326)
g["geometry"] = force_2d(g.geometry)
g["id_str"] = g["KDEPUM"].astype(str).str.replace(".", "", regex=False)
s = pd.read_csv("../data/S2_village_cvi_scores.csv", dtype={"id_desa": str})
s["id_str"] = s["id_desa"].str.replace(".", "", regex=False)
v = g[["id_str", "geometry"]].merge(s, on="id_str", how="inner")
print(f"matched {len(v)} of {len(s)} study villages")
assert len(v) == len(s)
v[["id_desa", "prov", "kab", "desa", "CVI", "CVI_class", "LISA", "consensus", "dom", "geometry"]].to_file("villages.gpkg")
print("wrote villages.gpkg")
