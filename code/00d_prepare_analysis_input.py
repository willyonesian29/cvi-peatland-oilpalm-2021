"""
Build the analysis input from the merged village table and the Ina-Geoportal boundaries.

Exposure must use the same boundaries in its numerator and denominator. Oil palm on
peat (peat_tp_ha, peat_ip_ha, peat_sp_ha) was summed over the Ina-Geoportal village
polygons, so village area is taken from the same polygons here (geodesic area on the
WGS 84 ellipsoid; identical to the layer's LUASWH field). The BPS-based columns in
final_merged_dataset.csv (VillageArea_Ha, Total/Ind/Small_Palm_Ha) are not used.

Inputs (not redistributed — see README):
  final_merged_dataset.csv   IDM 2021 indicators + oil palm on peat, joined by village code
  FinalGeoPortal.gpkg        Ina-Geoportal (BIG) village boundaries for the 901 study villages
Output:
  analysis_input.csv         input to scripts 00c and 01-04
"""
import geopandas as gpd
import pandas as pd
from pyproj import Geod

df = pd.read_csv("final_merged_dataset.csv", low_memory=False)
s = df[df["peat_tp_ha"].notna()].copy()                      # 923 rows (22 duplicate codes)
s = s.drop(columns=["VillageArea_Ha", "Total_Palm_Ha", "Ind_Palm_Ha", "Small_Palm_Ha"])
s["id_str"] = s["id_desa"].astype(str).str.replace(".", "", regex=False)

g = gpd.read_file("FinalGeoPortal.gpkg").to_crs(4326)
g["id_str"] = g["KDEPUM"].astype(str).str.replace(".", "", regex=False)
geod = Geod(ellps="WGS84")
g["VillageArea_Ha"] = [abs(geod.geometry_area_perimeter(geom)[0]) / 1e4 for geom in g.geometry]
print(f"boundary polygons: {len(g)}; duplicate codes: {g.id_str.duplicated().sum()}")

s = s.merge(g[["id_str", "VillageArea_Ha"]], on="id_str", how="left").drop(columns="id_str")
miss = s["VillageArea_Ha"].isna().sum()
print(f"rows: {len(s)}; without boundary area: {miss}")
assert miss == 0, "every study village must have an Ina-Geoportal polygon"
s.to_csv("analysis_input.csv", index=False)
print("wrote analysis_input.csv")
