# Composite Vulnerability Index of Villages with Oil Palm on Peatland in Indonesia, 2021

Data and code for the article:

> Siagian, W. O., & Putri, K. (2026). Composite Vulnerability Index of Villages with Oil Palm on Peatland in Indonesia in 2021: Integrating Remote Sensing and Village Development Index Data. *Jurnal Ilmiah Membangun Desa dan Pertanian (JIMDP)*. [DOI to be added]

The study builds a Composite Vulnerability Index (CVI) for 901 Indonesian villages whose territory contained oil palm on peatland in 2021. It uses an Exposure–Sensitivity–Adaptive Capacity framework. Exposure is the share of village area under oil palm on peat, from remote sensing over Ina-Geoportal village boundaries. Sensitivity and Adaptive Capacity come from the 2021 Village Development Index (Indeks Desa Membangun, IDM). The repository contains the village-level results and all code used to produce the tables and figures.

## Repository contents

```
cvi-peatland-oilpalm-2021/
├── README.md
├── LICENSE                     MIT licence (code)
├── CITATION.cff                citation metadata
├── requirements.txt            Python packages
├── data/
│   ├── S2_village_cvi_scores.csv   Supplementary Material S2 (901 villages)
│   └── LICENSE                 CC BY 4.0 (data)
└── code/
    ├── 00b_gee_oilpalm_on_peat_polygons.js    Earth Engine: oil palm on peat as polygons (original script)
    ├── 00b2_overlay_villages.py               overlay with village polygons -> oil palm on peat per village
    ├── 00d_prepare_analysis_input.py          builds analysis_input.csv (Ina-Geoportal village area)
    ├── 00c_check_zonal_consistency.py         consistency checks on the analysis input
    ├── 01_cvi_baseline_analysis.py            baseline CVI, Jenks classes, Tables 2, 5, 6
    ├── 02_robustness_fixed_effects_spatial.py Table 3, fixed effects, Moran's I, LISA, spatial error model
    ├── 03_priority_group_and_bias_checks.py   Table 4, Tables S8–S9, exports S2
    ├── 04_supplementary_tables_S5_S6.py       Tables S5 and S6
    ├── 05_prepare_map_layers.py               joins S2 to the Ina-Geoportal polygons
    ├── 06_make_maps_figure2_S1.py             Figure 2 and Figure S1
    ├── 07_make_figure1_workflow.py            Figure 1
    └── extra/
        └── national_oilpalm_area_by_village_BPS.ipynb   total oil palm area for all villages (not used in the CVI)
```

## Data dictionary: `data/S2_village_cvi_scores.csv`

| Column | Description |
| --- | --- |
| `prov`, `kab`, `kec`, `desa` | Province, regency, district and village name |
| `id_desa` | Kemendagri village code (KDEPUM without dots) |
| `VillageArea_Ha` | Village area (ha): geodesic area of the Ina-Geoportal village polygon |
| `peat_tp_ha`, `peat_ip_ha`, `peat_sp_ha` | Oil palm area on peatland: total, industrial, smallholder (ha), summed over the same polygon |
| `exp_ratio` | Exposure ratio = `peat_tp_ha` / `VillageArea_Ha`, capped at 1 |
| `E`, `S`, `AC` | Exposure, Sensitivity and Adaptive Capacity pillars, min-max scaled 0–1 |
| `CVI` | Composite Vulnerability Index = (E + S + (1 − AC)) / 3 |
| `CVI_class` | Jenks natural-breaks class: Very Low ≤ 0.186, Low ≤ 0.277, Moderate ≤ 0.363, High ≤ 0.467, Very High ≤ 0.713 |
| `dom` | Dominant plantation type of oil palm on peat: `S` smallholder, `I` industrial (larger area) |
| `top` | `True` for the 207 High or Very High villages (the priority group) |
| `consensus` | `True` if in the priority group under all eight alternative constructions and in ≥ 90% of 1,000 Monte Carlo weight draws (54 villages) |
| `mc_top_freq` | Share of Monte Carlo draws in which the village is among the 207 highest-ranked |
| `priority_type` | For priority villages: `exposure-led` (exposure ≥ 75th percentile), `development-deficit` (exposure < median), `intermediate`; otherwise `not priority` |
| `LISA` | Local Moran's I cluster (p < 0.05): High-High, Low-Low, High-Low, Low-High or Not significant |

## Source data (not redistributed)

The input datasets belong to their providers and are not included here. To reproduce the analysis, obtain them from:

| Dataset | Source | Used for |
| --- | --- | --- |
| Global oil palm extent and planting year, 2021 layer (10 m, Sentinel-1) | Descals et al. (2024), Zenodo: https://doi.org/10.5281/zenodo.11034131 | Oil palm area and plantation type |
| National peatland map (1:50,000) | Geospatial Information Agency (BIG), One Map Indonesia | Peatland extent |
| Village boundaries (`FinalGeoPortal.gpkg`: 901 study villages) | Ina-Geoportal, BIG village boundaries, 1:10,000 base-map series, edition 28 September 2023 (`rbi10k_administrasi_desa_20230928`) | Zonal statistics, village area, spatial weights, maps |
| Village Development Index (IDM) 2021, indicator level | Ministry of Villages, Disadvantaged Regions Development, and Transmigration (Kemendes PDTT) | Sensitivity and Adaptive Capacity indicators |
| Admin 0 countries and Admin 1 provinces, 1:10m | Natural Earth (public domain): https://www.naturalearthdata.com | Map background only |

## Spatial processing

1. **Oil palm on peat (steps 00b, 00b2).** `00b_gee_oilpalm_on_peat_polygons.js` (the original Code Editor script) rasterizes the peat polygons to a 10 m mask. It selects the industrial and smallholder oil palm pixels on peat, converts them to polygons at 10 m and exports two shapefiles. These were overlaid with the Ina-Geoportal village polygons in a desktop GIS, and their areas, calculated on the WGS 84 ellipsoid, were summed per village. `00b2_overlay_villages.py` reproduces this overlay with open-source tools and compares the result with the values used.
2. **Join.** The results are joined to the IDM 2021 records by village code. `KDEPUM` (with dots) equals the IDM `id_desa` without dots. This gives `final_merged_dataset.csv`: 923 villages with oil palm on peat, 22 of them duplicate codes.
3. **Village area (00d).** Village area is computed from the same Ina-Geoportal polygons, so the numerator and denominator of Exposure share one boundary set. The result is written to `analysis_input.csv`.
4. **Checks (00c).** Three villages have oil palm on peat larger than their polygon area, a geometry artefact of the overlay: Sumber Makmur and Makarti Jaya (Nagan Raya) and Pawalutan (Hulu Sungai Utara), with ratios of 1.56–2.18. They are capped at 1, and excluding them changes nothing material (rho = 0.999; manuscript Table 3).

`final_merged_dataset.csv` contains the raw IDM data, so it is available from the corresponding author with permission from Kemendes PDTT.

## How to reproduce

1. Install Python 3.10 or later, then the packages:
   ```
   pip install -r requirements.txt
   ```
2. Place `final_merged_dataset.csv` and `FinalGeoPortal.gpkg` in `code/` and run, from inside `code/`:
   ```
   python 00d_prepare_analysis_input.py
   python 00c_check_zonal_consistency.py
   python 01_cvi_baseline_analysis.py
   python 02_robustness_fixed_effects_spatial.py
   python 03_priority_group_and_bias_checks.py
   python 04_supplementary_tables_S5_S6.py
   ```
3. For the maps, also place the Natural Earth files `ne_10m_admin_0_countries.shp` and `ne_10m_admin_1_states_provinces.shp` in `code/`, then run:
   ```
   python 05_prepare_map_layers.py
   python 06_make_maps_figure2_S1.py
   python 07_make_figure1_workflow.py
   ```

Every script prints the numbers it produces and writes CSV tables. Jenks classification uses a fixed random seed (42), so the class counts reproduce exactly: 154 / 278 / 262 / 142 / 65.

## Licence

- Code: MIT licence (`LICENSE`).
- Data in `data/`: Creative Commons Attribution 4.0 International (`data/LICENSE`).

## Contact

Willy Onesimus Siagian, Urban and Regional Planning Master Programme Diponegoro University. [willy.onesian@gmail.com]
