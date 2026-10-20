"""
Consistency checks on the analysis input (run after 00d, before 01).

Checks that the plantation-type components of oil palm on peat add up, and lists
villages whose oil palm on peat exceeds the Ina-Geoportal village area (exposure
ratio > 1; capped at 1 in the analysis).

Input : analysis_input.csv (from 00d_prepare_analysis_input.py)
Output: zonal_consistency_flags.csv
"""
import pandas as pd

s = pd.read_csv("analysis_input.csv", low_memory=False).drop_duplicates("id_desa")
print(f"Study villages: {len(s)}")
d = (s.peat_tp_ha - s.peat_ip_ha.fillna(0) - s.peat_sp_ha.fillna(0)).abs().max()
print(f"max |peat_tp - peat_ip - peat_sp| = {d:.2e} ha")
s["exposure_ratio_raw"] = s.peat_tp_ha / s.VillageArea_Ha
f = s[s.exposure_ratio_raw > 1][["prov", "kab", "desa", "id_desa", "VillageArea_Ha", "peat_tp_ha",
                                 "exposure_ratio_raw"]].sort_values("exposure_ratio_raw", ascending=False)
print(f"exposure ratio > 1 (capped at 1): {len(f)} village(s)")
print(f.round(2).to_string(index=False))
f.round(3).to_csv("zonal_consistency_flags.csv", index=False)
print("Wrote zonal_consistency_flags.csv")
