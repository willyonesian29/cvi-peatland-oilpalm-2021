"""
Composite Vulnerability Index (CVI) of villages with oil palm on peatland in
Indonesia, 2021 — baseline analysis pipeline.

Reproduces the baseline results reported in the manuscript:
  - Data cleaning / analytic sample definition
  - Exposure, Sensitivity, Adaptive Capacity construction
  - Multicollinearity check
  - Equal-weight CVI + PCA-weight robustness check
  - Jenks natural-breaks classification
  - Province-level Kruskal-Wallis test
  - Smallholder- vs industrial-dominant Mann-Whitney test
  - Spearman correlations (pillars vs CVI, smallholder share vs CVI)
  - Summary tables as CSV

Input : analysis_input.csv (from 00d; IDM indicators, oil palm on peat and
        Ina-Geoportal village area; built from data not redistributed — see README)
Output: village_vulnerability_index.csv, table1.csv, table2.csv, table3.csv,
        figure1.png, figure2.png (exploratory plots, not the manuscript figures)

Requirements: pandas, numpy, scikit-learn, scipy, mapclassify, matplotlib
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import minmax_scale
from scipy.stats import spearmanr, kruskal, mannwhitneyu
import mapclassify

pd.set_option("display.width", 140)

INPUT_PATH = "analysis_input.csv"   # from 00d_prepare_analysis_input.py

# =================================================================
# 1. LOAD DATA & DEFINE ANALYTIC SAMPLE
# =================================================================
df = pd.read_csv(INPUT_PATH)

# Analytic sample: villages where oil palm on peat is actually present
sample = df[df["peat_tp_ha"].notna()].copy()

# Data cleaning: drop duplicate village codes (keep first occurrence)
n_before = len(sample)
sample = sample.drop_duplicates(subset="id_desa", keep="first")
print(f"Analytic sample: {len(sample)} villages "
      f"(dropped {n_before - len(sample)} duplicate id_desa)")

# =================================================================
# 2. EXPOSURE
# =================================================================
sample["exposure_ratio"] = sample["peat_tp_ha"] / sample["VillageArea_Ha"]
n_capped = (sample["exposure_ratio"] > 1).sum()
sample["exposure_ratio"] = sample["exposure_ratio"].clip(upper=1)
print(f"Exposure ratio capped at 1.0 for {n_capped} village(s) "
      f"(peat-palm area exceeded reported village area)")

sample["Exposure"] = minmax_scale(sample["exposure_ratio"])

# =================================================================
# 3. SENSITIVITY
#    Indicators where a LOW raw score = worse/more sensitive -> invert
# =================================================================
sensitivity_cols = ["air_minum", "air_mck", "jamban", "tps",
                     "konflik", "cemar", "bencana", "ragam_eko"]
for c in sensitivity_cols:
    sample[c] = sample[c].fillna(sample[c].median())

sample["Sensitivity_raw"] = sample[sensitivity_cols].mean(axis=1)
sample["Sensitivity"] = minmax_scale(1 - sample["Sensitivity_raw"])

# =================================================================
# 4. ADAPTIVE CAPACITY
#    Higher raw score = more capacity -> used directly
# =================================================================
capacity_cols = [
    "faskes_km", "nakes", "posyandu", "posyandu_aktif", "bpjs", "sd", "smp",
    "sma", "kursus", "listrik_kk", "sinyal", "inet_desa", "inet_warga",
    "mitigasi", "pusat_dagang", "pasar", "minimarket", "logistik", "bank",
    "kredit", "koperasi", "angkot", "jalan_mobil", "jalan_mutu",
]
for c in capacity_cols:
    sample[c] = sample[c].fillna(sample[c].median())

# BUMDes / farmer group presence as binary capacity boosters
sample["bumdes_bin"] = (sample["bumdes_tani"] == "Ada").astype(int)
sample["l_tani_bin"] = (sample["l_tani"] == "Ada").astype(int)

capacity_all = capacity_cols + ["bumdes_bin", "l_tani_bin"]
sample["Capacity_raw"] = sample[capacity_all].mean(axis=1)
sample["Capacity"] = minmax_scale(sample["Capacity_raw"])

# =================================================================
# 5. MULTICOLLINEARITY CHECK (flag pairs with |r| > 0.8)
# =================================================================
print("\n--- Multicollinearity check ---")
for name, cols in [("Sensitivity", sensitivity_cols), ("Capacity", capacity_cols)]:
    corr = sample[cols].corr()
    high_pairs = [(i, j, corr.loc[i, j]) for i in corr.columns for j in corr.columns
                  if i < j and abs(corr.loc[i, j]) > 0.8]
    print(f"{name}: {len(high_pairs)} high-correlation pair(s)", high_pairs or "")

# =================================================================
# 6. COMPOSITE VULNERABILITY INDEX
# =================================================================
# Primary: equal-weight average of Exposure, Sensitivity, (1 - Capacity)
sample["CVI_equal"] = (sample["Exposure"] + sample["Sensitivity"]
                        + (1 - sample["Capacity"])) / 3

# Robustness check: PCA-derived weights
pillars = sample[["Exposure", "Sensitivity", "Capacity"]].copy()
pillars["Capacity"] = 1 - pillars["Capacity"]  # flip so all 3 point "more vulnerable"
pca = PCA(n_components=1)
pc1 = pca.fit_transform(pillars.values)
loadings = pca.components_[0]
print(f"\nPCA loadings (Exposure, Sensitivity, InvCapacity): {loadings}")
print(f"Variance explained by PC1: {pca.explained_variance_ratio_[0]:.3f}")
sample["CVI_pca"] = minmax_scale(pc1.flatten())

rho_weights, p_weights = spearmanr(sample["CVI_equal"], sample["CVI_pca"])
print(f"Spearman (equal-weight vs PCA-weight CVI): rho={rho_weights:.3f}, p={p_weights:.4g}")

# =================================================================
# 7. CLASSIFICATION — Jenks natural breaks (5 classes) on equal-weight CVI
# =================================================================
labels = ["Very Low", "Low", "Moderate", "High", "Very High"]
# NaturalBreaks uses a k-means-style random search internally (10 restarts by
# default) without an exposed random_state, so results can vary slightly
# between runs unless the global numpy seed is fixed beforehand.
np.random.seed(42)
jenks = mapclassify.NaturalBreaks(sample["CVI_equal"], k=5)
sample["CVI_class"] = pd.cut(sample["CVI_equal"], bins=[-0.001] + list(jenks.bins),
                              labels=labels)

print("\n--- Jenks breaks (equal-weight CVI) ---")
print(jenks.bins)
print("\n--- Class distribution ---")
print(sample["CVI_class"].value_counts().reindex(labels))

# =================================================================
# 8. STATISTICAL TESTS
# =================================================================
# Pillar correlations with CVI
for pillar in ["Exposure", "Sensitivity", "Capacity"]:
    r, pv = spearmanr(sample[pillar], sample["CVI_equal"])
    print(f"Spearman {pillar} vs CVI_equal: rho={r:.3f}, p={pv:.4g}")

# Smallholder share vs CVI
sample["smallholder_share"] = sample["peat_sp_ha"] / sample["peat_tp_ha"].replace(0, np.nan)
rho_sh, p_sh = spearmanr(sample["smallholder_share"].fillna(0), sample["CVI_equal"])
print(f"Spearman smallholder share vs CVI: rho={rho_sh:.3f}, p={p_sh:.4g}")

# Kruskal-Wallis across provinces with n >= 10
prov_counts = sample["prov"].value_counts()
provs_keep = prov_counts[prov_counts >= 10].index
sub = sample[sample["prov"].isin(provs_keep)]
groups = [g["CVI_equal"].values for _, g in sub.groupby("prov")]
h_stat, p_kw = kruskal(*groups)
print(f"\nKruskal-Wallis CVI by province (k={len(groups)} provinces, n>=10 each): "
      f"H={h_stat:.3f}, p={p_kw:.4g}")

# Mann-Whitney: smallholder-dominant vs industrial-dominant villages
sample["dominant_type"] = np.where(sample["peat_sp_ha"] > sample["peat_ip_ha"],
                                    "Smallholder-dominant", "Industrial-dominant")
g_small = sample.loc[sample["dominant_type"] == "Smallholder-dominant", "CVI_equal"]
g_ind = sample.loc[sample["dominant_type"] == "Industrial-dominant", "CVI_equal"]
u_stat, p_mw = mannwhitneyu(g_small, g_ind)
print(f"Mann-Whitney smallholder-dominant vs industrial-dominant CVI: "
      f"U={u_stat:.0f}, p={p_mw:.4g}")
print(f"  median CVI (smallholder-dominant) = {g_small.median():.3f}")
print(f"  median CVI (industrial-dominant)  = {g_ind.median():.3f}")

# =================================================================
# 9. SUMMARY TABLES
# =================================================================
# Table 1: descriptive statistics
table1 = sample[["Exposure", "Sensitivity", "Capacity", "CVI_equal"]] \
    .agg(["mean", "std", "median", "min", "max"]).T
table1.to_csv("table1.csv")
print("\nTable 1 (descriptive stats):\n", table1)

# Table 2: class distribution
table2 = sample["CVI_class"].value_counts().reindex(labels).rename("n").to_frame()
table2["%"] = (table2["n"] / len(sample) * 100).round(1)
table2.to_csv("table2.csv")
print("\nTable 2 (class distribution):\n", table2)

# Table 3: province-level means, sorted descending
table3 = sample.groupby("prov").agg(
    n=("CVI_equal", "size"),
    mean_cvi=("CVI_equal", "mean"),
    sd_cvi=("CVI_equal", "std"),
).sort_values("mean_cvi", ascending=False)
table3.to_csv("table3.csv")
print("\nTable 3 (province means):\n", table3)

# =================================================================
# 10. SAVE FULL VILLAGE-LEVEL OUTPUT
# =================================================================
out_cols = [
    "prov", "kab", "kec", "id_desa", "desa", "VillageArea_Ha",
    "peat_tp_ha", "peat_sp_ha", "peat_ip_ha",
    "exposure_ratio", "Exposure", "Sensitivity", "Capacity",
    "CVI_equal", "CVI_pca", "CVI_class", "dominant_type",
]
sample[out_cols].to_csv("village_vulnerability_index.csv", index=False)
print("\nSaved: village_vulnerability_index.csv")

# =================================================================
# 11. EXPLORATORY FIGURES (not the manuscript figures; see 06 and 07)
# =================================================================
colors = {"Very Low": "#2ecc71", "Low": "#a3d977", "Moderate": "#f1c40f",
          "High": "#e67e22", "Very High": "#c0392b"}

fig, ax = plt.subplots(figsize=(7, 4.3))
counts = sample["CVI_class"].value_counts().reindex(labels)
bars = ax.bar(counts.index, counts.values, color=[colors[c] for c in counts.index])
for b, v in zip(bars, counts.values):
    ax.text(b.get_x() + b.get_width() / 2, v + 3, str(v), ha="center", fontsize=9)
ax.set_ylabel("Number of villages")
ax.set_title("Distribution of Villages by CVI Class (n = %d)" % len(sample))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("figure1.png", dpi=200, bbox_inches="tight")
plt.close()

fig, ax = plt.subplots(figsize=(7, 4.5))
for cls, sub_g in sample.groupby("CVI_class", observed=True):
    ax.scatter(sub_g["Exposure"], 1 - sub_g["Capacity"], s=14, alpha=0.65,
               label=cls, color=colors.get(cls, "grey"))
ax.set_xlabel("Exposure (normalized)")
ax.set_ylabel("1 - Adaptive Capacity")
ax.set_title("Exposure vs. Capacity Deficit by CVI Class")
ax.legend(fontsize=8, loc="upper right", title="CVI Class")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("figure2.png", dpi=200, bbox_inches="tight")
plt.close()

print("\nSaved: figure1.png, figure2.png")
print("\nDone.")
