"""
CVI second-round revision analyses (JIMDP #3227, reviewer comments of October 2026)
===================================================================================
Run in the same folder as final_merged_dataset.csv (the file used by
cvi_full_analysis.py and cvi_revision_robustness.py). It rebuilds the baseline
CVI exactly as before and prints every new number cited in the revised
manuscript. Each number carries a tag [V1], [V3], ... that matches a
yellow-highlighted placeholder in Manuscript_JIMDP_3227_revision2_marked.docx.
All tags are also written to revision2_numbers.txt.

Outputs
  revision2_numbers.txt        -> all [V..] values for the manuscript
  table4_priority_group.csv    -> new Table 4 (composition/stability of priority group)
  tableS8_bias_checks.csv      -> Table S8 (smallholder-industrial bias checks)
  tableS9_class_transitions.csv-> Table S9 (baseline x variant class transition matrices, kappa)
  S2_village_cvi_scores.csv    -> Supplementary Material S2 (village-level scores and classes;
                                  adds the LISA class if lisa_clusters.csv from script 02 is present)

Requirements
  pip install pandas numpy scipy scikit-learn mapclassify statsmodels
"""

import numpy as np
import pandas as pd
import mapclassify
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.metrics import cohen_kappa_score
from sklearn.preprocessing import minmax_scale

INPUT_PATH = "analysis_input.csv"   # from 00d_prepare_analysis_input.py
N_MC = 1000
SEED = 42
LABELS = ["Very Low", "Low", "Moderate", "High", "Very High"]
OUT = {}


def tag(key, value, label):
    OUT[key] = (value, label)
    print(f"[{key}] {label}: {value}")


# ------------------------------------------------------------------ baseline (identical to earlier scripts)
SENS = ["air_minum", "air_mck", "jamban", "tps", "konflik", "cemar", "bencana", "ragam_eko"]
CAP = ["faskes_km", "nakes", "posyandu", "posyandu_aktif", "bpjs", "sd", "smp", "sma", "kursus",
       "listrik_kk", "sinyal", "inet_desa", "inet_warga", "mitigasi", "pusat_dagang", "pasar",
       "minimarket", "logistik", "bank", "kredit", "koperasi", "angkot", "jalan_mobil", "jalan_mutu"]
BIN = ["bumdes_bin", "l_tani_bin"]

df = pd.read_csv(INPUT_PATH)
s = df[df["peat_tp_ha"].notna()].drop_duplicates("id_desa", keep="first").reset_index(drop=True)
for c in SENS + CAP:
    s[c] = s[c].fillna(s[c].median())
s["bumdes_bin"] = (s["bumdes_tani"] == "Ada").astype(int)
s["l_tani_bin"] = (s["l_tani"] == "Ada").astype(int)
s["exp_raw"] = s["peat_tp_ha"] / s["VillageArea_Ha"]
capped_idx = s.index[s["exp_raw"] > 1]
s["exp_ratio"] = s["exp_raw"].clip(upper=1)
s["dom"] = np.where(s["peat_sp_ha"] > s["peat_ip_ha"], "S", "I")   # dominant type of oil palm on peat
s["industrial"] = (s["dom"] == "I").astype(int)


def pillars(d, sens=SENS, cap=CAP + BIN, exposure="minmax"):
    if exposure == "minmax":
        E = minmax_scale(d["exp_ratio"])
    elif exposure == "rank":
        E = minmax_scale(d["exp_ratio"].rank(pct=True))
    else:  # winsor99
        E = minmax_scale(d["exp_ratio"].clip(upper=d["exp_ratio"].quantile(.99)))
    S = minmax_scale(1 - d[sens].mean(axis=1))
    AC = minmax_scale(d[cap].mean(axis=1))
    return np.asarray(E), np.asarray(S), np.asarray(AC)


def cvi(E, S, AC, w=(1/3, 1/3, 1/3), agg="arith"):
    V = 1 - AC
    if agg == "arith":
        return w[0]*E + w[1]*S + w[2]*V
    eps = 0.01
    return np.exp(w[0]*np.log(E+eps) + w[1]*np.log(S+eps) + w[2]*np.log(V+eps))


def jenks_classes(v):
    np.random.seed(SEED)
    jb = mapclassify.NaturalBreaks(v, k=5)
    cls = np.searchsorted(jb.bins, v, side="left")       # 0..4, upper bound inclusive
    return np.clip(cls, 0, 4), jb.bins


E0, S0, AC0 = pillars(s)
s["E"], s["S"], s["AC"] = E0, S0, AC0
s["CVI"] = cvi(E0, S0, AC0)
s["cls"], bins = jenks_classes(s["CVI"].values)
s["CVI_class"] = [LABELS[i] for i in s["cls"]]
rank0 = s["CVI"].rank(ascending=False, method="first")
TOPN = int((s["cls"] >= 3).sum())                 # High + Very High
s["top"] = rank0 <= TOPN
assert s["top"].sum() == TOPN
print(f"Priority group (High + Very High): {TOPN} villages")

# V1 Jenks bounds
tag("V1", "; ".join(f"{LABELS[i]} ≤ {b:.3f}" for i, b in enumerate(bins)), "Jenks class upper bounds")
for k, b in zip(["V1a", "V1b", "V1c", "V1d"], bins[:4]):
    tag(k, f"{b:.3f}", "upper bound of class")
print("Check class counts:", s["CVI_class"].value_counts().reindex(LABELS).tolist(),
      "(manuscript: 156, 260, 260, 151, 74)")

# V8 villages within 0.01 of the lower bound of High (= upper bound of Moderate)
tag("V8", int((np.abs(s["CVI"] - bins[2]) <= 0.01).sum()), "villages within ±0.01 of the High lower bound")

# ------------------------------------------------------------------ variants
pcs = np.column_stack([E0, S0, 1 - AC0])
pc1 = PCA(1).fit_transform(pcs).ravel()
if np.corrcoef(pc1, s["CVI"])[0, 1] < 0:
    pc1 = -pc1
Er, Sr, ACr = pillars(s, exposure="rank")
Ew, Sw, ACw = pillars(s, exposure="winsor99")
Eb, Sb, ACb = pillars(s, cap=CAP)
keep_idx = s.index.difference(capped_idx)
Ed, Sd, ACd = pillars(s.loc[keep_idx])

variants = {  # name -> (scores, index)
    "PCA-derived weights": (minmax_scale(pc1), s.index),
    "Exposure-weighted (0.50/0.25/0.25)": (cvi(E0, S0, AC0, (.5, .25, .25)), s.index),
    "Socioeconomic-weighted (0.20/0.40/0.40)": (cvi(E0, S0, AC0, (.2, .4, .4)), s.index),
    "Geometric aggregation": (cvi(E0, S0, AC0, agg="geo"), s.index),
    "Rank-normalized Exposure": (cvi(Er, Sr, ACr), s.index),
    "Exposure winsorized (99th pct)": (cvi(Ew, Sw, ACw), s.index),
    "Capped village excluded": (cvi(Ed, Sd, ACd), keep_idx),
    "Without BUMDes and farmer-group items": (cvi(Eb, Sb, ACb), s.index),
}

in_top_all = pd.Series(True, index=s.index)
trans_rows, kappas = [], {}
for name, (v, idx) in variants.items():
    v = pd.Series(np.asarray(v), index=idx)
    r = v.rank(ascending=False, method="first")
    top_v = r <= TOPN
    in_top_all.loc[idx] &= top_v
    cls_v, _ = jenks_classes(v.values)
    base_c = s.loc[idx, "cls"].values
    kap = cohen_kappa_score(base_c, cls_v, weights="quadratic")
    changed = (base_c != cls_v).mean() * 100
    multi = (np.abs(base_c - cls_v) > 1).mean() * 100
    kappas[name] = (kap, changed, multi)
    tm = pd.crosstab(pd.Categorical([LABELS[i] for i in base_c], LABELS),
                     pd.Categorical([LABELS[i] for i in cls_v], LABELS))
    for bl in LABELS:
        trans_rows.append({"variant": name, "baseline_class": bl,
                           **{f"variant_{c}": int(tm.loc[bl, c]) for c in LABELS},
                           "weighted_kappa": round(kap, 3), "pct_changed_class": round(changed, 1),
                           "pct_changed_more_than_one_class": round(multi, 1)})
    if name == "Geometric aggregation":
        geo_top = top_v

pd.DataFrame(trans_rows).to_csv("tableS9_class_transitions.csv", index=False)
kv = [k[0] for k in kappas.values()]
tag("V6a", f"{min(kv):.2f} and {max(kv):.2f}", "range of quadratic-weighted kappa across variants")
g = kappas["Geometric aggregation"]
tag("V6b", f"{g[1]:.1f}", f"% villages changing class under geometric aggregation "
                          f"(of which more than one class: {g[2]:.1f}%)")

# Monte Carlo frequency of being in the top 225
rng = np.random.default_rng(SEED)
freq = np.zeros(len(s))
for _ in range(N_MC):
    wv = rng.dirichlet([1, 1, 1])
    r = pd.Series(cvi(E0, S0, AC0, tuple(wv))).rank(ascending=False, method="first")
    freq += (r <= TOPN).values
s["mc_top_freq"] = freq / N_MC

s["consensus"] = s["top"] & in_top_all & (s["mc_top_freq"] >= 0.90)
n_cons = int(s["consensus"].sum())
tag("V3", n_cons, "consensus priority villages")
tag("V3p", f"{n_cons / TOPN * 100:.1f}", "consensus as % of the priority group")
tag("V4", TOPN - n_cons, "construction-dependent priority villages")

# leavers / entrants under geometric aggregation
leave = s["top"] & ~geo_top.reindex(s.index, fill_value=False)
enter = ~s["top"] & geo_top.reindex(s.index, fill_value=False)
tag("V14", int(leave.sum()), "villages leaving (= entering) the top group under geometric aggregation")
tag("V5a", f"{s.loc[leave, 'exp_ratio'].median()*100:.1f}%", "median exposure ratio, leavers")
tag("V5b", f"{s.loc[enter, 'exp_ratio'].median()*100:.1f}%", "median exposure ratio, entrants")

# exposure-led vs development-deficit
p75, p50 = s["exp_ratio"].quantile(.75), s["exp_ratio"].median()
exp_led = s["top"] & (s["exp_ratio"] >= p75)
dev_def = s["top"] & (s["exp_ratio"] < p50)
tag("V7a", int(exp_led.sum()), "exposure-led priority villages")
tag("V7b", int(dev_def.sum()), "development-deficit priority villages")
tag("V7c", f"{p75*100:.1f}", "sample 75th percentile of exposure ratio (% of village area)")


def profile(mask, keys):
    d = s.loc[mask]
    vals = [f"{d['exp_ratio'].median()*100:.1f}", f"{d['S'].median():.3f}",
            f"{(1-d['AC']).median():.3f}", f"{d['industrial'].mean()*100:.1f}"]
    for k, v in zip(keys, vals):
        if k:
            tag(k, v, "Table 4 cell")
    return [int(mask.sum())] + vals


groups = [
    ("All High/Very High villages", s["top"], ["V11a", "V11b", "V11c", "V11d"]),
    ("Consensus priority (all variants)", s["consensus"], ["V12a", "V12b", "V12c", "V12d"]),
    ("Construction-dependent", s["top"] & ~s["consensus"], ["V13a", "V13b", "V13c", "V13d"]),
    ("  of which leave top group under geometric aggregation", leave, [None, "V14b", "V14c", "V14d"]),
    ("Villages entering top group under geometric aggregation", enter, [None, "V15b", "V15c", "V15d"]),
    ("Exposure-led (exposure >= 75th percentile)", exp_led, ["V16a", "V16b", "V16c", "V16d"]),
    ("Development-deficit (exposure < median)", dev_def, ["V17a", "V17b", "V17c", "V17d"]),
]
t4 = [[name] + profile(m, keys) for name, m, keys in groups]
pd.DataFrame(t4, columns=["group", "n", "median_exposure_ratio_pct", "median_S", "median_1_minus_AC",
                          "industrial_dominated_pct"]).to_csv("table4_priority_group.csv", index=False)

# ------------------------------------------------------------------ smallholder-industrial bias checks
s["sh_share"] = s["peat_sp_ha"] / s["peat_tp_ha"].replace(0, np.nan)
cl = lambda d: {"groups": pd.factorize(d["kab"])[0]}
rows = []

mixed = s["sh_share"].between(0.4, 0.6)
tag("V10a", int(mixed.sum()), "mixed villages (smallholder share 40-60%)")
d = s.loc[~mixed]
m1 = smf.ols("CVI ~ industrial + C(prov)", data=d).fit(cov_type="cluster", cov_kwds=cl(d))
tag("V10b", f"{m1.params['industrial']:.3f}, SE = {m1.bse['industrial']:.3f}, p = {m1.pvalues['industrial']:.2g}",
    "industrial beta, mixed villages excluded, province FE")
rows.append(["Mixed villages (40-60%) excluded", len(d), m1.params["industrial"], m1.bse["industrial"],
             m1.pvalues["industrial"]])

d = s.loc[(s["sh_share"] >= 2/3) | (s["sh_share"] <= 1/3)].copy()
d["dom23"] = np.where(d["sh_share"] >= 2/3, "S", "I")
u, p = mannwhitneyu(d.loc[d.dom23 == "S", "CVI"], d.loc[d.dom23 == "I", "CVI"])
tag("V10c", f"{d.loc[d.dom23=='S','CVI'].median():.3f} vs {d.loc[d.dom23=='I','CVI'].median():.3f}, "
            f"p = {p:.2g}, n = {len(d)}", "two-thirds dominance threshold, smallholder vs industrial")
rows.append(["Two-thirds dominance threshold (Mann-Whitney)", len(d), np.nan, np.nan, p])

s["log_area"] = np.log(s["VillageArea_Ha"])
s["log_palm"] = np.log1p(s["peat_tp_ha"])
m2 = smf.ols("CVI ~ industrial + log_area + log_palm + C(prov)", data=s).fit(cov_type="cluster", cov_kwds=cl(s))
tag("V10d", f"{m2.params['industrial']:.3f}, SE = {m2.bse['industrial']:.3f}, p = {m2.pvalues['industrial']:.2g}",
    "industrial beta with log village area + log total oil palm area, province FE")
rows.append(["Controls: log village area, log total oil palm area", len(s), m2.params["industrial"],
             m2.bse["industrial"], m2.pvalues["industrial"]])
for pil in ["S", "AC"]:
    m3 = smf.ols(f"{pil} ~ industrial + log_area + log_palm + C(prov)", data=s).fit(
        cov_type="cluster", cov_kwds=cl(s))
    rows.append([f"Pillar {pil} with same controls", len(s), m3.params["industrial"], m3.bse["industrial"],
                 m3.pvalues["industrial"]])
pd.DataFrame(rows, columns=["check", "n", "beta_industrial", "SE", "p"]).to_csv("tableS8_bias_checks.csv",
                                                                               index=False)

# ------------------------------------------------------------------ S2 village-level file
s["priority_type"] = np.select([exp_led, dev_def, s["top"]],
                               ["exposure-led", "development-deficit", "intermediate"], "not priority")
s2 = s[["prov", "kab", "kec", "id_desa", "desa", "VillageArea_Ha", "peat_tp_ha", "peat_ip_ha",
         "peat_sp_ha", "exp_ratio", "E", "S", "AC", "CVI", "CVI_class", "dom", "top",
         "consensus", "mc_top_freq", "priority_type"]].copy()
import os
if os.path.exists("lisa_clusters.csv"):          # written by 02_robustness_fixed_effects_spatial.py
    lab = {"HH": "High-High", "LL": "Low-Low", "HL": "High-Low", "LH": "Low-High", "ns": "Not significant"}
    l = pd.read_csv("lisa_clusters.csv", dtype={"id_str": str})
    s2["id_str"] = s2["id_desa"].astype(str).str.replace(".", "", regex=False)
    s2 = s2.merge(l[["id_str", "LISA"]], on="id_str", how="left").drop(columns="id_str")
    s2["LISA"] = s2["LISA"].map(lab).fillna("Not linked to boundaries")
s2.to_csv("S2_village_cvi_scores.csv", index=False)

with open("revision2_numbers.txt", "w", encoding="utf8") as f:
    for k, (v, lab) in OUT.items():
        f.write(f"[{k}] {lab}: {v}\n")
print("\nWrote revision2_numbers.txt, table4_priority_group.csv, tableS8_bias_checks.csv, "
      "tableS9_class_transitions.csv, S2_village_cvi_scores.csv")
