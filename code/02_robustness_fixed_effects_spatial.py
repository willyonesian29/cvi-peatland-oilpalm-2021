"""
CVI robustness, fixed-effects and spatial analyses (JIMDP #3227)
================================================================
Rebuilds the baseline CVI exactly as in 01_cvi_baseline_analysis.py, then
produces the robustness (Table 3), leave-one-out (Table S5), fixed-effects and
spatial-autocorrelation results. Each number is printed with a tag [R1], [R2], ...
and all tags are written to revision_numbers.txt.

Analyses
  A. Missing-data audit                                    -> [R1]
  B. Weighting / aggregation / normalization sensitivity   -> Table 3, [R2]-[R6]
  C. Indicator-selection sensitivity (leave-one-out)       -> Table 3, [R7]
  D. Exposure-cap sensitivity (drop / winsorize)           -> Table 3
  E. Monte Carlo random weights (1,000 Dirichlet draws)    -> Table 3, [R8]-[R10]
  F. Pillar decomposition of smallholder/industrial gap    -> [R11]-[R13]
  G. Province fixed-effects model, regency-clustered SE    -> [R14]-[R15]
  H. Spatial autocorrelation: global Moran's I + LISA      -> [R16]-[R22]

Inputs
  analysis_input.csv   (from 00d_prepare_analysis_input.py)
  FinalGeoPortal.gpkg  Ina-Geoportal village boundaries (not redistributed — see README);
      SHP_ID_FIELD holds the village code (with dots) matching id_desa.

Requirements
  pip install pandas numpy scipy scikit-learn mapclassify statsmodels \
              geopandas libpysal esda spreg
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, kruskal
from sklearn.decomposition import PCA
from sklearn.preprocessing import minmax_scale

INPUT_PATH = "analysis_input.csv"       # from 00d_prepare_analysis_input.py
VILLAGE_SHP = "FinalGeoPortal.gpkg"      # Ina-Geoportal village boundaries (901 study villages)
SHP_ID_FIELD = "KDEPUM"                  # <-- set village-code field name
KNN_K = 8
N_MC = 1000
SEED = 42
rng = np.random.default_rng(SEED)
OUT = {}


def tag(key, value, label):
    OUT[key] = (value, label)
    print(f"[{key}] {label}: {value}")


# ---------------------------------------------------------------- baseline
SENS = ["air_minum", "air_mck", "jamban", "tps", "konflik", "cemar",
        "bencana", "ragam_eko"]
CAP = ["faskes_km", "nakes", "posyandu", "posyandu_aktif", "bpjs", "sd",
       "smp", "sma", "kursus", "listrik_kk", "sinyal", "inet_desa",
       "inet_warga", "mitigasi", "pusat_dagang", "pasar", "minimarket",
       "logistik", "bank", "kredit", "koperasi", "angkot", "jalan_mobil",
       "jalan_mutu"]
BIN = ["bumdes_bin", "l_tani_bin"]

df = pd.read_csv(INPUT_PATH)
s = df[df["peat_tp_ha"].notna()].drop_duplicates("id_desa", keep="first").copy()
s = s.reset_index(drop=True)

# A. missing-data audit (before imputation)
miss = s[SENS + CAP].isna().mean() * 100
tag("R1", f"{miss.max():.1f}% (max, {miss.idxmax()}); "
          f"{(miss > 0).sum()} of {len(miss)} indicators had any missing value",
    "Missing values imputed with median")

for c in SENS + CAP:
    s[c] = s[c].fillna(s[c].median())
s["bumdes_bin"] = (s["bumdes_tani"] == "Ada").astype(int)
s["l_tani_bin"] = (s["l_tani"] == "Ada").astype(int)

s["exp_raw"] = s["peat_tp_ha"] / s["VillageArea_Ha"]
capped_idx = s.index[s["exp_raw"] > 1]
s["exp_ratio"] = s["exp_raw"].clip(upper=1)


def pillars(d, sens=SENS, cap=CAP + BIN, exposure="minmax"):
    """Return E, S, AC (all 0-1) for data frame d."""
    if exposure == "minmax":
        E = minmax_scale(d["exp_ratio"])
    elif exposure == "rank":
        E = minmax_scale(d["exp_ratio"].rank(pct=True))
    elif exposure == "winsor99":
        E = minmax_scale(d["exp_ratio"].clip(upper=d["exp_ratio"].quantile(.99)))
    S = minmax_scale(1 - d[sens].mean(axis=1))
    AC = minmax_scale(d[cap].mean(axis=1))
    return np.asarray(E), np.asarray(S), np.asarray(AC)


def cvi(E, S, AC, w=(1/3, 1/3, 1/3), agg="arith"):
    V = 1 - AC
    if agg == "arith":
        return w[0]*E + w[1]*S + w[2]*V
    eps = 0.01  # avoids zeros collapsing the geometric mean
    return np.exp(w[0]*np.log(E+eps) + w[1]*np.log(S+eps) + w[2]*np.log(V+eps))


E0, S0, AC0 = pillars(s)
s["CVI"] = cvi(E0, S0, AC0)
s["dom"] = np.where(s["peat_sp_ha"] > s["peat_ip_ha"], "S", "I")          # dominant type of oil palm on peat
import mapclassify
np.random.seed(SEED)
_jb = mapclassify.NaturalBreaks(s["CVI"], k=5)
TOPN = int((s["CVI"] > _jb.bins[2]).sum())     # High + Very High in baseline
print(f"Priority group (High + Very High): {TOPN} villages")
base_rank = s["CVI"].rank(ascending=False)
base_top = set(s.index[base_rank <= TOPN])


def compare(v, idx=None, name=""):
    """Stability metrics of variant v against baseline (same rows)."""
    idx = s.index if idx is None else idx
    b = s.loc[idx, "CVI"].values
    rho = spearmanr(b, v)[0]
    rb = pd.Series(b, index=idx).rank(ascending=False)
    rv = pd.Series(v, index=idx).rank(ascending=False)
    top_b = set(rb.index[rb <= TOPN]); top_v = set(rv.index[rv <= TOPN])
    keep = len(top_b & top_v) / len(top_b) * 100
    dr = np.median(np.abs(rb - rv))
    g = s.loc[idx, "dom"].values
    u, p = mannwhitneyu(v[g == "S"], v[g == "I"])
    direction = "S<I" if np.median(v[g == "S"]) < np.median(v[g == "I"]) else "S>=I"
    return dict(scheme=name, rho=round(rho, 3), top225_retained_pct=round(keep, 1),
                median_abs_rank_shift=round(dr, 1),
                MW_direction=direction, MW_p=f"{p:.2g}")


rows = []
# B. weighting / aggregation / normalization
pcs = np.column_stack([E0, S0, 1 - AC0])
pca = PCA(1).fit(pcs)
pc1 = pca.transform(pcs).ravel()
if np.corrcoef(pc1, s["CVI"])[0, 1] < 0:
    pc1 = -pc1
rows.append(compare(minmax_scale(pc1), name="PCA-derived weights"))
rows.append(compare(cvi(E0, S0, AC0, (.5, .25, .25)),
                    name="Exposure-weighted (0.50/0.25/0.25)"))
rows.append(compare(cvi(E0, S0, AC0, (.2, .4, .4)),
                    name="Socioeconomic-weighted (0.20/0.40/0.40)"))
rows.append(compare(cvi(E0, S0, AC0, agg="geo"),
                    name="Geometric aggregation (equal weights)"))
Er, Sr, ACr = pillars(s, exposure="rank")
rows.append(compare(cvi(Er, Sr, ACr), name="Rank-normalized Exposure"))

# D. exposure cap
Ew, Sw, ACw = pillars(s, exposure="winsor99")
rows.append(compare(cvi(Ew, Sw, ACw), name="Exposure winsorized at 99th pct"))
keep_idx = s.index.difference(capped_idx)
sd = s.loc[keep_idx]
Ed, Sd, ACd = pillars(sd)
rows.append(compare(cvi(Ed, Sd, ACd), idx=keep_idx,
                    name=f"Capped village(s) excluded (n={len(keep_idx)})"))

# C. indicator selection
Eb, Sb, ACb = pillars(s, cap=CAP)
rows.append(compare(cvi(Eb, Sb, ACb), name="Without BUMDes/farmer-group binaries"))
loo = []
for c in SENS:
    E_, S_, A_ = pillars(s, sens=[x for x in SENS if x != c])
    loo.append(compare(cvi(E_, S_, A_), name=f"drop {c}"))
for c in CAP + BIN:
    E_, S_, A_ = pillars(s, cap=[x for x in CAP + BIN if x != c])
    loo.append(compare(cvi(E_, S_, A_), name=f"drop {c}"))
loo = pd.DataFrame(loo)
worst = loo.sort_values("rho").iloc[0]
rows.append(dict(scheme="Leave-one-indicator-out (34 runs), min-max range",
                 rho=f"{loo.rho.min():.3f}-{loo.rho.max():.3f}",
                 top225_retained_pct=f"{loo.top225_retained_pct.min():.1f}-"
                                     f"{loo.top225_retained_pct.max():.1f}",
                 median_abs_rank_shift=f"{loo.median_abs_rank_shift.max():.1f} (max)",
                 MW_direction=f"S<I in {(loo.MW_direction=='S<I').sum()}/34",
                 MW_p=f"max p = {pd.to_numeric(loo.MW_p).max():.2g}"))
tag("R7", f"lowest rho = {worst.rho:.3f} when {worst.scheme}",
    "Most influential single indicator")
loo.to_csv("leave_one_out_detail.csv", index=False)

# E. Monte Carlo random weights
mc = []
for _ in range(N_MC):
    w = rng.dirichlet([1, 1, 1])
    r = compare(cvi(E0, S0, AC0, tuple(w)))
    mc.append((r["rho"], r["top225_retained_pct"], r["MW_direction"],
               float(r["MW_p"])))
mc = pd.DataFrame(mc, columns=["rho", "keep", "dir", "p"])
rows.append(dict(scheme=f"Monte Carlo random weights ({N_MC} draws), median [5th pct]",
                 rho=f"{mc.rho.median():.3f} [{mc.rho.quantile(.05):.3f}]",
                 top225_retained_pct=f"{mc.keep.median():.1f} [{mc.keep.quantile(.05):.1f}]",
                 median_abs_rank_shift="-",
                 MW_direction=f"S<I in {(mc.dir=='S<I').mean()*100:.1f}% of draws",
                 MW_p=f"p<0.05 in {((mc.dir=='S<I')&(mc.p<.05)).mean()*100:.1f}%"))
tag("R8", f"{mc.rho.median():.3f} (5th pct {mc.rho.quantile(.05):.3f})",
    "MC median Spearman rho vs baseline")
tag("R9", f"{mc.keep.median():.1f}% (5th pct {mc.keep.quantile(.05):.1f}%)",
    "MC median retention of baseline top-225 villages")
tag("R10", f"{((mc.dir=='S<I')&(mc.p<.05)).mean()*100:.1f}%",
    "MC share of draws with smallholder<industrial at p<0.05")

tab3 = pd.DataFrame(rows)
tab3.to_csv("table3_robustness.csv", index=False)
print("\nTABLE 3\n", tab3.to_string(index=False))
r = tab3.set_index("scheme")
tag("R2", r.loc["Exposure-weighted (0.50/0.25/0.25)", "rho"], "rho exposure-weighted")
tag("R3", r.loc["Geometric aggregation (equal weights)", "rho"], "rho geometric")
tag("R4", r.loc["Rank-normalized Exposure", "rho"], "rho rank-normalized exposure")
cap_name = [x for x in r.index if x.startswith("Capped")][0]
tag("R5", r.loc[cap_name, "rho"], "rho with capped village excluded")
num = pd.to_numeric(tab3["rho"], errors="coerce")
tag("R6", f"{num.min():.3f}-{num.max():.3f}", "rho range across single alternatives")

# F. pillar decomposition by dominant type
g = s["dom"].values
for k, (arr, nm) in zip(["R11", "R12", "R13"],
                        [(E0, "Exposure"), (S0, "Sensitivity"), (1-AC0, "1-AC")]):
    u, p = mannwhitneyu(arr[g == "S"], arr[g == "I"])
    tag(k, f"median S-dom {np.median(arr[g=='S']):.3f} vs I-dom "
           f"{np.median(arr[g=='I']):.3f}, p = {p:.2g}", f"{nm} by dominant type")

# G. province FE + regency-clustered SE
import statsmodels.formula.api as smf
s["industrial"] = (s["dom"] == "I").astype(int)
m = smf.ols("CVI ~ industrial + C(prov)", data=s).fit(
    cov_type="cluster", cov_kwds={"groups": pd.factorize(s["kab"])[0]})
tag("R14", f"beta = {m.params['industrial']:.3f}, SE = {m.bse['industrial']:.3f}, "
           f"p = {m.pvalues['industrial']:.2g}",
    "Industrial-dominant effect, province FE, regency-clustered SE")
tag("R15", f"{s['kab'].nunique()} regencies", "Number of clusters")
s["resid"] = m.resid

# H. spatial autocorrelation
try:
    import geopandas as gpd
    from libpysal.weights import KNN
    from esda.moran import Moran, Moran_Local
    gdf = gpd.read_file(VILLAGE_SHP)
    gdf[SHP_ID_FIELD] = gdf[SHP_ID_FIELD].astype(str).str.replace(".", "", regex=False)
    s["id_str"] = s["id_desa"].astype(str).str.replace(".", "", regex=False)
    gdf = gdf.merge(s, left_on=SHP_ID_FIELD, right_on="id_str", how="inner")
    print(f"\nMatched {len(gdf)} of {len(s)} villages to boundaries")
    from shapely import force_2d
    gdf["geometry"] = force_2d(gdf.geometry)
    gdf = gdf.to_crs(3395)  # metric CRS for centroids
    pts = np.column_stack([gdf.geometry.centroid.x, gdf.geometry.centroid.y])
    W = KNN.from_array(pts, k=KNN_K); W.transform = "r"
    res = {}
    np.random.seed(SEED)
    for col in ["CVI", "exp_ratio", "resid"]:
        mi = Moran(gdf[col].values, W, permutations=999)
        res[col] = mi
    tag("R16", f"I = {res['CVI'].I:.3f}, pseudo-p = {res['CVI'].p_sim:.3f}",
        f"Global Moran's I, CVI (KNN k={KNN_K})")
    tag("R17", f"I = {res['exp_ratio'].I:.3f}, p = {res['exp_ratio'].p_sim:.3f}",
        "Global Moran's I, Exposure")
    tag("R18", f"I = {res['resid'].I:.3f}, p = {res['resid'].p_sim:.3f}",
        "Moran's I of FE-model residuals")
    lisa = Moran_Local(gdf["CVI"].values, W, permutations=999, seed=SEED)
    sig = lisa.p_sim < 0.05
    q = pd.Series(lisa.q).map({1: "HH", 2: "LH", 3: "LL", 4: "HL"})
    q[~sig] = "ns"
    gdf["LISA"] = q.values
    cnt = gdf["LISA"].value_counts()
    tag("R19", int(cnt.get("HH", 0)), "LISA High-High villages (p<0.05)")
    tag("R20", int(cnt.get("LL", 0)), "LISA Low-Low villages (p<0.05)")
    tag("R21", int(cnt.get("HL", 0) + cnt.get("LH", 0)), "LISA spatial outliers")
    hh = gdf.loc[gdf["LISA"] == "HH"].groupby(["prov", "kab"]).size() \
            .sort_values(ascending=False).head(5)
    tag("R22", "; ".join(f"{k[1]} ({k[0]}, n={v})" for k, v in hh.items()),
        "Regencies with most High-High villages")
    gdf[["id_str", "prov", "kab", "desa", "CVI", "LISA"]].to_csv(
        "lisa_clusters.csv", index=False)
    try:
        from spreg import ML_Error
        X = pd.get_dummies(gdf[["industrial", "prov"]], columns=["prov"],
                           drop_first=True).astype(float)
        me = ML_Error(gdf[["CVI"]].values, X.values, w=W,
                      name_x=list(X.columns), name_y="CVI")
        i = list(X.columns).index("industrial") + 1
        print(f"Spatial error model: industrial beta = {me.betas[i][0]:.3f}, "
              f"p = {me.z_stat[i][1]:.2g}, lambda = {me.lam:.3f}")
    except Exception as e:
        print("spreg spatial error model skipped:", e)
except Exception as e:
    print("\n[H] Spatial analysis skipped:", e)
    print("    Fallback: intraclass correlation of CVI within regencies")
    icc_m = smf.mixedlm("CVI ~ 1", s, groups=s["kab"]).fit()
    va = float(icc_m.cov_re.iloc[0, 0]); ve = icc_m.scale
    print(f"    ICC(regency) = {va/(va+ve):.3f}")

# ---------------------------------------------------------------- write out
with open("revision_numbers.txt", "w", encoding="utf8") as f:
    for k in sorted(OUT, key=lambda x: int(x[1:])):
        v, lab = OUT[k]
        f.write(f"[{k}] {lab}: {v}\n")
print("\nWrote revision_numbers.txt, table3_robustness.csv, "
      "leave_one_out_detail.csv, lisa_clusters.csv")
