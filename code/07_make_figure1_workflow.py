"""Figure 1. Methodological workflow (JIMDP #3227, revision 2)."""
import textwrap
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

mpl.rcParams.update({"font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "pdf.fonttype": 42, "svg.fonttype": "none"})
W, H = 16.0, 16.35                      # cm
fig = plt.figure(figsize=(W / 2.54, H / 2.54))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

ST = {  # stage: (fill, edge, band)
    "data": ("#e8f0f7", "#5b7896", "#5b7896"),
    "proc": ("#fbf3e2", "#a07b3a", "#a07b3a"),
    "pill": ("#e9f3e6", "#5b8a52", "#5b8a52"),
    "index": ("#f8e6e4", "#a3473f", "#a3473f"),
    "anal": ("#eeeaf6", "#6a5b96", "#6a5b96"),
    "out": ("#ffffff", "#333333", "#333333"),
}
X0, X1 = 1.45, 15.85
GAP = 0.28


def box(x, y, w, h, stage, title, lines, cpc=9.5, tsize=6.6, bsize=6.0, align="left"):
    fc, ec, _ = ST[stage]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.12", fc=fc, ec=ec,
                                lw=0.7, zorder=2))
    ty = y + h - 0.2
    ha, tx = ("left", x + 0.18) if align == "left" else ("center", x + w / 2)
    ax.text(tx, ty, title, fontsize=tsize, fontweight="bold", ha=ha, va="top", zorder=3, color="#1f1f1f")
    body = []
    for ln in lines:
        wrapped = textwrap.wrap(ln, width=int((w - 0.36) * cpc), break_on_hyphens=False) or [""]
        bullet = ln.startswith("• ")
        for k, s in enumerate(wrapped):
            body.append(("   " + s) if (bullet and k > 0) else s)
    ax.text(tx, ty - 0.42, "\n".join(body), fontsize=bsize, ha=ha, va="top", linespacing=1.32, zorder=3,
            color="#2b2b2b")
    return dict(x=x, y=y, w=w, h=h, cx=x + w / 2, top=y + h, bot=y)


def arrow(x0, y0, x1, y1, elbow=None):
    kw = dict(arrowstyle="-|>,head_width=0.16,head_length=0.28", lw=0.8, color="#444444",
              shrinkA=0, shrinkB=0, zorder=1)
    if elbow is None:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=kw)
    else:  # vertical-horizontal-vertical
        ax.plot([x0, x0, x1], [y0, elbow, elbow], color="#444444", lw=0.8, zorder=1, solid_capstyle="butt")
        ax.annotate("", xy=(x1, y1), xytext=(x1, elbow), arrowprops=kw)


def band(y0, y1, num, label, stage):
    c = ST[stage][2]
    ax.add_patch(FancyBboxPatch((0.15, y0), 0.95, y1 - y0, boxstyle="round,pad=0,rounding_size=0.12",
                                fc=c, ec="none", zorder=2))
    ax.text(0.625, (y0 + y1) / 2, f"{num}  {label}", rotation=90, ha="center", va="center", fontsize=6.6,
            fontweight="bold", color="white", zorder=3)


# ------------------------------------------------------------------ 1 data
y1t, h1 = 16.1, 2.1
w4 = (X1 - X0 - 3 * GAP) / 4
d = [box(X0 + i * (w4 + GAP), y1t - h1, w4, h1, "data", t, l) for i, (t, l) in enumerate([
    ("Oil palm map, 2021", ["Descals et al. (2024); Sentinel-1, 10 m", "Industrial vs smallholder; accuracy ≈ 91% and 71–72%"]),
    ("Peatland map", ["BIG, One Map Indonesia", "Scale 1:50,000; no depth classes"]),
    ("Village boundaries", ["Ina-Geoportal (BIG), 1:10,000, 2023 edition", "Kemendagri village codes (KDEPUM)"]),
    ("IDM 2021", ["Village Development Index, Kemendes PDTT", "50 indicators scored 0–5, plus BUMDes and farmer-group items"]),
])]
band(y1t - h1, y1t, 1, "Data", "data")

# ------------------------------------------------------------------ 2 processing
y2t = y1t - h1 - 0.6; h2 = 2.05
wl = 3 * w4 + 2 * GAP
p1 = box(X0, y2t - h2, wl, h2, "proc", "Spatial processing and sample definition (Google Earth Engine, GIS, Python)", [
    "• Earth Engine: oil palm on the 10 m peat mask, vectorized by plantation type; GIS overlay with Ina-Geoportal villages gives oil palm area on peat and village area",
    "• Joined to IDM records by village code: 923 villages",
    "• Cleaning: 22 duplicate codes removed; 3 exposure ratios > 1 capped at 1.0",
    "• Analytic sample: 901 villages in 11 provinces",
], cpc=10.2)
p2 = box(X0 + 3 * (w4 + GAP), y2t - h2, w4, h2, "proc", "Indicator selection", [
    "Conceptual link, clear direction, availability",
    "34 of 50 items kept (S1)",
    "+ 2 binaries",
], cpc=8.4)
band(y2t - h2, y2t, 2, "Processing", "proc")
for i in range(3):
    arrow(d[i]["cx"], d[i]["bot"], d[i]["cx"], p1["top"])
arrow(d[3]["cx"], d[3]["bot"], p2["cx"], p2["top"])

# ------------------------------------------------------------------ 3 pillars
y3t = y2t - h2 - 0.7; h3 = 1.9
w3 = (X1 - X0 - 2 * GAP) / 3
pe_ = box(X0, y3t - h3, w3, h3, "pill", "Exposure (E)", [
    "Oil palm area on peat ÷ village area (extent of drainage-based cultivation)",
    "Min-max scaled 0–1",
])
ps = box(X0 + w3 + GAP, y3t - h3, w3, h3, "pill", "Sensitivity (S): 8 indicators", [
    "Water and sanitation, waste, social conflict, pollution, disaster incidence, economic diversity",
    "Mean, reversed; min-max scaled 0–1",
])
pa = box(X0 + 2 * (w3 + GAP), y3t - h3, w3, h3, "pill", "Adaptive Capacity (AC): 26 items", [
    "Health, education, energy and connectivity, markets, finance, transport, BUMDes, farmer groups",
    "Mean; min-max scaled 0–1",
])
band(y3t - h3, y3t, 3, "Pillars", "pill")
arrow(pe_["cx"], p1["bot"], pe_["cx"], pe_["top"])
arrow(ps["cx"] - 0.6, p1["bot"], ps["cx"] - 0.6, ps["top"])
arrow(pa["cx"], p2["bot"], pa["cx"], pa["top"])
arrow(p2["x"] + 0.45, p2["bot"], ps["cx"] + 0.6, ps["top"], elbow=p2["bot"] - 0.38)

# ------------------------------------------------------------------ 4 index
y4t = y3t - h3 - 0.55; h4 = 1.4
ix = box(X0, y4t - h4, X1 - X0, h4, "index", "Composite Vulnerability Index (CVI) and classification", [
    "CVI = [E + S + (1 − AC)] / 3   (equal weights; baseline)",
    "Jenks natural breaks, k = 5:  Very Low ≤ 0.186  |  Low ≤ 0.277  |  Moderate ≤ 0.363  |  High ≤ 0.467  |  Very High ≤ 0.713",
], cpc=10, align="center", bsize=6.3)
band(y4t - h4, y4t, 4, "Index", "index")
for b in (pe_, ps, pa):
    arrow(b["cx"], b["bot"], b["cx"], ix["top"])

# ------------------------------------------------------------------ 5 analysis
y5t = y4t - h4 - 0.55; h5 = 3.45
a1 = box(X0, y5t - h5, w3, h5, "anal", "Robustness of construction", [
    "• 8 alternatives: PCA weights, two reweightings, geometric aggregation, rank-normalized and winsorized exposure, capped villages excluded, binaries removed",
    "• 34 leave-one-indicator-out runs; 1,000 Monte Carlo weight draws",
    "• Metrics: Spearman rho, priority-group retention, weighted kappa, consensus priority set",
])
a2 = box(X0 + w3 + GAP, y5t - h5, w3, h5, "anal", "Group comparisons", [
    "• Kruskal-Wallis: CVI by province (n ≥ 10)",
    "• Mann-Whitney U: smallholder- vs industrial-dominated villages (oil palm on peat), CVI and each pillar",
    "• OLS with province fixed effects, regency-clustered SE",
    "• Bias checks: mixed villages excluded, two-thirds threshold, area controls",
])
a3 = box(X0 + 2 * (w3 + GAP), y5t - h5, w3, h5, "anal", "Spatial analysis", [
    "• All 901 Ina-Geoportal village polygons",
    "• Global Moran's I and LISA (8 nearest neighbours, 999 permutations)",
    "• Regency intraclass correlation",
    "• Spatial error model with province fixed effects",
])
band(y5t - h5, y5t, 5, "Analysis", "anal")
for b in (a1, a2, a3):
    arrow(b["cx"], ix["bot"], b["cx"], b["top"])

# ------------------------------------------------------------------ 6 outputs
y6t = y5t - h5 - 0.55; h6 = 1.55
ou = box(X0, y6t - h6, X1 - X0, h6, "out", "Outputs", [
    "Village CVI classes and maps (Table 5; Figure 2, Figure S1)   •   Sensitivity to construction choices (Table 3)",
    "Priority-group typology: consensus, exposure-led and development-deficit villages (Table 4)   •   Implications for targeting village development and certification support",
], cpc=10, align="center", bsize=6.2)
band(y6t - h6, y6t, 6, "Outputs", "out")
for b in (a1, a2, a3):
    arrow(b["cx"], b["bot"], b["cx"], ou["top"])

ax.text(X0, y6t - h6 - 0.18, "S = Sensitivity; AC = Adaptive Capacity; IDM = Indeks Desa Membangun; LISA = local indicators of spatial "
        "association; SE = standard error.", fontsize=5.6, color="#555555", ha="left", va="top")

for ext in ("png", "pdf", "svg"):
    fig.savefig(f"Figure1_workflow.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", pad_inches=0.04)
fig.savefig("Figure1_workflow.tif", dpi=600, bbox_inches="tight", pad_inches=0.04, pil_kwargs={"compression": "tiff_lzw"})
print("ok")
