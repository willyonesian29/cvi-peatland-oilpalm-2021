"""Figure 2 (CVI classes) and Figure S1 (LISA clusters) for JIMDP #3227.

Run 05_prepare_map_layers.py first. Background layers are Natural Earth 1:10m
(public domain; https://www.naturalearthdata.com): Admin 0 - Countries and
Admin 1 - States, Provinces. Unzip them next to this script.

Requirements: geopandas, matplotlib, matplotlib-scalebar
"""
import geopandas as gpd, pandas as pd, numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle, FancyArrow
from matplotlib.ticker import FuncFormatter
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib.patheffects as pe

mpl.rcParams.update({"font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.4, "ytick.major.width": 0.4, "xtick.major.size": 2,
                     "ytick.major.size": 2, "pdf.fonttype": 42, "svg.fonttype": "none"})

SEA, LAND, LAND_EDGE = "#eaf1f6", "#e4e2dc", "#ffffff"
NEIGH = "#efede8"
v = gpd.read_file("villages.gpkg")
v["geometry"] = v.geometry.simplify(0.0005, preserve_topology=True)   # ~50 m; invisible at print scale
NE0 = "ne_10m_admin_0_countries.shp"
NE1 = "ne_10m_admin_1_states_provinces.shp"
ne = gpd.read_file(NE0)
ne.columns = [c.lower() if c.lower() in ("name", "admin") else c for c in ne.columns]
idn = ne[ne["admin"] == "Indonesia"]
neigh = ne[ne["admin"] != "Indonesia"]
adm1 = gpd.read_file(NE1)
adm1 = adm1[adm1["admin"] == "Indonesia"]
# keep only the map area and simplify backgrounds (keeps vector files small; invisible at print scale)
_bbox = (94.0, -7.0, 120.0, 7.0)
neigh = neigh.clip(_bbox); idn = idn.clip(_bbox); adm1 = adm1.clip(_bbox)
for _g in (neigh, idn, adm1):
    _g["geometry"] = _g.geometry.simplify(0.002, preserve_topology=True)
# regency label anchors: union of study villages per regency
kab = v.dissolve("kab").reset_index()[["kab", "geometry"]]
kab["nmkab"] = kab["kab"]

EXT = {"a": (94.6, 118.7, -6.3, 5.6), "b": (101.15, 103.75, -0.85, 2.30), "c": (108.55, 115.25, -3.85, 2.05)}

CVI_ORDER = ["Very Low", "Low", "Moderate", "High", "Very High"]
CVI_COL = dict(zip(CVI_ORDER, ["#fff3b0", "#fdc35a", "#f58a3c", "#d7301f", "#7f0000"]))
CVI_RNG = {"Very Low": "0.038–0.186", "Low": ">0.186–0.277", "Moderate": ">0.277–0.363",
           "High": ">0.363–0.467", "Very High": ">0.467–0.713"}
LISA_ORDER = ["High-High", "Low-Low", "High-Low", "Low-High", "Not significant"]
LISA_COL = dict(zip(LISA_ORDER, ["#c51b2a", "#2b6cb0", "#f4a582", "#92c5de", "#ffffff"]))

lon_fmt = FuncFormatter(lambda x, _: f"{abs(x):.0f}°{'E' if x >= 0 else 'W'}")
lat_fmt = FuncFormatter(lambda y, _: "0°" if abs(y) < 1e-9 else f"{abs(y):.0f}°{'N' if y > 0 else 'S'}")
halo = [pe.withStroke(linewidth=1.6, foreground="white")]


def base(ax, key, detail):
    x0, x1, y0, y1 = EXT[key]
    ax.set_facecolor(SEA)
    neigh.plot(ax=ax, color=NEIGH, edgecolor="#c9c6bf", linewidth=0.3, zorder=1)
    idn.plot(ax=ax, color=LAND, edgecolor="none", zorder=1.2)
    adm1.plot(ax=ax, facecolor="none", edgecolor=LAND_EDGE, linewidth=0.6 if key == "a" else 0.9, zorder=2)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
    step = 5 if key == "a" else 1
    ax.set_xticks(np.arange(np.ceil(x0 / step) * step, x1, step))
    ax.set_yticks(np.arange(np.ceil(y0 / step) * step, y1 + 1e-6, step if key != "a" else 2.5))
    if key == "a":
        ax.set_yticks([-5, -2.5, 0, 2.5, 5]); ax.yaxis.set_major_formatter(
            FuncFormatter(lambda y, _: "0°" if y == 0 else f"{abs(y):g}°{'N' if y > 0 else 'S'}"))
    else:
        ax.yaxis.set_major_formatter(lat_fmt)
    ax.xaxis.set_major_formatter(lon_fmt)
    ax.tick_params(labelsize=6, direction="out", pad=1.5)
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.grid(color="#ffffff", linewidth=0.3, zorder=0.5)
    ax.add_artist(ScaleBar(111320, "m", location="lower right", length_fraction=0.22, height_fraction=0.012,
                           box_alpha=0.85, font_properties={"size": 6}, border_pad=0.6, sep=2,
                           scale_loc="top", fixed_units="km" if key != "a" else None))
    ax.text(0.012, 0.985, key, transform=ax.transAxes, ha="left", va="top", fontsize=9, fontweight="bold",
            path_effects=halo, zorder=10)


def north(ax, x=0.965, y=0.86):
    ax.annotate("", xy=(x, y + 0.09), xytext=(x, y), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>,head_width=0.25,head_length=0.5", lw=0.8, color="#333"))
    ax.text(x, y + 0.105, "N", transform=ax.transAxes, ha="center", va="bottom", fontsize=7, fontweight="bold")


NUDGE = {  # (dx, dy) label offsets in degrees from the centroid of the study villages, to avoid overlaps
    "KAPUAS": (0.09, 0.60), "PULANG PISAU": (0.15, 0.77), "KATINGAN": (-0.09, 0.90),
    "KOTAWARINGIN TIMUR": (-0.45, 0.46),
}


def labels(ax, key):
    if key == "a":
        L = {"Aceh": (96.6, 4.4), "Sumatera Utara": (98.9, 2.6), "Riau": (101.0, -0.35),
             "Jambi": (102.9, -1.85), "Sumatera Selatan": (104.3, -3.6), "Kalimantan\nBarat": (110.6, 0.55),
             "Kalimantan\nTengah": (113.4, -1.05), "Kalimantan\nTimur": (116.4, 0.6),
             "Kalimantan\nUtara": (116.4, 3.2), "Kalimantan\nSelatan": (115.6, -3.3),
             "MALAYSIA": (103.0, 4.4), "MALAYSIA ": (113.6, 3.4)}
        for t, (x, y) in L.items():
            up = t.strip().isupper()
            ax.text(x, y, t.strip(), fontsize=5.6 if not up else 5.4, color="#6f6b63" if not up else "#9a968e",
                    style="italic" if not up else "normal", ha="center", va="center", path_effects=halo,
                    zorder=6, linespacing=0.95)
    else:
        names = {"b": ["BENGKALIS", "SIAK", "PELALAWAN", "INDRAGIRI HILIR", "INDRAGIRI HULU", "KEPULAUAN MERANTI",
                       "KOTA DUMAI", "ROKAN HILIR"],
                 "c": ["KAPUAS", "KOTAWARINGIN TIMUR", "PULANG PISAU", "KATINGAN", "KUBU RAYA", "KETAPANG",
                       "SINTANG", "KAPUAS HULU", "SERUYAN", "KOTAWARINGIN BARAT", "MEMPAWAH", "SAMBAS",
                       "SANGGAU", "MELAWI"]}[key]
        nudge = NUDGE
        x0, x1, y0, y1 = EXT[key]
        for n in names:
            r = kab[kab.nmkab == n]
            if r.empty:
                continue
            p = r.geometry.iloc[0].centroid
            dx, dy = nudge.get(n, (0, 0)); x, y = p.x + dx, p.y + dy
            if x0 + 0.1 < x < x1 - 0.1 and y0 + 0.1 < y < y1 - 0.1:
                ax.text(x, y, n.title().replace("Kota ", ""), fontsize=5.4, color="#4d4a44", style="italic",
                        ha="center", va="center", path_effects=halo, zorder=7)


def draw(kind, out):
    col, order = (CVI_COL, CVI_ORDER) if kind == "cvi" else (LISA_COL, LISA_ORDER)
    field = "CVI_class" if kind == "cvi" else "LISA"
    w_cm = 16.0
    fig = plt.figure(figsize=(w_cm / 2.54, 16.4 / 2.54))
    axa = fig.add_axes([0.055, 0.475, 0.93, 0.515])
    axb = fig.add_axes([0.055, 0.035, 0.335, 0.405])
    axc = fig.add_axes([0.445, 0.035, 0.54, 0.405])
    for ax, key in ((axa, "a"), (axb, "b"), (axc, "c")):
        base(ax, key, detail=key != "a")
        vv = v if kind == "cvi" else v[v[field].isin(order)]
        for cls in order[::-1]:
            sub = vv[vv[field] == cls]
            if len(sub):
                ec = "#5c5a55" if key != "a" else (("#8a877f") if col[cls] == "#ffffff" else col[cls])
                sub.plot(ax=ax, color=col[cls], edgecolor=ec,
                         linewidth=0.12 if key != "a" else (0.15 if col[cls] == "#ffffff" else 0.25), zorder=3 + (order.index(cls) == 0) * 0)
        labels(ax, key)
        ax.set_xlabel(""); ax.set_ylabel("")
    north(axa)
    for key, lab in (("b", "b"), ("c", "c")):
        x0, x1, y0, y1 = EXT[key]
        axa.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, lw=0.7, ec="#222", zorder=8))
        axa.text(x0 + 0.12, y1 - 0.12, lab, fontsize=7, fontweight="bold", ha="left", va="top",
                 path_effects=halo, zorder=9)
    counts = v[field].value_counts()
    if kind == "cvi":
        hand = [Patch(fc=col[c], ec="#5c5a55", lw=0.3, label=f"{c}: {CVI_RNG[c]} ({counts.get(c, 0)})")
                for c in order]
        title = "CVI class (villages)"
    else:
        hand = [Patch(fc=col[c], ec="#5c5a55", lw=0.3, label=f"{c} ({counts.get(c, 0)})") for c in order]
        title = "LISA cluster, p < 0.05 (villages)"
    hand.append(Patch(fc=LAND, ec="#bdbab3", lw=0.3, label="Other land"))
    leg = axa.legend(handles=hand, title=title, loc="lower left", fontsize=5.6, title_fontsize=5.9,
                     frameon=True, framealpha=0.92, edgecolor="#bdbab3", borderpad=0.6, handlelength=1.4,
                     handleheight=0.9, labelspacing=0.3, bbox_to_anchor=(0.006, 0.012), handletextpad=0.5)
    leg.get_frame().set_linewidth(0.4); leg._legend_box.align = "left"
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{out}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(f"{out}.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw"}, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(f"{out}.jpg", dpi=400, pil_kwargs={"quality": 95}, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


draw("cvi", "Figure2_CVI_classes")
draw("lisa", "FigureS1_LISA_clusters")
print("done", len(v), v.CVI_class.value_counts().to_dict(), v.LISA.value_counts().to_dict())
