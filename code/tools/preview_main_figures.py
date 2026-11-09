"""Review-only previews of five main figure roles from the integrated cohort.

These are not manuscript figures. Do not copy them into plots/ or the DOCX
without resolving the cohort and downstream estimands described in README.md.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/main_figure_preview_20261106"
INT = ROOT / "outputs/integrated_review_20261103"
WHOLE = ROOT / "outputs/whole_dataset_biomass_update_20261105"
OUT.mkdir(parents=True, exist_ok=True)

INPUTS = {
    "annual": INT / "annual_values.csv",
    "eligibility": INT / "model_eligibility.csv",
    "series_models": INT / "model_results.csv",
    "joint_series": INT / "joint_model_results.csv",
    "pooled": INT / "pooled_results.csv",
    "heatwave": WHOLE / "heatwave_state_pairs.csv",
    "shapefile": ROOT / "data/external/tom_bell_2026-10-20/02_MPA_CA.shp",
    "land": ROOT / "outputs/dataset_inventory_20261101/site_map_payload.json",
}
for key, path in INPUTS.items():
    if not path.exists():
        raise FileNotFoundError(f"{key}: {path}")

annual = pd.read_csv(INPUTS["annual"])
eligibility = pd.read_csv(INPUTS["eligibility"])
models = pd.read_csv(INPUTS["series_models"])
joint = pd.read_csv(INPUTS["joint_series"])
pooled = pd.read_csv(INPUTS["pooled"])
heat = pd.read_csv(INPUTS["heatwave"])
assert annual.series_id.nunique() == 175
assert int(annual.response_supported.sum()) == 2939
assert set(annual.Source) == {"KFM", "LTER", "PISCO", "Landsat"}

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
    }
)
COL = {
    "P. interruptus": "#D28E26",
    "S. pulcher": "#26799A",
    "S. purpuratus": "#7564A7",
    "M. franciscanus": "#B45F4F",
    "M. pyrifera": "#27765A",
}
LABEL = {
    "P. interruptus": "Lobster",
    "S. pulcher": "Sheephead",
    "S. purpuratus": "Purple urchin",
    "M. franciscanus": "Red urchin",
    "M. pyrifera": "Giant kelp",
}
ORDER = list(LABEL)


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=260, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def clean(ax):
    ax.grid(axis="y", color="#E9E9E9", linewidth=0.55)
    ax.set_axisbelow(True)


# Figure 1: current polygon map for the 22 protected areas represented in
# the corrected annual review. Boundaries are contemporary, not historical.
shapes = gpd.read_file(INPUTS["shapefile"])
names = sorted(annual.MPA.unique())
aliases = {
    "Anacapa Island SMR 2003": "Anacapa Island SMR",
    "Blue Cavern Onshore SMCA": "Blue Cavern Onshore SMCA (No-Take)",
    "Campus Point SMCA": "Campus Point SMCA (No-Take)",
    "Carrington Pt SMR": "Carrington Point SMR",
    "Judith Rk SMR": "Judith Rock SMR",
    "Point Vicente SMCA": "Point Vicente SMCA (No-Take)",
    "Skunk Pt SMR": "Skunk Point SMR",
    "Swamis SMCA": "Swami's SMCA",
}
geo_rows = []
for name in names:
    match = shapes.loc[shapes.NAME.eq(aliases.get(name, name))]
    if len(match) != 1:
        raise AssertionError(f"Unresolved map polygon: {name} ({len(match)} matches)")
    geo_rows.append(match.iloc[0])
selected = gpd.GeoDataFrame(geo_rows, geometry="geometry", crs=shapes.crs)
assert len(selected) == 22
coast_file = OUT / "california_coast.geojson"
if not coast_file.exists():
    r_code = '''library(sf)
utils::data("states50", package="rnaturalearthdata")
ca <- states50[states50$name == "California", ]
sf::st_write(sf::st_transform(ca, 4326), commandArgs(TRUE)[1], quiet=TRUE)'''
    subprocess.run(["Rscript", "-e", r_code, str(coast_file)], check=True)
coast = gpd.read_file(coast_file)
fig, ax = plt.subplots(figsize=(7.1, 4.35))
ax.set_facecolor("#EDF4F7")
coast.plot(ax=ax, facecolor="#F5F0E6", edgecolor="#A8A49E", linewidth=0.45, zorder=1)
for ring in json.loads(INPUTS["land"].read_text())["land_rings"]:
    ax.add_patch(Polygon(ring, closed=True, facecolor="#F5F0E6", edgecolor="#A8A49E", linewidth=0.45, zorder=2))
selected.plot(ax=ax, facecolor="#316D80", edgecolor="white", linewidth=0.6, alpha=0.95, zorder=3)
cent = selected.to_crs(3310).centroid.to_crs(4326)
ax.scatter(cent.x, cent.y, s=15, c="#234A5B", edgecolors="white", linewidths=0.5, zorder=4)
anchors = {
    "Naples SMCA": ("Naples", (8, 9)),
    "Scorpion SMR": ("Scorpion", (-54, -17)),
    "Harris Point SMR": ("Harris Point", (8, 11)),
    "Anacapa Island SMR 2003": ("Anacapa", (-25, -19)),
    "Santa Barbara Island SMR": ("Santa Barbara Island", (8, 9)),
}
for site, (label, offset) in anchors.items():
    point = cent.iloc[names.index(site)]
    ax.annotate(
        label, (point.x, point.y), xytext=offset, textcoords="offset points",
        fontsize=7.2, color="#234A5B", zorder=5,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.7},
    )
ax.set_xlim(-120.75, -117.05)
ax.set_ylim(32.65, 34.65)
ax.set_aspect(1 / np.cos(np.deg2rad(33.7)))
ax.set_xticks([-120.5, -119.5, -118.5, -117.5], ["120.5°W", "119.5°W", "118.5°W", "117.5°W"])
ax.set_yticks([33, 33.5, 34, 34.5], ["33°N", "33.5°N", "34°N", "34.5°N"])
ax.tick_params(length=2)
fig.tight_layout(pad=0.3)
save(fig, "figure_1_site_map_review")

# Figure 2: selected case histories, retaining each source as its own line.
case = ["Scorpion SMR", "Harris Point SMR", "Naples SMCA"]
groups = [
    ["P. interruptus", "S. pulcher"],
    ["S. purpuratus", "M. franciscanus"],
    ["M. pyrifera"],
]
fig, axes = plt.subplots(3, 3, figsize=(7.1, 6.2), sharex=True, sharey="row")
case_rows = []
for c, site in enumerate(case):
    for r, taxa in enumerate(groups):
        ax = axes[r, c]
        block = annual.loc[annual.MPA.eq(site) & annual.Taxa.isin(taxa)].copy()
        block = block.loc[block.Resp.eq("Bio" if r == 2 else "Den") & block.response_supported & np.isfinite(block.ln_ratio)]
        for (source, taxon), q in block.groupby(["Source", "Taxa"]):
            q = q.sort_values("year")
            case_rows.extend(q[["series_id", "year", "Source", "Taxa", "Resp", "ln_ratio"]].to_dict("records"))
            ax.plot(q.year, q.ln_ratio, color=COL[taxon], lw=0.85, alpha=0.67, marker="o", markersize=1.8)
        if not block.empty:
            ax.axvline(int(block.MPA_Start.iloc[0]), color="#777777", ls=(0, (2, 3)), lw=0.9)
        ax.axhline(0, color="#777777", lw=0.6)
        ax.text(0.02, 0.96, chr(97 + r * 3 + c), transform=ax.transAxes, ha="left", va="top", fontweight="bold")
        clean(ax)
        if r == 0:
            ax.set_title(site.replace(" SMR", "").replace(" SMCA", ""), fontsize=9.2)
        if c == 0:
            ax.set_ylabel(["Predator log ratio", "Urchin log ratio", "Kelp log ratio"][r])
        if r == 2:
            ax.set_xlabel("Calendar year")
        ax.set_xlim(1982, 2024)
fig.tight_layout(w_pad=0.65, h_pad=0.65)
save(fig, "figure_2_case_histories_review")
pd.DataFrame(case_rows).to_csv(OUT / "figure_2_plotted_values.csv", index=False)

# Figure 3: individual history estimates plus pooled intervals for the two
# distinct estimands. Zero has different interpretation in each panel.
series_lookup = annual[["series_id", "Taxa", "Resp"]].drop_duplicates()
model_points = joint.loc[joint.analysis.eq("estimated_AR1")].copy()
assert not model_points.duplicated(["series_id", "model", "term"]).any()
specs = [
    ("flat_before_after", "average_observed_post_change"),
    ("post_trend", "recovery_slope"),
]
fig, axes = plt.subplots(1, 2, figsize=(7.1, 4.35), sharey=True)
rows = [(t, resp) for t in ORDER for resp in (["Den", "Bio"] if t != "M. pyrifera" else ["Bio"])]
ys = {key: len(rows) - 1 - i for i, key in enumerate(rows)}
for panel, ((model, term), ax) in enumerate(zip(specs, axes)):
    sub = model_points.loc[model_points.model.eq(model) & model_points.term.eq(term)]
    pool = pooled.loc[pooled.Source.eq("All programs") & pooled.model.eq(model) & pooled.term.eq(term)]
    for i, (taxon, resp) in enumerate(rows):
        y = ys[(taxon, resp)]
        q = sub.loc[sub.Taxa.eq(taxon) & sub.Resp.eq(resp)]
        ax.scatter(q.estimate, np.full(len(q), y) + ((np.arange(len(q)) % 7) - 3) * 0.035,
                   s=8, color=COL[taxon], alpha=0.23, edgecolors="none", zorder=2)
        p = pool.loc[pool.Taxa.eq(taxon) & pool.Resp.eq(resp)]
        if len(p):
            z = p.iloc[0]
            ax.plot([z.low, z.high], [y, y], lw=2.4, color=COL[taxon], solid_capstyle="round", zorder=4)
            ax.scatter([z.estimate], [y], s=40, color=COL[taxon], edgecolors="white", linewidths=0.7, zorder=5)
    ax.axvline(0, color="#666666", lw=0.8)
    ax.set_yticks(list(ys.values()), [f"{LABEL[t]}  {'biomass' if r == 'Bio' else 'density'}" for t, r in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Before–after log ratio change" if panel == 0 else "Annual log ratio change")
    ax.text(0.01, 1.02, "a" if panel == 0 else "b", transform=ax.transAxes, fontweight="bold")
    ax.grid(axis="y", color="#EEEEEE", linewidth=0.5)
    ax.set_axisbelow(True)
axes[0].set_xlim(-3.2, 3.2)
axes[1].set_xlim(-0.52, 0.52)
fig.tight_layout(w_pad=1.5)
save(fig, "figure_3_separate_estimands_review")

# Figure 4: observed biomass ratios by elapsed calendar year. The thick
# summary is a descriptive annual median, not a fitted recovery model.
fig = plt.figure(figsize=(7.1, 6.25))
grid = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1], hspace=0.54, wspace=0.21)
axes = [fig.add_subplot(grid[i // 2, i % 2]) for i in range(4)] + [fig.add_subplot(grid[2, :])]
for i, taxon in enumerate(ORDER):
    ax = axes[i]
    q = annual.loc[annual.Taxa.eq(taxon) & annual.Resp.eq("Bio") & annual.response_supported & np.isfinite(annual.ln_ratio)].copy()
    q = q.loc[q.year.ge(q.MPA_Start)].copy()
    q["elapsed"] = q.year - q.MPA_Start
    for _, g in q.groupby("series_id"):
        g = g.sort_values("elapsed")
        ax.plot(g.elapsed, g.ln_ratio, lw=0.7, alpha=0.23, color=COL[taxon])
        ax.scatter(g.elapsed, g.ln_ratio, s=3.3, alpha=0.28, color=COL[taxon], edgecolors="none")
    if not q.empty:
        med = q.groupby("elapsed").ln_ratio.agg(["median", "count"]).reset_index()
        med = med[med["count"].ge(3)]
        ax.plot(med.elapsed, med["median"], lw=2, color=COL[taxon])
    ax.axhline(0, color="#777777", lw=0.6)
    ax.text(0.02, 0.96, chr(97 + i), transform=ax.transAxes, ha="left", va="top", fontweight="bold")
    ax.set_title(LABEL[taxon], fontsize=9)
    ax.set_xlim(0, 23)
    ax.set_ylim(-5.5, 7)
    clean(ax)
fig.supylabel("Log protected/reference biomass", x=0.02, fontsize=9)
fig.supxlabel("Calendar years since designation", y=0.01, fontsize=9)
fig.subplots_adjust(left=0.11, right=0.99, bottom=0.09, top=0.96)
save(fig, "figure_4_observed_biomass_review")

# Figure 5: source-specific paired changes relative to each side's own
# 2010–2013 kelp baseline. Individual source rows stay separate.
kelp = heat.loc[heat.Taxa.eq("M. pyrifera") & heat.Source.ne("Landsat") & heat.period.isin(["during", "recent"])].copy()
kelp = kelp.loc[(kelp.inside_ratio > 0) & (kelp.outside_ratio > 0)]
sites = sorted(set(zip(kelp.MPA, kelp.Source)))
fig, axes = plt.subplots(1, 2, figsize=(7.1, 5.15), sharey=True)
for j, period in enumerate(["during", "recent"]):
    ax = axes[j]
    z = kelp.loc[kelp.period.eq(period)]
    for i, (site, source) in enumerate(sites):
        row = z.loc[z.MPA.eq(site) & z.Source.eq(source)]
        if len(row) == 0:
            continue
        a, b = float(row.inside_ratio.iloc[0]), float(row.outside_ratio.iloc[0])
        y = len(sites) - 1 - i
        ax.plot([a, b], [y, y], color="#A9A9A9", lw=0.75, zorder=1)
        ax.scatter([a], [y], color="#23769A", s=18, zorder=3)
        ax.scatter([b], [y], color="#D28E26", s=18, zorder=3)
    ax.axvline(1, color="#777777", lw=0.8, ls=(0, (3, 3)))
    ax.set_xscale("log")
    ax.set_xlim(0.004, 500)
    ax.set_xticks([0.01, 0.1, 1, 10, 100], ["0.01", "0.1", "1", "10", "100"])
    ax.set_xlabel("Biomass / 2010–2013 baseline")
    ax.text(0.02, 1.02, "a" if j == 0 else "b", transform=ax.transAxes, fontweight="bold")
    ax.text(0.5, 1.02, "2014–2016" if j == 0 else "2020–2023", transform=ax.transAxes, ha="center")
    ax.grid(axis="x", color="#EEEEEE", linewidth=0.55)
axes[0].set_yticks(range(len(sites)), [f"{s.replace(' SMR', '').replace(' SMCA', '')} · {src}" for s, src in sites[::-1]])
axes[0].tick_params(axis="y", length=0)
fig.tight_layout(w_pad=1.0)
save(fig, "figure_5_kelp_baseline_pairs_review")
kelp.to_csv(OUT / "figure_5_plotted_values.csv", index=False)

INPUTS["coast"] = coast_file
manifest = {key: {"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for key, path in INPUTS.items()}
(OUT / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(f"Wrote five review previews to {OUT}")
