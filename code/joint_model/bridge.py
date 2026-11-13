"""Step-by-step comparison of the paper's Table 2 with the joint model (rules: BRIDGE_PLAN.md). Read-only on the paper.

  python3 code/joint_model/bridge.py [uses parquet] [validation folder] [output folder]

Reads tables/table_02_meta_analysis.csv, data/sumstats_final.csv, data/harmonized/harmonized_response_ratios.csv (paper)
and the joint model's survey table and fits. Writes REPORT_bridge.md and validation/bridge_steps.csv (pooled values only).
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))

USES = sys.argv[1] if len(sys.argv) > 1 else "outputs/integrated_review_20261103/private/supported_measurement_uses.parquet"
V = Path(sys.argv[2] if len(sys.argv) > 2 else "outputs/joint_model_20261109/validation")
O = Path(sys.argv[3] if len(sys.argv) > 3 else "outputs/joint_model_20261109")
ROOT = Path(sys.argv[4]) if len(sys.argv) > 4 else Path(".")
TAXA = {"lob": ("P. interruptus", "Panulirus interruptus", "Spiny lobster"), "shp": ("S. pulcher", "Semicossyphus pulcher", "Sheephead"),
        "pur": ("S. purpuratus", "Strongylocentrotus purpuratus", "Purple urchin"), "red": ("M. franciscanus", "Mesocentrotus franciscanus", "Red urchin"),
        "kelp": ("M. pyrifera", "Macrocystis pyrifera", "Giant kelp")}
rows = []


def pool(est, se):
    """Plain mean, median and DerSimonian-Laird random-effects mean with its standard error."""
    est = np.asarray(est, float); se = np.asarray(se, float); k = len(est)
    if k == 0: return dict(k=0)
    out = dict(k=k, mean=float(est.mean()), median=float(np.median(est)))
    ok = np.isfinite(se) & (se > 0)
    if ok.sum() >= 2:
        y, v = est[ok], se[ok] ** 2; w = 1 / v; mu = (w * y).sum() / w.sum(); Q = (w * (y - mu) ** 2).sum()
        tau2 = max(0.0, (Q - (len(y) - 1)) / (w.sum() - (w ** 2).sum() / w.sum())); w = 1 / (v + tau2)
        out.update(re=float((w * y).sum() / w.sum()), re_se=float(np.sqrt(1 / w.sum())))
    return out


def add(step, label, g, resp, p, note=""):
    if p.get("k", 0):
        rows.append(dict(step=step, what=label, taxon=g, response=resp, histories=p["k"], mean=p.get("mean"), median=p.get("median"),
                         random_effects=p.get("re"), se=p.get("re_se"), note=note))


def three(d, tcol="age", min_before=3, min_after=3):
    """The three uniform quantities for one history: d has columns year, lr, after (bool), age."""
    a, b = d[d.after], d[~d.after]; out = {}
    if len(a) >= min_after:
        out["2a"] = (a.lr.mean(), a.lr.std(ddof=1) / np.sqrt(len(a)))
        if a[tcol].nunique() >= 3:
            x = a[tcol].values - a[tcol].mean(); sl = (x * a.lr.values).sum() / (x ** 2).sum(); res = a.lr.values - a.lr.mean() - sl * x
            out["2c"] = (10 * sl, 10 * np.sqrt((res ** 2).sum() / (len(a) - 2) / (x ** 2).sum()))
        if len(b) >= min_before:
            out["2b"] = (a.lr.mean() - b.lr.mean(), np.sqrt(a.lr.var(ddof=1) / len(a) + b.lr.var(ddof=1) / len(b)))
    return out


# ---------------------------------------------------------------- step 0: Table 2
t2 = pd.read_csv(ROOT / "tables/table_02_meta_analysis.csv")
for g, (short, full, name) in TAXA.items():
    for resp, R in (("Density", "Den"), ("Biomass", "Bio")):
        x = t2[(t2.Taxa == short) & (t2.Response == resp)]
        if len(x):
            rows.append(dict(step="0", what="Paper Table 2 (multilevel model)", taxon=g, response=R, histories=int(x.k.iloc[0]), random_effects=float(x.Estimate.iloc[0]), se=float(x.SE.iloc[0])))

# ---------------------------------------------------------------- step 1: paper's per-history effects
ss = pd.read_csv(ROOT / "data/sumstats_final.csv")
for g, (short, full, name) in TAXA.items():
    for R in ("Den", "Bio"):
        x = ss[(ss.Taxa == short) & (ss.Resp == R)]
        if not len(x): continue
        add("1", "Paper's histories, pooled simply", g, R, pool(x.Mean, x.SE))
        for lab, m in (("mean after (siting included)", x.Model == "Mean"), ("line over after years", (x.Model == "Linear") & (x.BA == "N")),
                       ("before-and-after (step, line or curve)", x.BA == "Y")):
            add("1", f"  of which: {lab}", g, R, pool(x.Mean[m], x.SE[m]))
        for s_ in sorted(x.Source.unique()):
            add("1", f"  of which: {s_}", g, R, pool(x.Mean[x.Source == s_], x.SE[x.Source == s_]))

# ---------------------------------------------------------------- step 2: paper's annual ratios, uniform quantities
rr = pd.read_csv(ROOT / "data/harmonized/harmonized_response_ratios.csv")
start = ss.drop_duplicates("MPA").set_index("MPA").MPA_Start.to_dict()
keep = set(zip(ss.MPA, ss.Taxa, ss.Resp, ss.Source))
LAB = {"2a": "inside against outside, after years", "2b": "after minus before", "2c": "slope over after years x 10"}
for g, (short, full, name) in TAXA.items():
    for R in ("Den", "Bio"):
        acc = {k: [] for k in LAB}
        for (mpa, src), d in rr[(rr.y == full) & (rr.resp == R)].groupby(["CA_MPA_Name_Short", "source"]):
            if (mpa, short, R, src) not in keep or mpa not in start: continue
            d = d[np.isfinite(d.lnDiff)].assign(lr=lambda z: z.lnDiff, after=lambda z: z.BA == "After", age=lambda z: z.year - start[mpa])
            for k, v in three(d).items(): acc[k].append(v)
        for k, v in acc.items():
            if v: add(k, f"Paper's annual ratios: {LAB[k]}", g, R, pool([a for a, _ in v], [b for _, b in v]))

# ---------------------------------------------------------------- step 3: joint model's survey data, same quantities
from data import build
for g, (short, full, name) in TAXA.items():
    D = build(USES, short, regions="three"); c = D["cells"].copy()
    c["dens"] = c.y / (c.area * c.n)
    sy = c.groupby(["MPA", "Source", "station", "status", "year", "T"], as_index=False).agg(y=("y", "sum"), A=("area", lambda a: float((a * c.loc[a.index, "n"]).sum())))
    sy["dens"] = sy.y / sy.A
    c0 = 0.5 * sy.dens[sy.dens > 0].min()
    for mult, tag in ((1.0, ""), (2.0, " (offset doubled)"), (0.5, " (offset halved)")):
        acc = {k: [] for k in LAB}
        for (mpa, src), d in sy.groupby(["MPA", "Source"]):
            w = d.groupby(["year", "status"]).dens.mean().unstack()
            if not {"mpa", "reference"} <= set(w.columns): continue
            w = w.dropna(subset=["mpa", "reference"])
            if not len(w): continue
            T = int(d["T"].iloc[0])
            h = pd.DataFrame(dict(year=w.index, lr=np.log((w.mpa + mult * c0) / (w.reference + mult * c0)).values))
            h["after"] = h.year >= T; h["age"] = h.year - T
            for k, v in three(h).items(): acc[k].append(v)
        for k, v in acc.items():
            if v: add("3" + k[1], f"Joint model's surveys as ratios: {LAB[k]}{tag}", g, "Den", pool([a for a, _ in v], [b for _, b in v]))

# ---------------------------------------------------------------- step 4-5: joint model
def load(p):
    try: return json.loads(Path(p).read_text())
    except Exception: return None


for g in TAXA:
    for c, cn in (("rw", "trend carries forward"), ("pw", "each step near zero")):
        r = load(V / "m_fix" / f"mx_shp_{c}_s1.json") if g == "shp" else load(V / "m_mass_fit" / f"mm_{g}_{c}_s1.json")
        if not r: continue
        q = lambda k: (r[f"{g}:{k}"]["mean"], (r[f"{g}:{k}"]["95.0%"] - r[f"{g}:{k}"]["5.0%"]) / 3.29) if f"{g}:{k}" in r else None
        for k, lab, R in (("mu_lambda", "siting difference at designation (average over reserves, log scale per station; not a ratio of means)", "Den"), ("E10_total", "effect at ten years", "Den"),
                          ("E10_mass", "effect at ten years, mean mass", "Mass"), ("E10_biomass", "effect at ten years", "Bio")):
            v = q(k)
            if v: rows.append(dict(step="4", what=f"Joint model ({cn}): {lab}", taxon=g, response=R, histories=len(r["mpas"][g]) if isinstance(r["mpas"], dict) else None, random_effects=v[0], se=v[1]))
        for ev, lab in (("ba", "before-and-after reserves only"), ("post", "post-only reserves only")):
            e_ = load(V / "v2_split" / f"e2_{g}_{c}_{ev}.json")
            if e_: rows.append(dict(step="5", what=f"Joint model ({cn}): effect at ten years, {lab}", taxon=g, response="Den", histories=len(e_["mpas"][g]) if isinstance(e_["mpas"], dict) else None,
                                    random_effects=e_[f"{g}:E10_total"]["mean"], se=(e_[f"{g}:E10_total"]["95.0%"] - e_[f"{g}:E10_total"]["5.0%"]) / 3.29))

T = pd.DataFrame(rows); T.to_csv(V / "bridge_steps.csv", index=False)
f = lambda x: "" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:+.2f}"
NOTES = (V / "notes_bridge.md").read_text().strip().splitlines() + [""] if (V / "notes_bridge.md").exists() else []
L = ["# From the paper's Table 2 to the joint model, one change at a time", "",
     "Plan and limits: `BRIDGE_PLAN.md`, written before these numbers were computed. Descriptive; nothing is adopted and the paper is unchanged. "
     "All values are log ratios or log changes. \"Pooled\" is the paper's multilevel estimate at step 0, a simple random-effects average at steps 1-3, "
     "and the model's posterior mean at steps 4-5.", ""] + NOTES
for g, (short, full, name) in TAXA.items():
    for R, rn in (("Den", "density"), ("Bio", "biomass"), ("Mass", "mean mass")):
        x = T[(T.taxon == g) & (T.response == R)]
        if not len(x) or (R == "Mass"): continue
        if R == "Bio": x = pd.concat([x, T[(T.taxon == g) & (T.response == "Mass")]])
        L += [f"## {name}, {rn}", "", "| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |", "|---|---|---|---|---|---|---|"]
        for _, r in x.iterrows():
            L.append(f"| {r.step} | {r.what} | {'' if pd.isna(r.histories) else int(r.histories)} | {f(r.get('mean'))} | {f(r.get('median'))} | {f(r.random_effects)} | {'' if pd.isna(r.se) else format(r.se, '.2f')} |")
        L.append("")
(O / "REPORT_bridge.md").write_text("\n".join(L) + "\n")
print(f"{len(T)} rows")
