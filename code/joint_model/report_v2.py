"""Validation round 2 report (rules in VALIDATION_PLAN_v2.md, fixed before the runs).

  python3 code/joint_model/report_v2.py [validation folder] [output folder]

Reads validation/{v2_perm,v2_split,v2_pre,v2_sim,m_mass_fit,m_fix}/*.json and evidence_summary.json; writes
REPORT_validation_round2.md and validation/check_status_round2.json. Effects are log changes at ten years, 90% intervals.
"""
import json, sys
from pathlib import Path
import numpy as np

V = Path(sys.argv[1] if len(sys.argv) > 1 else "outputs/joint_model_20261109/validation")
O = Path(sys.argv[2] if len(sys.argv) > 2 else "outputs/joint_model_20261109")
TAXA = {"lob": "Lobster", "shp": "Sheephead", "pur": "Purple urchin", "red": "Red urchin", "kelp": "Kelp"}
PRI = {"rw": "random walk", "pw": "independent increments"}


def load(p):
    try:
        return json.loads(Path(p).read_text())
    except Exception:
        return None


def main(g, c):
    """Main fit: sheephead uses the fit that includes KFM roving-diver sizes (m_fix/mx_), as in report_mass.py."""
    return load(V / "m_fix" / f"mx_shp_{c}_s1.json") if g == "shp" else load(V / "m_mass_fit" / f"mm_{g}_{c}_s1.json")


def ok(r):
    return r is not None and r["worst_rhat"] <= 1.05 and r["divergences"] <= 0.01 * r["chains"] * r["samples"]


def e(r, g, k="E10_total"):
    x = r.get(f"{g}:{k}") if r else None
    return None if x is None else (x["mean"], x["5.0%"], x["95.0%"])


def f(x, d=2):
    return "—" if x is None else f"{x[0]:+.{d}f} [{x[1]:+.{d}f}, {x[2]:+.{d}f}]"


def conv(r):
    return "—" if r is None else f"{'yes' if ok(r) else 'NO'} ({r['worst_rhat']:.3f} / {r['divergences']})"


L = ["# Validation round 2: shuffled labels, evidence split, pre-trend", "",
     "Rules were fixed in `VALIDATION_PLAN_v2.md` before the runs. Review only; nothing is adopted. Effects are log changes at "
     "ten years; intervals are 90%.", ""]
status = {}

# ---------------------------------------------------------------- R1
L += ["## R1. Shuffled-label placebo", "",
      "Within each MPA and its reference stations the protected label is dealt at random to the same number of stations; all "
      "stations and surveys are kept. If the model is unbiased the shuffled effects centre on zero. Density only.", "",
      "| Taxon | Prior | Shuffles (usable) | Mean of shuffled effects | SD | Mean ÷ SE | 90% intervals covering zero | Real effect | Shuffles at least as extreme | Randomization p |",
      "|---|---|---|---|---|---|---|---|---|---|"]
r1 = {}
for g, name in TAXA.items():
    for c, cn in PRI.items():
        fits = [load(p) for p in sorted((V / "v2_perm").glob(f"p2_{g}_{c}_k*.json"))]
        fits = [r for r in fits if r is not None]
        use = [r for r in fits if ok(r)]
        real = e(main(g, c), g)
        if not fits:
            L.append(f"| {name} | {cn} | not finished | | | | | {f(real)} | | |"); r1[f"{g}_{c}"] = "INCOMPLETE"; continue
        m = np.array([e(r, g)[0] for r in use]); lo = np.array([e(r, g)[1] for r in use]); hi = np.array([e(r, g)[2] for r in use])
        n = len(m); se = m.std(ddof=1) / np.sqrt(n); cov = float(((lo <= 0) & (hi >= 0)).mean())
        ext = int((np.abs(m) >= abs(real[0])).sum()) if real else None
        p = (1 + ext) / (n + 1) if real else None
        a = abs(m.mean()) <= 2 * se; b = 0.80 <= cov <= 0.97
        done = len(fits) >= 40
        r1[f"{g}_{c}"] = dict(n=n, fitted=len(fits), mean=float(m.mean()), sd=float(m.std(ddof=1)), se=float(se), cover_zero=cov,
                              real=real[0] if real else None, as_extreme=ext, p=p, R1a="PASS" if a else "FAIL", R1b="PASS" if b else "FAIL",
                              complete=done)
        L.append(f"| {name} | {cn} | {len(fits)}{'' if done else ' of 40'} ({n}) | {m.mean():+.2f}{'' if a else ' ✗'} | {m.std(ddof=1):.2f} | {m.mean() / se:+.1f} | "
                 f"{cov:.0%}{'' if b else ' ✗'} | {f(real)} | {ext if real else '—'} of {n} | {p:.3f} |" if real else
                 f"| {name} | {cn} | {len(fits)} ({n}) | {m.mean():+.2f} | {m.std(ddof=1):.2f} | {m.mean() / se:+.1f} | {cov:.0%} | — | — | — |")
L += ["", "✗ marks a failed rule: R1a, the mean of the shuffled effects is within 2 standard errors of zero; R1b, between 80% and 97% of "
      "the shuffled 90% intervals cover zero. The randomization p is reported, not gated.", ""]
full = [v for v in r1.values() if isinstance(v, dict) and v["complete"]]
status["R1a bias"] = ("INCOMPLETE" if len(full) < 10 else "PASS" if all(v["R1a"] == "PASS" for v in full) else "FAIL") if full else "INCOMPLETE"
status["R1b calibration"] = ("INCOMPLETE" if len(full) < 10 else "PASS" if all(v["R1b"] == "PASS" for v in full) else "FAIL") if full else "INCOMPLETE"
fa = [k for k, v in r1.items() if isinstance(v, dict) and v["complete"] and v["R1a"] == "FAIL"]
fb = [k for k, v in r1.items() if isinstance(v, dict) and v["complete"] and v["R1b"] == "FAIL"]
if fa and status["R1a bias"] == "INCOMPLETE": status["R1a bias"] = "FAIL"
if fb and status["R1b calibration"] == "INCOMPLETE": status["R1b calibration"] = "FAIL"
L += [f"**R1a: {status['R1a bias']}**" + (f" (fails: {', '.join(fa)})" if fa else "") + f". **R1b: {status['R1b calibration']}**" + (f" (fails: {', '.join(fb)})" if fb else "") +
      f". {len(full)} of 10 taxon-prior sets complete.", ""]

# ---------------------------------------------------------------- R2
ev = load(V / "evidence_summary.json") or {}
L += ["## R2. Effects by kind of evidence", "",
      "Each taxon fitted on all MPAs (main), on MPAs with before-and-after surveys only, and on post-only MPAs only. Counts are animals "
      "counted at protected stations before designation, at the before-and-after MPAs.", "",
      "| Taxon | Prior | MPAs used | Density | Mean mass | Biomass | Converged (R-hat / divergences) |", "|---|---|---|---|---|---|---|"]
r2 = {}
for g, name in TAXA.items():
    s = ev.get(g, {})
    for c, cn in PRI.items():
        rows = [("all", main(g, c), s.get("mpas")),
                ("before-and-after", load(V / "v2_split" / f"e2_{g}_{c}_ba.json"), s.get("before_after_mpas")),
                ("post-only", load(V / "v2_split" / f"e2_{g}_{c}_post.json"), s.get("post_only_mpas"))]
        for lab, r, nm in rows:
            L.append(f"| {name} | {cn} | {lab} ({nm}) | {f(e(r, g))} | {f(e(r, g, 'E10_mass'))} | {f(e(r, g, 'E10_biomass'))} | {conv(r)} |"
                     if r else f"| {name} | {cn} | {lab} ({nm}) | not finished | | | |")
        ba = e(rows[1][1], g)
        r2[f"{g}_{c}"] = dict(ba=ba, post=e(rows[2][1], g), all=e(rows[0][1], g))
L += ["", "| Taxon | Before-and-after MPAs | Counted before (protected stations) | Largest MPA's share | Of those MPAs, with fewer than 20 counted before | Label under the proposed rule |",
      "|---|---|---|---|---|---|"]
for g, name in TAXA.items():
    s = ev.get(g)
    if not s: continue
    enough = s["before_after_mpas"] - s["before_after_mpas_with_under_20_counted_before"] >= 3
    labs = []
    for c in PRI:
        ba = r2[f"{g}_{c}"]["ba"]
        labs.append(f"{c}: " + ("not finished" if ba is None else "before-and-after supported" if (enough and (ba[1] > 0 or ba[2] < 0)) else
                                "rests on the after-designation contrast"))
    L.append(f"| {name} | {s['before_after_mpas']} of {s['mpas']} | {s['protected_counted_before_at_before_after_mpas']:,} | "
             f"{s['largest_single_mpa_share_of_that']:.0%} | {s['before_after_mpas_with_under_20_counted_before']} | {'; '.join(labs)} |")
L += ["", "Proposed label rule (Adrian to confirm): \"before-and-after supported\" needs the before-and-after-only density interval to exclude "
      "zero and at least three of those MPAs with 20 or more animals counted before designation at protected stations.", ""]
status["R2 evidence split"] = "REPORTED" if all(v["ba"] and v["post"] for v in r2.values()) else "INCOMPLETE"

# ---------------------------------------------------------------- R3
L += ["## R3. Pre-trend sensitivity", "",
      "`--pretrend` lets protected stations at before-and-after MPAs trend against their references before designation and extends that "
      "line through the after years. The flat baseline stays the main specification.", "",
      "### Time placebos (KFM years before the real designation, fake designation year)", "",
      "| Taxon | Fake year | Flat baseline | With pre-trend | Pre-trend per decade | Converged flat / pre-trend |", "|---|---|---|---|---|---|"]
r3 = {}
for g in ("shp", "pur", "red"):
    for y in (1990, 1993, 1996, 1999):
        a = load(V / "v2_pre" / f"t2_{g}_rw_y{y}.json"); b = load(V / "v2_pre" / f"t2_{g}_rw_y{y}_pre.json")
        r3[f"{g}_{y}"] = dict(flat=e(a, g), pre=e(b, g), pre_covers=None if not b else bool(e(b, g)[1] <= 0 <= e(b, g)[2]),
                              flat_covers=None if not a else bool(e(a, g)[1] <= 0 <= e(a, g)[2]))
        L.append(f"| {TAXA[g]} | {y} | {f(e(a, g))} | {f(e(b, g))} | {f(e(b, g, 'mu_pretrend'))} | {conv(a)} / {conv(b)} |")
cs = [r3[f"shp_{y}"]["pre_covers"] for y in (1990, 1993, 1996, 1999)]
status["R3a sheephead time placebo with pre-trend"] = "INCOMPLETE" if None in cs else "PASS" if sum(cs) >= 3 else "FAIL"
fl = [r3[f"shp_{y}"]["flat_covers"] for y in (1990, 1993, 1996, 1999)]
L += ["", f"Sheephead: flat-baseline interval covers zero at {sum(bool(x) for x in fl)} of {sum(x is not None for x in fl)} fake years; with pre-trend at "
      f"{sum(bool(x) for x in cs)} of {sum(x is not None for x in cs)} (rule R3a: at least 3 of 4). **{status['R3a sheephead time placebo with pre-trend']}.**", "",
      "### Real data: flat baseline against pre-trend", "",
      "| Taxon | Prior | Density, flat | Density, pre-trend | Biomass, flat | Biomass, pre-trend | Pre-trend per decade | Depends on the baseline? | Converged |",
      "|---|---|---|---|---|---|---|---|---|"]
dep = []
for g, name in TAXA.items():
    for c, cn in PRI.items():
        a = main(g, c); b = load(V / "v2_pre" / f"t2_{g}_{c}_pre.json")
        if not b:
            L.append(f"| {name} | {cn} | {f(e(a, g))} | not finished | {f(e(a, g, 'E10_biomass'))} | | | | |"); continue
        x, y = e(a, g), e(b, g)
        d = (x[0] * y[0] < 0) or ((y[2] - y[1]) > 2 * (x[2] - x[1]))
        if d: dep.append(f"{name} ({c})")
        L.append(f"| {name} | {cn} | {f(x)} | {f(y)} | {f(e(a, g, 'E10_biomass'))} | {f(e(b, g, 'E10_biomass'))} | {f(e(b, g, 'mu_pretrend'))} | {'YES' if d else 'no'} | {conv(b)} |")
L += ["", "\"Depends on the baseline\" (rule R3b): the two density estimates differ in sign, or the pre-trend interval is more than twice as wide.",
      f"Flagged: {', '.join(dep) if dep else 'none so far'}.", ""]
status["R3b flat vs pre-trend"] = "REPORTED" if all(load(V / "v2_pre" / f"t2_{g}_{c}_pre.json") for g in TAXA for c in PRI) else "INCOMPLETE"

L += ["### Simulated data with a true pre-trend (+0.5 per decade, SD 0.3 across MPAs)", "",
      "| Taxon | Fitted | Data sets (usable) | True effect | Mean error | 90% intervals covering the truth |", "|---|---|---|---|---|---|"]
sim = {}
for g in ("shp", "pur"):
    for tag, lab in (("", "flat baseline"), ("_pre", "with pre-trend")):
        fits = [load(p) for p in sorted((V / "v2_sim").glob(f"s2_{g}_s??{tag}.json"))]; fits = [r for r in fits if r]
        use = [r for r in fits if ok(r)]
        if not use:
            L.append(f"| {TAXA[g]} | {lab} | {len(fits)} (0) | | | |"); continue
        tr = np.array([r["truth"][g]["E10_total"] for r in use]); m = np.array([e(r, g)[0] for r in use])
        cv = float(np.mean([(e(r, g)[1] <= r["truth"][g]["E10_total"] <= e(r, g)[2]) for r in use]))
        sim[f"{g}{tag}"] = dict(n=len(use), fitted=len(fits), bias=float((m - tr).mean()), coverage=cv)
        L.append(f"| {TAXA[g]} | {lab} | {len(fits)} ({len(use)}) | {tr.mean():+.2f} | {(m - tr).mean():+.2f} | {cv:.0%} |")
pre = [v for k, v in sim.items() if k.endswith("_pre")]
status["R3 simulation (pre-trend option recovers the effect)"] = ("INCOMPLETE" if (len(pre) < 2 or any(v["fitted"] < 20 for v in pre)) else
                                                                 "PASS" if all(0.80 <= v["coverage"] <= 0.99 for v in pre) else "FAIL")
L += ["", f"Rule: coverage with the pre-trend option between 80% and 99%. **{status['R3 simulation (pre-trend option recovers the effect)']}.**", ""]

# ---------------------------------------------------------------- A1 lognormal size distribution
L += ["## A1. Size distribution of the mass part: Gamma against lognormal", "",
      "Rules as written for the Gamma family (M7: no contrast statistic outside its 90% replicate interval and at most 10% of all "
      "statistics outside; M8: 85-95% of held-out group means inside the 90% predictive interval; M3: coverage on simulated data at least 0.84).", "",
      "| Taxon | Family | Density | Mean mass | Biomass | Predictive statistics outside | Contrast statistics outside | Hold-out coverage (2 splits) | Converged |",
      "|---|---|---|---|---|---|---|---|---|"]


def ppc_out(r, g):
    p = (r or {}).get("mass_ppc", {}).get(g)
    if not p: return None
    st_ = {k: v for k, v in p.items() if isinstance(v, dict) and "observed" in v}
    out = [k for k, v in st_.items() if not (v["rep_q05"] <= v["observed"] <= v["rep_q95"])]
    return len(out), len(st_), [k for k in out if k.startswith("contrast")]


def hold(files, g):
    c = []
    for f_ in files:
        h = (load(f_) or {}).get("mass_holdout", {})
        h = h.get(g, h) if isinstance(h, dict) else {}
        v = h.get("cover90", h.get("coverage90")) if isinstance(h, dict) else None
        if v is not None: c.append(v)
    return c


a1 = {}
for g, name in TAXA.items():
    for fam, r, hf in (("Gamma", load(V / "v_chk" / f"v_ppc_{g}.json"), [V / "v_chk" / f"v_hold_{g}_s{k}.json" for k in (1, 2)]),
                       ("lognormal", load(V / "v3_logn" / f"l3_{g}_rw_ppc.json"), [V / "v3_logn" / f"l3_{g}_hold_s{k}.json" for k in (1, 2)])):
        po = ppc_out(r, g); hc = hold(hf, g)
        m7 = None if po is None else (not po[2] and po[0] <= 0.10 * po[1]); m8 = None if len(hc) < 2 else all(0.85 <= x <= 0.95 for x in hc)
        a1[f"{g}_{fam}"] = dict(M7=m7, M8=m8, outside=po[0] if po else None, stats=po[1] if po else None, holdout=hc)
        L.append(f"| {name} | {fam} | {f(e(r, g))} | {f(e(r, g, 'E10_mass'))} | {f(e(r, g, 'E10_biomass'))} | "
                 f"{'not finished' if po is None else f'{po[0]} of {po[1]}' + ('' if m7 else ' ✗')} | {'' if po is None else (', '.join(po[2]) or 'none')} | "
                 f"{'not finished' if len(hc) < 2 else ', '.join(f'{x:.0%}' for x in hc) + ('' if m8 else ' ✗')} | {conv(r)} |")
ls = {}
for g in TAXA:
    fits = [r for r in (load(p_) for p_ in sorted((V / "v3_lsim").glob(f"l3_sim_std_{g}_s*.json"))) if r and ok(r)]
    for k in ("E10_mass", "E10_biomass"):
        if fits:
            ls[f"{g}_{k}"] = [bool(e(r, g, k)[1] <= r["truth"][g][k] <= e(r, g, k)[2]) for r in fits]
tot = sum(len(v) for v in ls.values()); cov = sum(sum(v) for v in ls.values()) / tot if tot else None
L += ["", f"Simulated data generated from the Gamma model and fitted with the lognormal one: {tot} intervals (mass and biomass) from "
      f"{sum(len(v) for k, v in ls.items() if k.endswith('E10_mass'))} data sets, coverage {'not finished' if cov is None else f'{cov:.0%}'} (rule: at least 84%; 50 data sets planned).", ""]
for fam in ("Gamma", "lognormal"):
    v7 = [a1[f"{g}_{fam}"]["M7"] for g in TAXA]; v8 = [a1[f"{g}_{fam}"]["M8"] for g in TAXA]
    status[f"A1 {fam}: predictive check (M7)"] = "INCOMPLETE" if None in v7 and all(x in (None, True) for x in v7) else "PASS" if all(v7) else "FAIL"
    status[f"A1 {fam}: hold-out (M8)"] = "INCOMPLETE" if None in v8 and all(x in (None, True) for x in v8) else "PASS" if all(v8) else "FAIL"
status["A1 lognormal: calibration on Gamma-generated data"] = "INCOMPLETE" if tot < 100 else "PASS" if cov >= 0.84 else "FAIL"

# ---------------------------------------------------------------- A2 sheephead local trend
L += ["## A2. Sheephead shuffles with a local trend per MPA cluster", "",
      "The same shuffles as R1 with `--cluster-drift` (a random walk shared by an MPA's protected and reference stations). "
      "Diagnostic; the option is not adopted.", "",
      "| Model | Shuffles (usable) | Mean of shuffled effects | SD | Mean ÷ SE | 90% intervals covering zero | Real effect (same sampler settings) | Shuffles at least as extreme | Randomization p |",
      "|---|---|---|---|---|---|---|---|---|"]
a2 = {}
for lab, pat, real_f in (("without local trend", "v2_perm/p2_shp_rw_k*.json", "v3_clus/c3_shp_rw_real_base.json"),
                         ("with local trend", "v3_clus/c3_shp_rw_k*.json", "v3_clus/c3_shp_rw_real_cluster.json")):
    fits = [r for r in (load(p_) for p_ in sorted(V.glob(pat))) if r]; use = [r for r in fits if ok(r)]
    real = e(load(V / real_f), "shp")
    if not use:
        L.append(f"| {lab} | not finished | | | | | {f(real)} | | |"); continue
    m = np.array([e(r, "shp")[0] for r in use]); lo = np.array([e(r, "shp")[1] for r in use]); hi = np.array([e(r, "shp")[2] for r in use])
    n = len(m); se = m.std(ddof=1) / np.sqrt(n); cv = float(((lo <= 0) & (hi >= 0)).mean())
    ext = int((np.abs(m) >= abs(real[0])).sum()) if real else None
    a2[lab] = dict(n=n, fitted=len(fits), mean=float(m.mean()), se=float(se), cover_zero=cv, R1a=bool(abs(m.mean()) <= 2 * se), R1b=bool(0.80 <= cv <= 0.97),
                   real=real[0] if real else None, as_extreme=ext)
    L.append(f"| {lab} | {len(fits)}{'' if len(fits) >= 40 else ' of 40'} ({n}) | {m.mean():+.2f}{'' if a2[lab]['R1a'] else ' ✗'} | {m.std(ddof=1):.2f} | {m.mean() / se:+.1f} | "
             f"{cv:.0%}{'' if a2[lab]['R1b'] else ' ✗'} | {f(real)} | {'—' if ext is None else f'{ext} of {n}'} | {'—' if ext is None else f'{(1 + ext) / (n + 1):.3f}'} |")
w = a2.get("with local trend")
status["A2 sheephead shuffles with local trend (R1a, R1b)"] = "INCOMPLETE" if (not w or w["fitted"] < 40) else "PASS" if (w["R1a"] and w["R1b"]) else "FAIL"
L += [""]

head = ["## Scorecard", "", "| Check | Result |", "|---|---|"] + [f"| {k} | {v} |" for k, v in status.items()] + [""]
L[4:4] = head
(O / "REPORT_validation_round2.md").write_text("\n".join(L) + "\n")
(V / "check_status_round2.json").write_text(json.dumps(dict(status=status, R1=r1, R3_placebo=r3, simulation=sim, baseline_dependent=dep, A1=a1, A2=a2), indent=1, default=str))
for k, v in status.items():
    print(f"{v:<11}{k}")
