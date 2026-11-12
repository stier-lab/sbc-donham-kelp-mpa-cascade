"""Landsat canopy response: report against the rules in CANOPY_PLAN.md (fixed before the full runs).

  python3 code/joint_model/report_canopy.py [validation folder] [output folder]

Reads validation/{c_fit,c_perm,c_sim,c_loo}/*.json, prep_canopy_checks.json and test_multi_real.log; writes
REPORT_canopy.md and validation/check_status_canopy.json. Effects are log changes in canopy density ten years after
designation; intervals are 90%.
"""
import json, sys
from pathlib import Path
import numpy as np

V = Path(sys.argv[1] if len(sys.argv) > 1 else "outputs/joint_model_20261109/validation")
O = Path(sys.argv[2] if len(sys.argv) > 2 else "outputs/joint_model_20261109")
PRI = {"rw": "random walk", "pw": "independent increments"}


def load(p):
    try:
        return json.loads(Path(p).read_text())
    except Exception:
        return None


def ok(r, strict=False):
    if r is None: return False
    if strict: return r["worst_rhat"] <= 1.01 and r["min_neff"] >= 400 and r["divergences"] == 0
    return r["worst_rhat"] <= 1.05 and r["divergences"] <= 0.01 * r["chains"] * r["samples"]


def e(r, k="E10_total"):
    x = r.get(f"can:{k}") if r else None
    return None if x is None else (x["mean"], x["5.0%"], x["95.0%"])


def f(x):
    return "not finished" if x is None else f"{x[0]:+.2f} [{x[1]:+.2f}, {x[2]:+.2f}]"


def conv(r):
    return "—" if r is None else f"R-hat {r['worst_rhat']:.3f}, {r['divergences']} div., ESS {r['min_neff']:.0f}"


fit = lambda n: load(V / "c_fit" / f"can_{n}.json")
status = {}; L = []
chk = load(V / "prep_canopy_checks.json") or {}
log = (V / "test_multi_real.log").read_text() if (V / "test_multi_real.log").exists() else ""
status["C1 independent arithmetic"] = "PASS" if "ALL CHECKS PASSED (1-4)" in log else "INCOMPLETE (fetch out/test_multi_real.log)"

main = {c: fit(f"spring_{c}") for c in PRI}
status["C2 convergence of the main fits"] = "INCOMPLETE" if None in main.values() else "PASS" if all(ok(r, True) for r in main.values()) else "FAIL"
L += ["## Main result: spring canopy (January-June quarters), all depths, flat baseline", "",
      "| Curve prior | Effect at ten years | Effect at four years | Effect at eight years | Quarters without canopy | Converged |", "|---|---|---|---|---|---|"]
for c, cn in PRI.items():
    r = main[c]; g = r["can:g_knots"] if r else None
    k = lambda i: "—" if g is None else f"{g['mean'][i]:+.2f} [{g['5.0%'][i]:+.2f}, {g['95.0%'][i]:+.2f}]"
    L.append(f"| {cn} | {f(e(r))} | {k(0)} | {k(1)} | {'' if r is None else format(r['canopy']['can']['zero_fraction'], '.0%')} | {conv(r)} |")
L += ["", "MPAs designated in 2012 have nine years of canopy data after designation; the ten-year effect leans on the 2003 MPAs and on the curve.", ""]

L += ["## C8. Sensitivities, reported side by side (none chosen by its result)", "",
      "| Variant | Random walk | Independent increments |", "|---|---|---|"]
for lab, n in (("Spring, designation year after (main)", "spring_{c}"), ("Annual maximum (the earlier satellite input)", "max_{c}"),
               ("July-December quarters", "autumn_{c}"), ("All four quarters", "all_{c}"), ("Spring, designation year counted as before", "spring_{c}_next"),
               ("Spring, pixels to 10 m deep", "spring_{c}_to10m"), ("Spring, pixels deeper than 10 m", "spring_{c}_over10m"),
               ("Spring, with pre-trend", "spring_{c}_pre")):
    L.append(f"| {lab} | {f(e(fit(n.format(c='rw'))))} | {f(e(fit(n.format(c='pw'))))} |")
pre = {c: e(fit(f"spring_{c}_pre"), "mu_pretrend") for c in PRI}
L += ["", f"Estimated pre-trend of protected against reference polygons, per decade: random walk {f(pre['rw'])}; independent increments {f(pre['pw'])}.", ""]

s1, s2 = e(main["rw"]), e(fit("spring_rw_s2"))
status["C6 seeds agree"] = "INCOMPLETE" if not (s1 and s2) else "PASS" if abs(s1[0] - s2[0]) <= 0.05 else "FAIL"
L += ["## C6. Seeds", "", f"Seed 1 {f(s1)}; seed 2 {f(s2)}. **{status['C6 seeds agree']}.**", ""]

L += ["## C3. Simulated data on the real design, noise at real-data levels", "",
      "| Scenario | Data sets (usable) | True effect | Mean error | Rule | Result |", "|---|---|---|---|---|---|"]
c3 = {}
for n, lab in (("std", "effect present"), ("null", "no effect")):
    fits = [r for r in (load(p) for p in sorted((V / "c_sim").glob(f"can_sim_{n}_s*.json"))) if r]; use = [r for r in fits if ok(r)]
    if not use:
        L.append(f"| {lab} | {len(fits)} (0) | | | | not finished |"); c3[n] = None; continue
    t = np.array([r["truth"]["can"]["E10_total"] for r in use]); m = np.array([e(r)[0] for r in use]); lo = np.array([e(r)[1] for r in use]); hi = np.array([e(r)[2] for r in use])
    cov = int(((lo <= t) & (hi >= t)).sum()); n_ = len(use)
    passed = cov >= 16 if n == "std" else (n_ - cov) <= 4
    c3[n] = dict(n=n_, fitted=len(fits), covered=cov, bias=float((m - t).mean()), passed=bool(passed), complete=len(fits) >= 20)
    L.append(f"| {lab} | {len(fits)} ({n_}) | {t.mean():+.2f} | {(m - t).mean():+.2f} | "
             f"{'at least 16 of 20 intervals cover the truth' if n == 'std' else 'at most 4 of 20 intervals exclude zero'} | {cov} of {n_} cover |")
status["C3 calibration"] = ("INCOMPLETE" if (None in c3.values() or not all(v["complete"] for v in c3.values())) else "PASS" if all(v["passed"] for v in c3.values()) else "FAIL")
L += ["", f"**{status['C3 calibration']}.**", ""]

L += ["## C4. Shuffled labels", "", "The protected label is given at random to one polygon of each comparison (the MPA or one of its reference areas).", "",
      "| Prior | Shuffles (usable) | Mean of shuffled effects | SD | Mean ÷ SE | 90% intervals covering zero | Real effect | Shuffles at least as extreme | Randomization p |", "|---|---|---|---|---|---|---|---|---|"]
c4 = {}
for c, cn in PRI.items():
    fits = [r for r in (load(p) for p in sorted((V / "c_perm").glob(f"can_perm_{c}_k*.json"))) if r]; use = [r for r in fits if ok(r)]; real = e(main[c])
    if not use:
        L.append(f"| {cn} | not finished | | | | | {f(real)} | | |"); c4[c] = None; continue
    m = np.array([e(r)[0] for r in use]); lo = np.array([e(r)[1] for r in use]); hi = np.array([e(r)[2] for r in use])
    n = len(m); se = m.std(ddof=1) / np.sqrt(n); cv = float(((lo <= 0) & (hi >= 0)).mean()); ext = int((np.abs(m) >= abs(real[0])).sum()) if real else None
    c4[c] = dict(n=n, fitted=len(fits), mean=float(m.mean()), se=float(se), cover_zero=cv, R1a=bool(abs(m.mean()) <= 2 * se), R1b=bool(0.80 <= cv <= 0.97), as_extreme=ext)
    L.append(f"| {cn} | {len(fits)}{'' if len(fits) >= 40 else ' of 40'} ({n}) | {m.mean():+.2f}{'' if c4[c]['R1a'] else ' ✗'} | {m.std(ddof=1):.2f} | {m.mean() / se:+.1f} | "
             f"{cv:.0%}{'' if c4[c]['R1b'] else ' ✗'} | {f(real)} | {'—' if ext is None else f'{ext} of {n}'} | {'—' if ext is None else f'{(1 + ext) / (n + 1):.3f}'} |")
full = [v for v in c4.values() if v and v["fitted"] >= 40]
status["C4 shuffled labels: bias (R1a)"] = "INCOMPLETE" if len(full) < 2 else "PASS" if all(v["R1a"] for v in full) else "FAIL"
status["C4 shuffled labels: interval calibration (R1b)"] = "INCOMPLETE" if len(full) < 2 else "PASS" if all(v["R1b"] for v in full) else "FAIL"
L += ["", "✗ marks a failed rule (mean within 2 standard errors of zero; 80-97% of intervals cover zero). "
      f"**Bias: {status['C4 shuffled labels: bias (R1a)']}. Interval calibration: {status['C4 shuffled labels: interval calibration (R1b)']}.**", ""]

L += ["## C5. Time placebo (years before designation only, fake designation year)", "", "| Fake year | Random walk | Independent increments |", "|---|---|---|"]
c5 = {c: [] for c in PRI}
for y in (1990, 1993, 1996, 1999):
    row = []
    for c in PRI:
        x = e(fit(f"time_{c}_y{y}")); row.append(f(x)); c5[c].append(None if x is None else bool(x[1] <= 0 <= x[2]))
    L.append(f"| {y} | {row[0]} | {row[1]} |")
status["C5 time placebo"] = ("INCOMPLETE" if any(None in v for v in c5.values()) else "PASS" if all(sum(v) >= 3 for v in c5.values()) else "FAIL")
L += ["", f"Intervals covering zero: random walk {sum(bool(x) for x in c5['rw'])} of 4, independent increments {sum(bool(x) for x in c5['pw'])} of 4 "
      f"(rule: at least 3 of 4 each). **{status['C5 time placebo']}.**", ""]

L += ["## C7. Leave one MPA out (random walk)", "", "| MPA left out | Effect at ten years | Shift from the main fit |", "|---|---|---|"]
shifts = []; rev = False
for p in sorted((V / "c_loo").glob("can_loo_*.json")):
    r = load(p)
    if not r: continue
    x = e(r); d = x[0] - s1[0]; shifts.append(d)
    if (s1[1] > 0 or s1[2] < 0) and x[0] * s1[0] < 0: rev = True
    L.append(f"| {', '.join(r['omit_mpa'])} | {f(x)} | {d:+.2f}{' (above 0.25)' if abs(d) > 0.25 else ''} |")
status["C7 leave one MPA out"] = "INCOMPLETE" if len(shifts) < 18 else "FAIL" if rev else "PASS"
L += ["", f"{len(shifts)} of 18 done; largest shift {max(map(abs, shifts)) if shifts else float('nan'):.2f}. **{status['C7 leave one MPA out']}.**", ""]

head = ["# Landsat kelp canopy as its own response", "",
        "Rules and limits on inference were fixed in `CANOPY_PLAN.md` before the full runs. Review only; nothing is adopted. Effects are log "
        "changes in spring canopy density ten years after designation (18 frozen comparisons, 48 polygons, 1984-2021); intervals are 90%.", "",
        "What this response can be used to say: surface canopy biomass in spring over whole reserves and their reference areas. Not plant numbers, "
        "size structure or autumn kelp; weaker below about 10 m (Reed, Rassweiler, Nelson and Kui 2026).", "",
        f"Input check: the quarterly table reproduces Tom Bell's delivered annual-maximum table for {chk.get('non_empty', '?')} non-empty polygons with largest "
        f"difference {chk.get('annual_reconstruction_max_abs_difference', '?')}.", "",
        "## Scorecard", "", "| Check | Result |", "|---|---|"] + [f"| {k} | {v} |" for k, v in status.items()] + [""]
(O / "REPORT_canopy.md").write_text("\n".join(head + L) + "\n")
(V / "check_status_canopy.json").write_text(json.dumps(dict(status=status, C3=c3, C4=c4, C5=c5, loo_shifts=shifts), indent=1, default=str))
for k, v in status.items():
    print(f"{v:<11}{k}")
