"""Validation report for the joint density, mean mass and biomass model.

  python3 code/joint_model/report_mass_validation.py

Applies the rules of outputs/joint_model_20261109/MASS_VALIDATION_PLAN.md (written before the runs) to
validation/{m_mass_fit,m_fix,v_sim,v_chk,v_big} and writes REPORT_mass_validation.md and
validation/check_status_mass_validation.json. A check with missing runs is reported as incomplete, never as passed.
"""
import glob, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("outputs/joint_model_20261109")
V = ROOT / "validation"
NAME = {"lob": "Lobster", "shp": "Sheephead", "pur": "Purple urchin", "red": "Red urchin", "kelp": "Kelp"}
ORDER = ["lob", "shp", "pur", "red", "kelp"]
QS = [("E10_total", "density"), ("E10_mass", "mean mass"), ("E10_biomass", "biomass")]
NOTES = json.loads((V / "notes_mass_validation.json").read_text()) if (V / "notes_mass_validation.json").exists() else {}


def loads(*gates):
    out = {}
    for gate in gates:
        for p in sorted(glob.glob(str(V / gate / "*.json"))):
            try:
                out[Path(p).stem] = json.loads(Path(p).read_text())
            except Exception:
                pass
    return out


def iv(d):
    return f"{d['mean']:+.2f} [{d['5.0%']:+.2f}, {d['95.0%']:+.2f}]"


def excl(d):
    return d["5.0%"] > 0 or d["95.0%"] < 0


def ok(r):
    return r["worst_rhat"] < 1.01 and r["min_neff"] > 400 and r["divergences"] == 0


def conv(r):
    return ("yes" if ok(r) else "NO") + f" ({r['worst_rhat']:.3f} / {r['min_neff']:.0f} / {r['divergences']})"


def two_se(p, n):
    return 2 * np.sqrt(p * (1 - p) / max(n, 1))


def main_fit(F, g, c="rw", seed=1):
    return F.get(f"mx_{g}_{c}_s{seed}") if g == "shp" else F.get(f"mm_{g}_{c}_s{seed}")


def verdict(status, key, passed, complete=True):
    status[key] = "PASS" if (passed and complete) else ("INCOMPLETE" if not complete else "FAIL")
    return f"**{status[key]}.**"


def main():
    F = loads("m_mass_fit", "m_fix"); S = loads("v_sim"); C = loads("v_chk"); B = loads("v_big")
    status = {}
    L = ["# Is the density, mass and biomass model sound? Validation report", "",
         "Rules were fixed in `MASS_VALIDATION_PLAN.md` before the runs. Review only; nothing is adopted. "
         "Effects are log changes at ten years; intervals are 90%.", ""]
    L += NOTES.get("summary", []) + [""]
    body = []

    # M1, M2
    body += ["## M1. Convergence", "", "| Taxon | Random-walk prior, seed 1 | seed 2 | Independent-increment prior |", "|---|---|---|---|"]
    allr = []
    for g in ORDER:
        rs = [main_fit(F, g, "rw", 1), main_fit(F, g, "rw", 2), main_fit(F, g, "pw", 1)]; allr += rs
        body.append(f"| {NAME[g]} | " + " | ".join(conv(r) if r else "not finished" for r in rs) + " |")
    done = all(r is not None for r in allr)
    body += ["", "(R-hat / effective draws / divergent transitions.)", "", verdict(status, "M1 convergence", done and all(ok(r) for r in allr if r), done), ""]
    body += ["## M2. Reproducibility", ""]
    d = {}
    for g in ORDER:
        a, b = main_fit(F, g, "rw", 1), main_fit(F, g, "rw", 2)
        if a and b:
            d[g] = max(abs(a[f"{g}:{q}"]["mean"] - b[f"{g}:{q}"]["mean"]) for q, _ in QS)
    body += ["Largest difference between two independent fits: " + "; ".join(f"{NAME[g]} {v:.3f}" for g, v in d.items()) + ".", "",
             verdict(status, "M2 seeds agree", bool(d) and max(d.values()) < 0.05, len(d) == len(ORDER)), ""]

    # M3, M4
    def sim_table(scs):
        rows, pooled = [], {}
        for sc in scs:
            for g in ORDER:
                rs = [r for k, r in S.items() if k.startswith(f"v_sim_{sc}_{g}_")]
                if not rs:
                    continue
                row = [sc, NAME[g], str(len(rs)), str(sum(r["worst_rhat"] < 1.05 for r in rs))]
                for q in ("E10_mass", "E10_biomass"):
                    tr = np.array([r["truth"][g][q] for r in rs]); est = np.array([r[f"{g}:{q}"]["mean"] for r in rs])
                    lo = np.array([r[f"{g}:{q}"]["5.0%"] for r in rs]); hi = np.array([r[f"{g}:{q}"]["95.0%"] for r in rs])
                    c = (lo <= tr) & (tr <= hi); ex = (lo > 0) | (hi < 0)
                    pooled.setdefault((sc, q), dict(c=[], ex=[], bias={}))
                    pooled[(sc, q)]["c"] += c.tolist(); pooled[(sc, q)]["ex"] += ex.tolist(); pooled[(sc, q)]["bias"][g] = float((est - tr).mean())
                    row += [f"{tr.mean():+.2f}", f"{(est - tr).mean():+.2f}", f"{c.sum()}/{len(rs)}", f"{ex.sum()}/{len(rs)}"]
                rows.append("| " + " | ".join(row) + " |")
        head = ["| Scenario | Taxon | Data sets | R-hat < 1.05 | Mass: truth | mean error | covered | excludes zero | Biomass: truth | mean error | covered | excludes zero |", "|" + "---|" * 12]
        return head + rows, pooled

    body += ["## M3. Calibration when the model is right", "",
             "Simulated on each taxon's real survey and measurement design, then refitted (2 chains, 1,500 + 1,500).", ""]
    t, pooled = sim_table(["std", "null", "trend"]); body += t + [""]
    okk, comp, msgs = True, True, []
    for sc in ("std", "null", "trend"):
        for q in ("E10_mass", "E10_biomass"):
            p = pooled.get((sc, q))
            if not p or len(p["c"]) < 100:
                comp = False
            if p:
                n = len(p["c"]); c = float(np.mean(p["c"])); thr = 0.90 - two_se(.9, n); okk &= c >= thr
                msgs.append(f"{sc}, {q.replace('E10_', '')}: coverage {c:.2f} of {n} (rule ≥ {thr:.2f})")
        p = pooled.get((sc, "E10_mass"))
        if p and sc in ("null", "trend"):
            n = len(p["ex"]); x = float(np.mean(p["ex"])); thr = 0.10 + two_se(.1, n); okk &= x <= thr
            msgs.append(f"{sc}: {100 * x:.0f}% of mass intervals exclude zero (rule ≤ {100 * thr:.0f}%)")
        if p:
            bad = [f"{NAME[g]} {v:+.2f}" for g, v in p["bias"].items() if abs(v) >= 0.10]
            okk &= not bad
            if bad:
                msgs.append(f"{sc}: mean error of the mass effect at or above 0.10 for " + ", ".join(bad))
    body += ["- " + m for m in msgs] + ["", verdict(status, "M3 calibration", okk, comp), ""]

    body += ["## M4. Robustness when an assumption is wrong", "",
             "`bimodal`: two size cohorts instead of a Gamma. `visitnoise`: extra visit-to-visit noise in mean size. "
             "`coupled`: mean size falls with station density while the model assumes independence.", ""]
    t, pooled = sim_table(["bimodal", "visitnoise", "coupled"]); body += t + [""]
    for sc in ("bimodal", "visitnoise", "coupled"):
        okk, comp, msgs = True, True, []
        for q in ("E10_mass", "E10_biomass"):
            p = pooled.get((sc, q))
            if not p or len(p["c"]) < 100:
                comp = False
            if p:
                n = len(p["c"]); c = float(np.mean(p["c"])); thr = 0.90 - two_se(.9, n); okk &= c >= thr
                msgs.append(f"{q.replace('E10_', '')}: coverage {c:.2f} of {n} (rule ≥ {thr:.2f})")
        body += [f"`{sc}` — " + "; ".join(msgs) + ". " + verdict(status, f"M4 robust to {sc}", okk, comp), ""]

    # M5
    body += ["## M5. Placebo in time", "", "KFM surveys before the real designation with a fake designation in 1992. A sound model finds nothing.", "",
             "| Taxon | Prior | Density | Mean mass | Biomass | Converged |", "|---|---|---|---|---|---|"]
    okk, comp = True, True
    for g in ("pur", "red", "kelp"):
        for c in ("rw", "pw"):
            r = C.get(f"v_time_{g}_{c}")
            if r is None:
                comp = False; body.append(f"| {NAME[g]} | {c} | not finished | | | |"); continue
            has_m = f"{g}:E10_mass" in r
            okk &= not excl(r[f"{g}:E10_biomass"] if has_m else r[f"{g}:E10_total"]) and (not has_m or not excl(r[f"{g}:E10_mass"]))
            body.append(f"| {NAME[g]} | {c} | {iv(r[f'{g}:E10_total'])} | {iv(r[f'{g}:E10_mass']) if has_m else 'no size data in window'} | {iv(r[f'{g}:E10_biomass']) if has_m else '—'} | {conv(r)} |")
    body += ["", verdict(status, "M5 time placebo", okk, comp), ""]

    # M6
    body += ["## M6. Placebo in space", "", "Protected stations dropped; one reference station per MPA relabelled as protected; ten random relabellings.", "",
             "| Taxon | Runs | Density: mean, excluding zero | Mean mass: mean, excluding zero | Biomass: mean, excluding zero |", "|---|---|---|---|---|"]
    okk, comp = True, True
    for g in ORDER:
        rs = [r for k, r in B.items() if k.startswith(f"v_space_{g}_")]
        if len(rs) < 10:
            comp = False
        if not rs:
            body.append(f"| {NAME[g]} | 0 | | | |"); continue
        row = [NAME[g], str(len(rs))]
        for q, _ in QS:
            if f"{g}:{q}" not in rs[0]:
                row.append("—"); continue
            m = np.array([r[f"{g}:{q}"]["mean"] for r in rs]); ex = np.array([excl(r[f"{g}:{q}"]) for r in rs])
            bad = ex.mean() > 0.25 or (abs(m.mean()) > 0.22 and (np.all(m > 0) or np.all(m < 0)))
            if q != "E10_total":
                okk &= not bad
            row.append(f"{m.mean():+.2f}, {ex.sum()}/{len(rs)}{' ✗' if bad else ''}")
        body.append("| " + " | ".join(row) + " |")
    body += ["", "✗ marks a failure of the rule. The rule is applied to mass and biomass; the density column repeats yesterday's density placebo under today's fits.", "",
             verdict(status, "M6 space placebo (mass, biomass)", okk, comp), ""]

    # M7
    body += ["## M7. Predictive check of the mass part", "",
             "Mean sizes replicated from the fitted model, compared with the data: spread and tails of group means by program, and two "
             "protected-versus-reference contrasts.", "", "| Taxon | Statistics outside their 90% replicate interval | Contrast statistics |", "|---|---|---|"]
    tot, out, cbad, comp = 0, 0, 0, True
    for g in ORDER:
        r = C.get(f"v_ppc_{g}")
        if r is None or "mass_ppc" not in r:
            comp = False; body.append(f"| {NAME[g]} | not finished | |"); continue
        fl = []; cs = []
        for k, v in r["mass_ppc"][g].items():
            if v["observed"] != v["observed"]:
                continue
            f = v["p_rep_ge_obs"] < 0.05 or v["p_rep_ge_obs"] > 0.95; tot += 1; out += f
            if f:
                fl.append(k)
            if k.startswith("contrast"):
                cbad += f; cs.append(f"{k.replace('contrast_', '')}: observed {v['observed']:+.2f}, replicates [{v['rep_q05']:+.2f}, {v['rep_q95']:+.2f}]{' ✗' if f else ''}")
        body.append(f"| {NAME[g]} | {len(fl)}: {', '.join(fl) or 'none'} | {'; '.join(cs)} |")
    body += ["", f"{out} of {tot} statistics outside (rule: at most 10%, and no contrast statistic outside).", "",
             verdict(status, "M7 predictive check of mass", tot > 0 and out <= 0.10 * tot and cbad == 0, comp), ""]

    # M8
    body += ["## M8. Predicting size data the model was not shown", "", "| Taxon | Split | Held-out groups | 90% intervals cover | 50% intervals cover | By program |", "|---|---|---|---|---|---|"]
    okk, comp = True, True
    for g in ORDER:
        for k in (1, 2):
            r = C.get(f"v_hold_{g}_s{k}")
            if r is None or g not in r.get("mass_holdout", {}):
                comp = False; body.append(f"| {NAME[g]} | {k} | not finished | | | |"); continue
            h = r["mass_holdout"][g]; okk &= 0.85 <= h["coverage90"] <= 0.95
            body.append(f"| {NAME[g]} | {k} | {h['n_groups']:,} | {h['coverage90']:.2f} | {h['coverage50']:.2f} | " + "; ".join(f"{p} {v['coverage90']:.2f} (n={v['n']})" for p, v in h["by_program"].items()) + " |")
    body += ["", verdict(status, "M8 hold-out prediction", okk, comp), ""]

    # M9
    body += ["## M9. Leave one MPA out", "", "| Taxon | Runs | Largest shift: density | mean mass | biomass | Sign reversals of clear main effects | Shifts above 0.25 |", "|---|---|---|---|---|---|---|"]
    okk, comp = True, True
    for g in ORDER:
        m = main_fit(F, g); rs = {k: r for k, r in B.items() if k.startswith(f"v_loo_{g}_")}
        if m is None or not rs:
            comp = False; body.append(f"| {NAME[g]} | {len(rs)} | | | | | |"); continue
        if len(rs) < len(m["mpas"][g]):
            comp = False
        row = [NAME[g], f"{len(rs)} of {len(m['mpas'][g])}"]; rev, big = [], []
        for q, lab in QS:
            sh = {r["omit_mpa"][0]: r[f"{g}:{q}"]["mean"] - m[f"{g}:{q}"]["mean"] for r in rs.values()}
            w = max(sh, key=lambda x: abs(sh[x])); row.append(f"{sh[w]:+.2f} ({w})")
            big += [f"{lab}: {k} {v:+.2f}" for k, v in sh.items() if abs(v) > 0.25]
            if excl(m[f"{g}:{q}"]):
                rev += [f"{lab}: {r['omit_mpa'][0]}" for r in rs.values() if np.sign(r[f"{g}:{q}"]["mean"]) != np.sign(m[f"{g}:{q}"]["mean"])]
        okk &= not rev
        body.append("| " + " | ".join(row + ["; ".join(rev) or "none", "; ".join(big) or "none"]) + " |")
    body += ["", verdict(status, "M9 leave one MPA out", okk, comp), ""]

    # M10, M11
    t = (V / "test_multi_real.log")
    body += ["## M10. Independent arithmetic", ""]
    if t.exists():
        txt = t.read_text(); p = "ALL CHECKS PASSED (1-3)" in txt
        body += ["`test_multi.py` on the real inputs:", "", "```", *[x for x in txt.strip().splitlines() if x[:2] in ("1.", "2.", "3.") or "PASSED" in x][-14:], "```", "", verdict(status, "M10 independent arithmetic", p), ""]
    else:
        body += [verdict(status, "M10 independent arithmetic", False, False), ""]
    pc = V / "prep_mass_checks.json"
    body += ["## M11. Input integrity", ""]
    if pc.exists():
        c = json.loads(pc.read_text())
        p = c["pools"]["unmatched"] == 0 and c["pools"]["max_rel_diff_mean"] < 1e-9 and all(v["max_rel_diff_mean"] < 1e-8 for v in c["same_record_units"].values()) and c["kfm_sheephead"]["unmatched"] == 0
        body += [f"{c['pools']['pools']:,} pools, {sum(v['units'] for v in c['same_record_units'].values()):,} same-record units and {c['kfm_sheephead']['review_exact_visits']} KFM sheephead visits "
                 "reproduce the integrated review's values from source rows.", "", verdict(status, "M11 input integrity", p), ""]
    else:
        body += [verdict(status, "M11 input integrity", False, False), ""]

    L += ["## Scorecard", "", "| Check | Result |", "|---|---|"] + [f"| {k} | {v} |" for k, v in status.items()] + [""]
    L += ["Known failures of the density part from 2026-11-10 that this battery does not re-test: sheephead time placebo; urchin count "
          "layer misfit; calibration narrowly below its rule; effect size depends on the curve prior.", ""]
    L += body
    (V / "check_status_mass_validation.json").write_text(json.dumps(status, indent=1))
    (ROOT / "REPORT_mass_validation.md").write_text("\n".join(L) + "\n")
    print("\n".join(f"{v:10s} {k}" for k, v in status.items()))


if __name__ == "__main__":
    main()
