"""Report for the multi-taxon model (stage 2), 2026-11-11.

  python3 code/joint_model/report_m.py     # reads outputs/joint_model_20261109/validation/{m_sim,m_fit}

Writes outputs/joint_model_20261109/REPORT_multi_taxon.md. Review only; nothing is adopted. The pass rules below were
written before any result existed:
  simulated data, 90% intervals: coverage of the truth no lower than 0.90 minus two binomial standard errors;
      share of intervals that exclude zero when the truth is zero no higher than 0.10 plus two standard errors;
      mean error of each path coefficient below 0.05 in absolute value.
  real data: R-hat < 1.01, effective sample size > 400, no divergences (as in the single-taxon reports);
      seeds agree within 0.05; the no-links fit reproduces the single-taxon estimates within 0.05.
"""
import glob, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("outputs/joint_model_20261109")
V = ROOT / "validation"
NAME = {"lob": "Lobster", "shp": "Sheephead", "pur": "Purple urchin", "red": "Red urchin", "kelp": "Kelp"}
ORDER = ["lob", "shp", "pur", "red", "kelp"]
# single-taxon estimates under the same structure (REPORT_recommended_structure.md and section 8 of the structure report)
SINGLE = {"rw": {"pur": -0.95, "red": -0.63, "lob": 1.74, "shp": 0.21, "kelp": 0.93},
          "piecewise": {"pur": -0.30, "red": -0.28, "lob": 1.82, "shp": 0.15}}
RHAT, NEFF = 1.01, 400
NOTES = json.loads((V / "notes_m.json").read_text()) if (V / "notes_m.json").exists() else {}


def loads(gate, pattern):
    out = {}
    for p in sorted(glob.glob(str(V / gate / f"{pattern}.json"))):
        try:
            out[Path(p).stem] = json.loads(Path(p).read_text())
        except Exception:
            pass
    return out


def iv(d, i=None):
    m, lo, hi = d["mean"], d["5.0%"], d["95.0%"]
    if i is not None:
        m, lo, hi = m[i], lo[i], hi[i]
    return f"{m:+.2f} [{lo:+.2f}, {hi:+.2f}]"


def excl(d, i=None):
    lo, hi = (d["5.0%"], d["95.0%"]) if i is None else (d["5.0%"][i], d["95.0%"][i])
    return lo > 0 or hi < 0


def conv(r):
    ok = r["worst_rhat"] < RHAT and r["min_neff"] > NEFF and r["divergences"] == 0
    return ("yes" if ok else "NO") + f" ({r['worst_rhat']:.3f} / {r['min_neff']:.0f} / {r['divergences']})"


def edge_label(e):
    k, j = e.split("<-")
    return f"{NAME[j]} → {NAME[k]}"


def two_se(p, n):
    return 2 * np.sqrt(p * (1 - p) / max(n, 1))


def sim_section(L):
    sims = loads("m_sim", "ms_*")
    L += ["## 2. Simulated data: does the model recover known links?", ""]
    if not sims:
        L += ["No simulation results yet.", ""]; return {}
    status = {}
    L += ["Counts were simulated on the real survey designs of all five taxa from the linked model, with known link "
          "strengths, then fitted with the same code as the real data (4 chains, 1,500 warm-up and 1,500 kept draws; "
          "shorter than the real-data fits). Intervals are 90%.", ""]
    for sc, title in (("paths", "Links present"), ("nopaths", "No links"), ("confounded", "No links, but a local driver shared by all taxa")):
        rs = {k: r for k, r in sims.items() if r["scenario"] == sc}
        if not rs:
            continue
        edges = next(iter(rs.values()))["edges"]
        n = len(rs); nconv = sum(r["worst_rhat"] < 1.05 for r in rs.values())
        L += [f"### {title} ({n} simulated data sets; {nconv} with R-hat below 1.05)", "",
              "| Link | Truth | Mean estimate | Mean error | Interval covers truth | Interval excludes zero |", "|---|---|---|---|---|---|"]
        cov_all, exc_all, bias_bad = [], [], []
        for i, e in enumerate(edges):
            tr = np.array([r["truth"]["beta"][i] for r in rs.values()])
            est = np.array([r["beta"]["mean"][i] for r in rs.values()])
            lo = np.array([r["beta"]["5.0%"][i] for r in rs.values()]); hi = np.array([r["beta"]["95.0%"][i] for r in rs.values()])
            cov = (lo <= tr) & (tr <= hi); ex = (lo > 0) | (hi < 0)
            cov_all += cov.tolist(); exc_all += ex.tolist()
            if abs((est - tr).mean()) >= 0.05:
                bias_bad.append(edge_label(e))
            L.append(f"| {edge_label(e)} | {tr[0]:+.2f} | {est.mean():+.2f} | {(est - tr).mean():+.3f} | {cov.sum()}/{n} | {ex.sum()}/{n} |")
        L.append("")
        # effects at ten years
        L += ["| Taxon | True direct | Direct: mean error, coverage | True total | Total: mean error, coverage |", "|---|---|---|---|---|"]
        ecov = []
        for g in ORDER:
            row = [NAME[g]]
            for q in ("E10_direct", "E10_total"):
                tr = np.array([r["truth"][g][q] for r in rs.values()])
                est = np.array([r[f"{g}:{q}"]["mean"] for r in rs.values()])
                lo = np.array([r[f"{g}:{q}"]["5.0%"] for r in rs.values()]); hi = np.array([r[f"{g}:{q}"]["95.0%"] for r in rs.values()])
                c = (lo <= tr) & (tr <= hi); ecov += c.tolist()
                row += [f"{tr.mean():+.2f}", f"{(est - tr).mean():+.2f}, {c.sum()}/{n}"]
            L.append("| " + " | ".join(row) + " |")
        L.append("")
        nb = len(cov_all); c = float(np.mean(cov_all)); x = float(np.mean(exc_all)); ec = float(np.mean(ecov))
        if sc == "paths":
            ok = c >= 0.90 - two_se(0.9, nb) and not bias_bad and ec >= 0.90 - two_se(0.9, len(ecov))
            status["links recovered"] = ok
            L += [f"**{'PASS' if ok else 'FAIL'}.** Link coverage {c:.2f} over {nb} intervals (rule: at least {0.90 - two_se(0.9, nb):.2f}); "
                  f"effect coverage {ec:.2f} over {len(ecov)} (rule: at least {0.90 - two_se(0.9, len(ecov)):.2f}); "
                  + ("no link with mean error of 0.05 or more." if not bias_bad else "mean error of 0.05 or more for: " + ", ".join(bias_bad) + "."),
                  f" With true links present, {100 * x:.0f}% of link intervals excluded zero (the power at these link strengths).", ""]
        elif sc == "nopaths":
            ok = x <= 0.10 + two_se(0.1, nb) and ec >= 0.90 - two_se(0.9, len(ecov))
            status["no false links"] = ok
            L += [f"**{'PASS' if ok else 'FAIL'}.** {100 * x:.0f}% of {nb} link intervals excluded zero when there was no link "
                  f"(rule: at most {100 * (0.10 + two_se(0.1, nb)):.0f}%); effect coverage {ec:.2f} (rule: at least {0.90 - two_se(0.9, len(ecov)):.2f}).", ""]
        else:
            L += [f"**Not a pass/fail check.** The simulated driver (SD {next(iter(rs.values()))['truth']['common_sd']}, persistent, same sign for "
                  f"every taxon) is not in the fitted model. {100 * x:.0f}% of {nb} link intervals excluded zero although no link "
                  "existed. This is the size of error an unmeasured local driver can cause.", ""]
    return status


def real_section(L):
    fits = loads("m_fit", "m_*")
    L += ["## 3. Real data", ""]
    if not fits:
        L += ["No real-data fits yet.", ""]; return {}
    status = {}
    main = fits.get("m_rw_s1") or fits.get("m_quick_rw")
    which = "m_rw_s1" if "m_rw_s1" in fits else "m_quick_rw"
    L += ["### 3.1 Convergence", "", "| Fit | What it is | Converged (R-hat / effective draws / divergences) | Hours |", "|---|---|---|---|"]
    what = {"m_quick_rw": "short first look (1,000 + 1,500)", "m_rw_s1": "main, random-walk curve prior, seed 1", "m_rw_s2": "same, seed 2",
            "m_pw_s1": "independent-increment curve prior", "m_rw_nopaths": "no links (taxa side by side)", "m_pw_nopaths": "no links, independent-increment prior",
            "m_rw_lag0": "links within the same year", "m_rw_lag2": "links two years back", "m_rw_beta1": "link prior twice as wide",
            "m_rw_noPD": "without the two Point Dume MPAs", "m_rw_reverse": "links reversed (diagnostic)"}
    allok = True
    for k in what:
        if k in fits:
            r = fits[k]; c = conv(r); allok &= c.startswith("yes") or k == "m_quick_rw"
            L.append(f"| `{k}` | {what[k]} | {c} | {r['seconds'] / 3600:.1f} |")
    L.append("")
    status["convergence"] = allok
    if main is None:
        return status
    edges = main["edges"]
    L += [f"### 3.2 Links (from `{which}`)", "",
          "A link is an elasticity: the change in the log density of the second taxon that follows a change of one log unit "
          "in the first taxon at the same station a year earlier. A value of −0.2 means a doubling is followed by a "
          "13% drop. Intervals are 90%.", "",
          "| Link | Estimate | Chance it is negative | Station-years linked |", "|---|---|---|---|"]
    for i, e in enumerate(edges):
        pn = main["intervals95"]["beta"]["prob_negative"][i]
        lk = main["links"][e]
        L.append(f"| {edge_label(e)} | {iv(main['beta'], i)} | {100 * pn:.0f}% | {lk['n_linked']} of {lk['n_target_station_years']} |")
    L.append("")
    L += [f"### 3.3 Effect of protection at ten years (log scale; from `{which}`)", "",
          "Direct = the part not carried by the measured links. Total = direct plus what the links carry.", "",
          "| Taxon | Direct | Carried by links | Total | Single-taxon estimate |", "|---|---|---|---|---|"]
    for g in ORDER:
        via = [f"{NAME[k.split('_via_')[1]].lower()} {main[k]['mean']:+.2f}" for k in main if k.startswith(f"{g}:E10_via_")]
        L.append(f"| {NAME[g]} | {iv(main[f'{g}:E10_direct'])} | {'; '.join(via) or '—'} | {iv(main[f'{g}:E10_total'])} | {SINGLE['rw'].get(g, float('nan')):+.2f} |")
    L.append("")
    # agreement and sensitivities
    L += ["### 3.4 Do the answers hold?", "", "Links (mean estimate):", "", "| Fit | " + " | ".join(edge_label(e) for e in edges) + " |", "|---|" + "---|" * len(edges)]
    for k in ("m_rw_s1", "m_rw_s2", "m_pw_s1", "m_rw_lag0", "m_rw_lag2", "m_rw_beta1", "m_rw_noPD", "m_quick_rw"):
        if k in fits and fits[k].get("edges") == edges:
            L.append(f"| `{k}` | " + " | ".join(f"{x:+.2f}{'*' if excl(fits[k]['beta'], i) else ''}" for i, x in enumerate(fits[k]["beta"]["mean"])) + " |")
    L += ["", "\\* 90% interval excludes zero.", "", "Total effect at ten years (mean):", "", "| Fit | " + " | ".join(NAME[g] for g in ORDER) + " |", "|---|" + "---|" * len(ORDER)]
    for k in ("m_rw_s1", "m_rw_s2", "m_rw_nopaths", "m_pw_s1", "m_pw_nopaths", "m_rw_lag0", "m_rw_lag2", "m_rw_beta1", "m_rw_noPD", "m_quick_rw"):
        if k in fits:
            L.append(f"| `{k}` | " + " | ".join(f"{fits[k][f'{g}:E10_total']['mean']:+.2f}{'*' if excl(fits[k][f'{g}:E10_total']) else ''}" for g in ORDER) + " |")
    L.append("| single-taxon, random-walk prior | " + " | ".join(f"{SINGLE['rw'][g]:+.2f}" for g in ORDER) + " |")
    L.append("| single-taxon, independent-increment prior | " + " | ".join(f"{SINGLE['piecewise'][g]:+.2f}" if g in SINGLE["piecewise"] else "—" for g in ORDER) + " |")
    L.append("")
    if "m_rw_s1" in fits and "m_rw_s2" in fits:
        a, b = fits["m_rw_s1"], fits["m_rw_s2"]
        d = max([abs(x - y) for x, y in zip(a["beta"]["mean"], b["beta"]["mean"])] + [abs(a[f"{g}:E10_total"]["mean"] - b[f"{g}:E10_total"]["mean"]) for g in ORDER])
        status["seeds agree"] = d < 0.05
        L += [f"Seeds: largest difference between the two main fits {d:.3f} (rule: below 0.05) — **{'PASS' if d < 0.05 else 'FAIL'}**.", ""]
    if "m_rw_nopaths" in fits:
        d = max(abs(fits["m_rw_nopaths"][f"{g}:E10_total"]["mean"] - SINGLE["rw"][g]) for g in ORDER)
        status["no-links fit reproduces single-taxon fits"] = d < 0.05
        L += [f"No-links fit against the single-taxon fits: largest difference {d:.3f} (rule: below 0.05) — **{'PASS' if d < 0.05 else 'FAIL'}**.", ""]
    if "m_rw_reverse" in fits:
        r = fits["m_rw_reverse"]
        L += ["### 3.5 Reversed links (diagnostic)", "",
              "Same model with every link turned around (kelp last year → urchins; urchins last year → predators). "
              "If a reversed link is as strong as the forward one and has the same sign, the data cannot tell the direction, "
              "and a shared local driver is as good an explanation as predation or grazing.", "",
              "| Reversed link | Estimate | Forward link | Estimate |", "|---|---|---|---|"]
        for i, e in enumerate(r["edges"]):
            k, j = e.split("<-"); f = f"{j}<-{k}"
            fw = iv(main["beta"], edges.index(f)) if f in edges else "—"
            L.append(f"| {edge_label(e)} | {iv(r['beta'], i)} | {edge_label(f)} | {fw} |")
        L.append("")
    return status


def main():
    L = ["# Multi-taxon model: all five taxa in one model (stage 2)", "",
         "2026-11-11. Review only. Nothing here is adopted, and the paper and the frozen pipeline are untouched.", ""]
    L += NOTES.get("summary", []) + [""]
    L += ["## 1. What the model is", "",
          "One model for lobster, sheephead, purple urchin, red urchin and kelp diver counts. Each taxon keeps the "
          "single-taxon model under the recommended structure (regional random-walk trends, three regions, wider prior on "
          "the first response increment). The taxa are tied together by six links at the station level: lobster and "
          "sheephead last year → each urchin this year; each urchin last year → kelp this year.", "",
          "The predictor in a link is the source taxon's latent log density at that station, relative to that station's "
          "own long-run level and with the source's regional year effect removed. So a link is learned from within-station "
          "change over time and from differences between MPAs in how strongly the source responded to protection. "
          "It is not learned from region-wide swings (a heatwave moving every taxon at once) or from one station having "
          "more of everything than another.", "",
          "Because the source's protection effect is part of the predictor, the model splits each downstream effect of "
          "protection into a direct part and a part carried by the links.", "",
          "Code checks (`code/joint_model/test_multi.py`): with the links switched off, each taxon's part of the model "
          "has the same log density as the single-taxon code to 1e-16, for both curve priors; every link points at the "
          "same station one year earlier.", ""]
    L += NOTES.get("after_model", [])
    status = sim_section(L)
    status.update(real_section(L))
    L += ["## 4. Limits", ""] + NOTES.get("limits", [
        "- The links are associations within stations over time. They are consistent with predation and grazing; they do not prove them.",
        "- The split into direct and carried parts assumes the elasticity seen in station-level change also applies to the change caused by protection.",
        "- Kelp is diver counts (plants for KFM and PISCO, stipes for LTER), not biomass or canopy.",
        "- The count layer misfits urchin counts (single-taxon check E); the multi-taxon model inherits that.",
        "- The single-taxon placebo failures (sheephead in time; kelp in space) are not addressed by linking taxa."]) + [""]
    (V / "check_status_multi.json").write_text(json.dumps({k: bool(v) for k, v in status.items()}, indent=1))
    (ROOT / "REPORT_multi_taxon.md").write_text("\n".join(L) + "\n")
    print("\n".join(f"{'PASS' if v else 'FAIL'}  {k}" for k, v in status.items()))


if __name__ == "__main__":
    main()
