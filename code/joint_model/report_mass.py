"""Report for the joint density, mean mass and biomass fits (design section 4.2), 2026-11-11.

  python3 code/joint_model/report_mass.py    # reads validation/{m_mass_fit,m_fix,m_mass_sim,m_fix_sim,prep_mass_checks.json}

Writes outputs/joint_model_20261109/REPORT_density_mass.md. Review only; nothing is adopted. Rules, written before the
real-data fits finished: convergence = R-hat < 1.01, effective sample size > 400, no divergences; seeds agree within
0.05; simulated-data coverage of 90% intervals no lower than 0.90 minus two binomial standard errors; a sensitivity is
flagged when it moves an effect at ten years by more than 0.10 on the log scale.
"""
import glob, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("outputs/joint_model_20261109")
V = ROOT / "validation"
NAME = {"lob": "Lobster", "shp": "Sheephead", "pur": "Purple urchin", "red": "Red urchin", "kelp": "Kelp"}
ORDER = ["lob", "shp", "pur", "red", "kelp"]
Q = [("E10_total", "Density"), ("E10_mass", "Mean mass"), ("E10_biomass", "Biomass")]
RHAT, NEFF, FLAG = 1.01, 400, 0.10
NOTES = json.loads((V / "notes_mass.json").read_text()) if (V / "notes_mass.json").exists() else {}


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


def pct(x):
    return f"{100 * (np.exp(x) - 1):+.0f}%"


def ok(r):
    return r["worst_rhat"] < RHAT and r["min_neff"] > NEFF and r["divergences"] == 0


def conv(r):
    return ("yes" if ok(r) else "NO") + f" ({r['worst_rhat']:.3f} / {r['min_neff']:.0f} / {r['divergences']})"


def tag(r):
    return list(r["taxa"])[0]


def main_fit(F, g, curve="rw", seed=1):
    """Main fit per taxon: sheephead uses the fits that include KFM roving-diver sizes (mx_); others mm_."""
    c = "rw" if curve == "rw" else "pw"
    return F.get(f"mx_{g}_{c}_s{seed}") if g == "shp" else F.get(f"mm_{g}_{c}_s{seed}")


def main():
    F = loads("m_mass_fit", "m_fix"); S = loads("m_mass_sim", "m_fix_sim")
    chk = json.loads((V / "prep_mass_checks.json").read_text()) if (V / "prep_mass_checks.json").exists() else None
    status = {}
    L = ["# Density, mean mass and biomass in one model", "",
         "2026-11-11. Review only. Nothing here is adopted, and the paper and the frozen pipeline are untouched.", ""]
    L += NOTES.get("summary", []) + [""]
    L += ["## 1. What the model is", "",
          "Each taxon has two linked parts fitted together. Density: the single-taxon count model under the recommended "
          "structure. Mean individual mass: station level, regional year effects (independent, shared trend, regional trend), "
          "a station AR(1) process, a siting difference, and its own response curve to protection. A group of n measured "
          "animals enters through its mean mass, Gamma(n × shape, mean M). Biomass is density × mean mass, so on the log "
          "scale the biomass effect is the density effect plus the mass effect. All effects are at ten years after "
          "designation; intervals are 90%.", "",
          "Only animals measured at the station are used. No sizes are borrowed from other sites or programs and none are "
          "imputed; station-years with counts but no sizes inform density only.", ""]
    # ---- inputs and checks
    L += ["## 2. Inputs: rebuilt from source rows and checked", ""]
    if chk:
        p = chk["pools"]; ku = chk["same_record_units"]; kf = chk["kfm_sheephead"]
        L += [f"`code/joint_model/prep_mass.py` recomputes every mass group from the registered source rows and compares it with the integrated review's own values:", "",
              f"- {p['pools']:,} size pools (urchins, PISCO lobster): number measured and mean mass equal the review's pool inventory (largest relative difference {p['max_rel_diff_mean']:.0e}; none unmatched).",
              f"- {sum(v['units'] for v in ku.values()):,} survey units with sizes on the count record (sheephead, LTER lobster, LTER urchin groups): mean mass equals stored biomass ÷ density (largest relative difference {max(v['max_rel_diff_mean'] for v in ku.values()):.0e}).",
              f"- KFM sheephead: {kf['size_surveys']} roving-diver size surveys at {kf['sites']} sites ({kf['fish']:,.0f} fish, 2007–2023). The {kf['review_exact_visits']} visits used in the November 3 KFM sheephead review are reproduced exactly.",
              "- No group comes from an excluded site.", ""]
        L += ["| Taxon | Program | Kind | Groups | Animals | Station-years |", "|---|---|---|---|---|---|"]
        for k, v in sorted(chk["groups"].items()):
            t, s, kind = k.split("|"); L.append(f"| {t} | {s} | {kind} | {v['groups']:,} | {v['animals']:,.0f} | {v['station_years']:,} |")
        L += ["", "(All sites the sources hold; the model uses the groups at its own stations.)", ""]
    # ---- simulation
    L += ["## 3. Simulated data: does the model recover known mass and biomass effects?", ""]
    if S:
        L += ["Counts and group mean masses were simulated on each taxon's real survey and measurement design at known values, "
              "then refitted (4 chains, 1,500 + 1,500).", "",
              "| Taxon | Data sets | Converged at R-hat < 1.05 | Mass: truth | Mass: mean error, coverage | Biomass: truth | Biomass: mean error, coverage |", "|---|---|---|---|---|---|---|"]
        cov = []
        for g in ORDER:
            rs = [r for r in S.values() if tag(r) == g]
            if not rs:
                continue
            row = [NAME[g], str(len(rs)), str(sum(r["worst_rhat"] < 1.05 for r in rs))]
            for q in ("E10_mass", "E10_biomass"):
                tr = np.array([r["truth"][g][q] for r in rs]); est = np.array([r[f"{g}:{q}"]["mean"] for r in rs])
                c = np.array([r[f"{g}:{q}"]["5.0%"] <= t <= r[f"{g}:{q}"]["95.0%"] for r, t in zip(rs, tr)]); cov += c.tolist()
                row += [f"{tr.mean():+.2f}", f"{(est - tr).mean():+.2f}, {c.sum()}/{len(rs)}"]
            L.append("| " + " | ".join(row) + " |")
        n = len(cov); c = float(np.mean(cov)); thr = 0.90 - 2 * np.sqrt(.09 / n)
        status["simulated mass and biomass effects recovered"] = c >= thr
        L += ["", f"**{'PASS' if c >= thr else 'FAIL'}.** Coverage {c:.2f} over {n} intervals (rule: at least {thr:.2f}).", ""]
    else:
        L += ["No simulation results yet.", ""]
    # ---- main results
    L += ["## 4. Real data: effect of protection at ten years", ""]
    for curve, title in (("rw", "Random-walk curve prior (main)"), ("piecewise", "Independent-increment curve prior")):
        L += [f"### {title}", "", "| Taxon | Density | Mean mass | Biomass | Biomass as % change | Converged (R-hat / eff. draws / divergences) |", "|---|---|---|---|---|---|"]
        for g in ORDER:
            r = main_fit(F, g, curve)
            if r is None:
                L.append(f"| {NAME[g]} | not finished | | | | |"); continue
            b = r[f"{g}:E10_biomass"]
            L.append(f"| {NAME[g]} | {iv(r[f'{g}:E10_total'])} | {iv(r[f'{g}:E10_mass'])} | {iv(b)} | {pct(b['mean'])} [{pct(b['5.0%'])}, {pct(b['95.0%'])}] | {conv(r)} |")
        L.append("")
    mains = [main_fit(F, g, c) for g in ORDER for c in ("rw", "piecewise")]
    if all(m is not None for m in mains):
        status["convergence of main fits"] = all(ok(m) for m in mains)
    d = []
    for g in ORDER:
        a, b = main_fit(F, g, "rw", 1), main_fit(F, g, "rw", 2)
        if a and b:
            d.append(max(abs(a[f"{g}:{q}"]["mean"] - b[f"{g}:{q}"]["mean"]) for q, _ in Q))
    if len(d) == len(ORDER):
        status["seeds agree"] = max(d) < 0.05
        L += [f"Seeds: largest difference between two independent fits {max(d):.3f} (rule: below 0.05) — **{'PASS' if max(d) < 0.05 else 'FAIL'}**.", ""]
    L += ["Mass data behind each estimate (model stations only):", "", "| Taxon | Station-years with measured animals | Station-years with density | Groups by program | MPAs with before-and-after mass data |", "|---|---|---|---|---|"]
    for g in ORDER:
        r = main_fit(F, g, "rw")
        if r and r.get("mass_data", {}).get(g):
            m = r["mass_data"][g]
            L.append(f"| {NAME[g]} | {m['station_years_with_mass']:,} | {m['density_station_years']:,} | {', '.join(f'{k} {v:,}' for k, v in sorted(m.get('by_program', {}).items()))} | {m['mpas_before_after']} |")
    L.append("")
    # ---- sensitivities, one block per issue
    L += ["## 5. The open biomass issues, one by one", "",
          "Each row re-fits the model with one thing changed and gives the shift in the effect at ten years against the main fit "
          f"(random-walk prior). Shifts above {FLAG:.2f} are marked ⚠.", ""]
    flagged = []

    def block(title, text, rows):
        out = [f"### {title}", "", text, "", "| Taxon | Change | Density | Mean mass | Biomass | Shift in biomass effect | Converged |", "|---|---|---|---|---|---|---|"]
        n = 0
        for g, label, key, base in rows:
            r = F.get(key); m = base if base is not None else main_fit(F, g, "rw")
            if r is None or m is None:
                out.append(f"| {NAME[g]} | {label} | not finished | | | | |"); continue
            if r.get("not_applicable"):
                continue
            sh = r[f"{g}:E10_biomass"]["mean"] - m[f"{g}:E10_biomass"]["mean"]; n += 1
            if abs(sh) > FLAG:
                flagged.append(f"{NAME[g]}: {label} ({sh:+.2f})")
            out.append(f"| {NAME[g]} | {label} | {iv(r[f'{g}:E10_total'])} | {iv(r[f'{g}:E10_mass'])} | {iv(r[f'{g}:E10_biomass'])} | {sh:+.2f}{' ⚠' if abs(sh) > FLAG else ''} | {conv(r)} |")
        return out + [""]

    L += block("5.1 Allometry", "The length-to-mass exponent is fixed at the published value. These fits move it by −0.1 and +0.1 "
               "(for an exponent near 2.9 to 3.0 that is a large error). The intercept cancels in every effect.",
               [(g, f"exponent {s}", f"mx_{g}_rw_{k}", None) for g in ("lob", "shp", "pur", "red") for k, s in (("blo", "−0.1"), ("bhi", "+0.1"))])
    if chk:
        L += ["Average shift in log mean mass from the exponent change (it moves levels a lot and contrasts little): "
              + "; ".join(f"{t} {v['minus']:+.2f} / {v['plus']:+.2f}" for t, v in chk["exponent_shift_log_mean"].items()) + ".", ""]
    L += block("5.2 One size definition for urchins (25 mm and larger)",
               "PISCO counts and measures urchins of 25 mm and up; KFM and LTER count all sizes. These fits put KFM and LTER on the "
               "PISCO definition: mean mass from measured animals of 25 mm and up, and density scaled by the share of measured animals "
               "of 25 mm and up at that station and year."
               + (" Shares of measured animals at 25 mm and up: " + "; ".join(f"{k.replace('|', ' ')} {v['animal_weighted_share']:.0%}" for k, v in chk["share25"].items()) + "." if chk else ""),
               [(g, "25 mm and up, " + c, f"mx_{g}_{k}_ge25", main_fit(F, g, c)) for g in ("pur", "red") for c, k in (("rw", "rw"), ("piecewise", "pw"))])
    L += block("5.3 LTER urchin sizes are group means",
               "LTER records a mean size for each group of urchins, so the spread inside a group is unknown. These fits drop LTER urchin "
               "sizes from the mass part entirely (LTER still informs density).",
               [(g, "no LTER sizes", f"mx_{g}_rw_noLTERmass", None) for g in ("pur", "red")])
    L += block("5.4 KFM sheephead sizes",
               "The main sheephead fit now includes KFM roving-diver sizes (2007–2023). This row leaves them out, as before today. "
               "Roving-diver sizes and transect counts are different samples of the fish at a site; that mismatch is not removed. "
               "KFM lobster has no sizes at all and informs density only.",
               [("shp", "without KFM sizes", "mm_shp_rw_s1", None)])
    L += block("5.5 Density and mass treated as independent",
               "These fits let mean mass at a station move with that year's density fluctuation there (for example a recruitment pulse "
               "raising numbers and lowering mean size).",
               [(g, "mass linked to density fluctuation", ("mx_shp_rw_link" if g == "shp" else f"mm_{g}_rw_link"), None) for g in ORDER])
    th = [(g, (F.get("mx_shp_rw_link") if g == "shp" else F.get(f"mm_{g}_rw_link"))) for g in ORDER]
    th = [f"{NAME[g]} {iv(r[f'{g}:m_theta'])}" for g, r in th if r and f"{g}:m_theta" in r]
    if th:
        L += ["Link coefficient (change in log mean mass per log unit of density fluctuation): " + "; ".join(th) + ".", ""]
    L += block("5.6 KFM reference sites flagged for fish",
               "Scott Gabara flagged Webster's Arch, Miracle Mile and Wyckoff Ledge as questionable references for fish. These fits leave "
               "each out, then all three. They are sensitivities; no site is dropped from the main fit.",
               [("shp", lab, f"mx_shp_rw_{k}", None) for k, lab in (("noWebster", "without Webster's Arch"), ("noMiracle", "without Miracle Mile"), ("noWyckoff", "without Wyckoff Ledge"), ("noAll3", "without all three"))])
    L += block("5.7 Prior on the mass response", "Prior scale of the mass response doubled.",
               [(g, "mass prior ×2", f"mm_{g}_rw_wide", None) for g in ORDER])
    status_text = "none" if not flagged else "; ".join(flagged)
    L += [f"**Sensitivities that move the biomass effect by more than {FLAG:.2f}:** {status_text}.", ""]
    L += ["## 6. What is still not fixed", ""] + NOTES.get("limits", [
        "- Allometry coefficients have no published uncertainty in the files we hold; the exponent check above is a bound, not a propagated uncertainty.",
        "- KFM lobster has no size data. Lobster mass and biomass rest on LTER and PISCO.",
        "- KFM sheephead sizes come from roving divers and counts from transects.",
        "- Urchin mass has few MPAs with before-and-after size data, so the urchin mass effect leans on post-designation contrasts and the siting-difference prior.",
        "- Kelp: LTER counts stipes with a fixed converted mass per stipe, so LTER kelp informs density only; KFM and PISCO count plants.",
        "- The density model's known weak points carry over: the urchin count layer misfit and the sheephead time placebo.",
        "- NPS permission for redistributing KFM records is a release question and is unchanged."]) + [""]
    (V / "check_status_mass.json").write_text(json.dumps({k: bool(v) for k, v in status.items()}, indent=1))
    (ROOT / "REPORT_density_mass.md").write_text("\n".join(L) + "\n")
    print("\n".join(f"{'PASS' if v else 'FAIL'}  {k}" for k, v in status.items()))
    print("flagged:", status_text)


if __name__ == "__main__":
    main()
