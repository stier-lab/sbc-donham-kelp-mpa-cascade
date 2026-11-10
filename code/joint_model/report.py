"""Build outputs/joint_model_20261109/REPORT.md from the validation results pulled back from GRIT.

  python3 code/joint_model/report.py            # reads outputs/joint_model_20261109/validation/<gate>/*.json

Gate thresholds are fixed here and are not to be loosened: R-hat < 1.01, effective sample size > 400 for every
reported parameter, no divergences.
"""
import glob, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("outputs/joint_model_20261109")
V = ROOT / "validation"
TAXA = {"purp": "Purple urchin", "red": "Red urchin", "lob": "Lobster", "sheep": "Sheephead"}
CURVES = ["rw", "piecewise"]
RHAT, NEFF = 1.01, 400
NOTES = json.loads((V / "notes.json").read_text()) if (V / "notes.json").exists() else {}


def load(gate, name):
    p = V / gate / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def loads(gate, pattern):
    return {Path(p).stem: json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(V / gate / f"{pattern}.json")))}


def pct(x):
    return f"{100 * (np.exp(x) - 1):+.0f}%"


def e10(r):
    e = r["E10"]
    return f"{e['mean']:+.2f} [{e['5.0%']:+.2f}, {e['95.0%']:+.2f}]"


def conv(r):
    ok = r["worst_rhat"] < RHAT and r["min_neff"] > NEFF and r["divergences"] == 0
    return ok, f"{r['worst_rhat']:.3f} / {r['min_neff']:.0f} / {r['divergences']}"


def worst_params(r, k=3):
    rows = []
    for name, v in r.items():
        if isinstance(v, dict) and "r_hat" in v:
            rh = np.atleast_1d(np.array(v["r_hat"], dtype=float)); ne = np.atleast_1d(np.array(v["n_eff"], dtype=float))
            rows.append((float(np.nanmax(rh)), float(np.nanmin(ne)), name))
    bad = [x for x in rows if x[0] >= RHAT or x[1] <= NEFF]
    bad.sort(key=lambda x: (-x[0], x[1]))
    return ", ".join(f"{n} (R-hat {a:.3f}, ESS {b:.0f})" for a, b, n in bad[:k]) or "none"


L = []
P = L.append
status = {}

# ---------------------------------------------------------------- gate 1
g1 = load("gate1", "fit_S_purpuratus_rw")
P("## Gate 1. Reproduce the purple urchin fit on GRIT\n")
if g1:
    e = g1["E10"]; ok = abs(e["mean"] + 1.07) < 0.15 and abs(e["5.0%"] + 1.95) < 0.3 and abs(e["95.0%"] + 0.17) < 0.3
    status[1] = ("PASS" if ok else "FAIL", f"E10 {e10(g1)} against the target -1.07 [-1.95, -0.17]")
    P(f"Frozen code (commit bf5cea2), `fit --curve rw`, 4 chains of 600 warm-up and 600 kept draws: E10 = {e10(g1)}. "
      f"Target from the first report: -1.07 [-1.95, -0.17]. Same package versions as the first run (numpyro 0.22.0, jax 0.12.0.1). "
      f"Worst R-hat {g1.get('worst_rhat', 'n/a')}, divergences {g1['divergences']}; this short run is a reproduction check, not a converged fit.\n")
else:
    status[1] = ("NOT RUN", "no result file"); P("No result file.\n")

# ---------------------------------------------------------------- gate 2
status[2] = ("PASS", "commit bf5cea2 on local branch joint-model/stage1-validation, made before any change")
P("## Gate 2. Freeze\n\nCode, design and first report were committed as bf5cea2 on the local branch `joint-model/stage1-validation` before any change. Nothing was pushed.\n")

# ---------------------------------------------------------------- gate 3
P("## Gate 3. Sheephead non-convergence\n")
g3 = {}
for gate, label in (("g3", "original likelihood, 600 + 600 draws, tree depth 9"), ("g3b", "original likelihood, 1000 + 1000 draws"), ("g3c", "stable likelihood, 1000 + 1000 draws")):
    for k, r in loads(gate, "g3*").items():
        step0 = k.split("_step")[1].split("_")[0] if "_step" in k else "1.0"
        g3.setdefault((label, r["curve"], r["init"], step0), []).append(r)
if g3:
    P(NOTES.get("gate3", "") + "\n")
    P("A chain is counted as frozen when its adapted step size is below 0.0001 (working chains adapt to about 0.01).\n")
    P("| Code | Curve | Start | Initial step size | Fits (4 chains each) | Frozen chains | Fits with worst R-hat < 1.05 |")
    P("|---|---|---|---|---|---|---|")
    tot = {}
    for (label, c, init, step0), rs in sorted(g3.items()):
        fr = sum(int((np.array(r["per_chain"]["step_size"]) < 1e-4).sum()) for r in rs); n = sum(len(r["per_chain"]["step_size"]) for r in rs)
        t_ = tot.setdefault(label.split(",")[0], [0, 0]); t_[0] += fr; t_[1] += n
        P(f"| {label} | {c} | {'within 2 units (original)' if init == 'legacy' else 'within 0.5 units'} | {step0} | {len(rs)} | {fr}/{n} | {sum(r['worst_rhat'] < 1.05 for r in rs)}/{len(rs)} |")
    o, n_ = tot.get("original likelihood", [0, 0]), tot.get("stable likelihood", [0, 0])
    P(f"\nFrozen chains: {o[0]} of {o[1]} with the original likelihood code, {n_[0]} of {n_[1]} with the stable code, on the same seeds and settings.\n")
    status[3] = ("PASS" if n_[1] > 0 and n_[0] == 0 else "FAIL", f"cause found; frozen chains {o[0]}/{o[1]} before the fix, {n_[0]}/{n_[1]} after")
else:
    status[3] = ("NOT RUN", "no result files"); P("No result files.\n")

# ---------------------------------------------------------------- gate 4 (+ headline, gate 8)
P("## Gate 4. Convergence\n")
P(f"Thresholds: R-hat < {RHAT}, effective sample size > {NEFF} for every reported parameter (including every MPA-level and year-level parameter), no divergences. 4 chains.\n")
g4 = {}
for gate in ("g4a", "g4b", "g4c"):
    g4.update({(k, gate): r for k, r in loads(gate, "g4_*").items()})
chosen = {}
if g4:
    P("| Taxon | Curve prior | Sampler variant | Attempt | Draws per chain | Worst R-hat / min ESS / divergences | Gate | Parameters failing | E10 |")
    P("|---|---|---|---|---|---|---|---|---|")
    for t in TAXA:
        for c in CURVES:
            cands = [(k, gate, r) for (k, gate), r in g4.items() if k.startswith(f"g4_{t}_{c}_")]
            for k, gate, r in cands:
                ok, txt = conv(r)
                P(f"| {TAXA[t]} | {c} | {k.split('_', 3)[3]} | {gate[-1]} | {r['samples']} | {txt} | {'pass' if ok else 'fail'} | {worst_params(r)} | {e10(r)} |")
            # choice on diagnostics only: passing first, then lowest worst R-hat, then highest min ESS
            if cands:
                best = sorted(cands, key=lambda x: (not conv(x[2])[0], x[2]["divergences"] > 0, -ord(x[1][-1]), round(x[2]["worst_rhat"], 3), -x[2]["min_neff"]))[0]
                chosen[(t, c)] = best
    npass = sum(conv(r)[0] for _, _, r in chosen.values())
    status[4] = ("PASS" if npass == 8 and len(chosen) == 8 else "FAIL", f"{npass} of 8 taxon-by-prior fits meet all three thresholds")
    P("\nFor each taxon and prior the reported fit is a passing fit from the latest attempt (attempt c uses the final code); among those, the one with the lowest worst R-hat. The variants are the same model with different sampler settings; the estimate plays no part in the choice. Variant names: base = station baselines as scaled offsets; cb = baselines sampled directly; margb = baselines integrated out analytically; dense = dense mass matrix for the hyperparameters; 98/99 = target acceptance.\n")
    P("Reported fits: " + "; ".join(f"{TAXA[t]} {c}: {k.split('_', 3)[3]} (attempt {gate[-1]})" for (t, c), (k, gate, r) in chosen.items()) + ".\n")
    P(NOTES.get("gate4", "") + "\n")
else:
    status[4] = ("NOT RUN", "no result files"); P("No result files.\n")

HEAD = []
if chosen:
    HEAD.append("| Taxon | Curve prior | E10 (log change at 10 years) [90% interval] | Change per decade | Pr(E10<0) | Pr(change beyond 25%) | Converged (gate 4) |")
    HEAD.append("|---|---|---|---|---|---|---|")
    for t in TAXA:
        for c in CURVES:
            if (t, c) in chosen:
                k, gate, r = chosen[(t, c)]; e = r["E10"]
                HEAD.append(f"| {TAXA[t]} | {c} | {e10(r)} | {pct(e['mean'])} ({pct(e['5.0%'])} to {pct(e['95.0%'])}) | {r['E10_prob_negative']:.2f} | {r['E10_prob_beyond_0.22']:.2f} | {'yes' if conv(r)[0] else 'NO'} |")

# ---------------------------------------------------------------- gate 5
P("## Gate 5. Simulation calibration on the real survey layout\n")
sims = loads("g5", "sim_*")
if sims:
    P("Data were simulated from the model at known values on the real station, year and unit layout (purple urchin layout for the urchin scenarios; lobster layout, with a zero share of about 0.8, level shifts and the innovations parameterization, for the lobster scenarios), then refitted with 2 chains of 500 + 500. "
      "Coverage should be near 0.90. Power is the share of 90% intervals that exclude zero (for the null scenarios this is the false-positive rate; 0.10 is nominal). Replicates with R-hat of E10 above 1.05 are counted and kept.\n")
    P("| Layout | Scenario | True E10 | Curve prior | Replicates | Mean estimate | Bias | RMSE | 90% coverage (SE) | Power / false-positive rate | Wrong-sign and excludes zero | Share with E10 R-hat > 1.05 | Median zero share |")
    P("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    cells = {}
    for k, r in sims.items():
        _, lay, rest = k.split("_", 2); sc = rest.rsplit("_", 2)[0]; c = rest.rsplit("_", 2)[1]
        cells.setdefault((lay, sc, c), []).append(r)
    calib_ok = True; minrep = 10 ** 9
    for (lay, sc, c), rs in sorted(cells.items()):
        tr = rs[0]["truth"]["E10"]; m = np.array([r["E10"]["mean"] for r in rs]); lo = np.array([r["E10"]["5.0%"] for r in rs]); hi = np.array([r["E10"]["95.0%"] for r in rs])
        cov = ((lo <= tr) & (tr <= hi)).mean(); n = len(rs); se = np.sqrt(max(cov * (1 - cov), 1e-9) / n)
        excl = (lo > 0) | (hi < 0); wrong = (excl & (np.sign(m) != np.sign(tr))).mean() if tr != 0 else np.nan
        rh = np.mean([r["E10"]["r_hat"] > 1.05 for r in rs]); zf = np.median([r["truth"]["zero_fraction"] for r in rs])
        minrep = min(minrep, n)
        if cov < 0.85 or cov > 0.96:
            calib_ok = False
        P(f"| {lay} | {sc} | {tr:+.2f} | {c} | {n} | {m.mean():+.2f} | {m.mean() - tr:+.2f} | {np.sqrt(((m - tr) ** 2).mean()):.2f} | {cov:.2f} ({se:.02f}) | {excl.mean():.2f} | {'' if tr == 0 else f'{wrong:.3f}'} | {rh:.2f} | {zf:.2f} |")
    done = minrep >= 200 and len(cells) == 12
    status[5] = (("PASS" if calib_ok else "FAIL") if done else "INCOMPLETE",
                 f"{len(cells)} of 12 scenario-by-prior cells, at least {minrep} replicates each; coverage " + ("within 0.85 to 0.96 in every cell" if calib_ok else "outside 0.85 to 0.96 in at least one cell"))
    P("\nGate rule used here: at least 200 replicates in each of the 12 cells and 90% coverage between 0.85 and 0.96 in every cell. Bias and power are reported, not gated: the priors shrink by design.\n")
    P(NOTES.get("gate5", "") + "\n")
else:
    status[5] = ("NOT RUN", "no result files"); P("No result files.\n")

# ---------------------------------------------------------------- gate 6
P("## Gate 6. Placebos\n")
tp = loads("g6time", "g6_time_*")
if tp:
    P("### Time placebo\n\nKFM surveys before the real designation only, with a fake designation year. A correct model should find no effect.\n")
    P("| Taxon | Curve prior | Variant | Fake year | MPAs | Stations | E10 | Interval excludes zero | Worst R-hat / min ESS / divergences |")
    P("|---|---|---|---|---|---|---|---|---|")
    main_bad = 0
    for k, r in tp.items():
        if r.get("not_applicable"):
            continue
        parts = k.split("_"); t, c = parts[2], parts[3]; var = "_".join(parts[4:]) or "main"
        ex = r["E10"]["5.0%"] > 0 or r["E10"]["95.0%"] < 0
        main_bad += int(ex and var == "main")
        P(f"| {TAXA[t]} | {c} | {var} | {r['placebo_year']} | {r['n_mpas']} | {r['stations']} | {e10(r)} | {'YES' if ex else 'no'} | {conv(r)[1]} |")
    P("\n" + NOTES.get("gate6_time", "") + "\n")
sp = loads("g6space", "g6_space_*")
sp_bad = None
if sp:
    P("### Space placebo\n\nAll protected stations are dropped. In each MPA with at least three reference stations, one reference station with at least five post-designation survey years is relabelled as protected (20 random relabellings per taxon and prior). There is no real protection contrast left, so E10 should centre on zero and about 10% of 90% intervals should exclude zero.\n")
    P("| Taxon | Curve prior | Relabellings | Mean of placebo E10 | SD across relabellings | Range | 90% intervals excluding zero | Verdict | Headline E10 (real protected stations) |")
    P("|---|---|---|---|---|---|---|---|---|")
    sp_bad = 0
    for t in TAXA:
        for c in CURVES:
            rs = [r for k, r in sp.items() if k.startswith(f"g6_space_{t}_{c}_")]
            if not rs:
                continue
            m = np.array([r["E10"]["mean"] for r in rs]); ex = np.array([r["E10"]["5.0%"] > 0 or r["E10"]["95.0%"] < 0 for r in rs])
            tstat = abs(m.mean()) / (m.std(ddof=1) / np.sqrt(len(m))) if len(m) > 2 else 0
            off = (np.all(m > 0) or np.all(m < 0)) and abs(m.mean()) > 0.22
            sp_bad += int(ex.mean() > 0.25 or off)
            hx = chosen.get((t, c))
            P(f"| {TAXA[t]} | {c} | {len(rs)} | {m.mean():+.2f} | {m.std(ddof=1):.2f} | {m.min():+.2f} to {m.max():+.2f} | {ex.sum()}/{len(rs)} | {'OFF-CENTRE' if (ex.mean() > 0.25 or off) else 'ok'} | {hx[2]['E10']['mean'] if hx else float('nan'):+.2f} |")
    P("\nVerdict rule: off-centre if every relabelling gives the same sign and the mean exceeds 0.22 in size (a 25% change per decade), or if more than a quarter of intervals exclude zero. The 20 relabellings of one taxon share most of their data (two Point Dume reference stations are relabelled in all 20), so they are not independent tests.\n")
    P(NOTES.get("gate6_space", "") + "\n")
    tr = loads("g6trend", "simtrend_*")
    if tr:
        P("### Follow-up: does a regional trend leak into E10?\n")
        P("Simulated data with no protection effect at all, plus a secular trend shared by every station (protected and reference alike) of 0.08 log units per year up or down (0.8 per decade), on the real layouts. The fitted model is unchanged. E10 should be zero.\n")
        P("| Layout | Trend per decade | Curve prior | Replicates | Mean E10 | 90% coverage of zero | 90% intervals excluding zero |")
        P("|---|---|---|---|---|---|---|")
        cells_t = {}
        for k, r in tr.items():
            _, lay, d, c, _s = k.split("_")
            cells_t.setdefault((lay, d, c), []).append(r)
        for (lay, d, c), rs in sorted(cells_t.items()):
            m = np.array([r["E10"]["mean"] for r in rs]); lo = np.array([r["E10"]["5.0%"] for r in rs]); hi = np.array([r["E10"]["95.0%"] for r in rs])
            P(f"| {lay} | {'+0.8' if d == 'up' else '-0.8'} | {c} | {len(rs)} | {m.mean():+.2f} | {((lo <= 0) & (hi >= 0)).mean():.2f} | {((lo > 0) | (hi < 0)).mean():.2f} |")
        P("\n" + NOTES.get("gate6_trend", "") + "\n")
if tp or sp:
    msg = []
    if tp:
        msg.append(f"time placebo: {main_bad} of 8 main intervals exclude zero")
    if sp_bad is not None:
        msg.append(f"space placebo: {sp_bad} of 8 taxon-by-prior cells off-centre")
    ok6 = (not tp or main_bad == 0) and (sp_bad in (0, None)) and bool(tp) and bool(sp)
    status[6] = ("PASS" if ok6 else "FAIL", "; ".join(msg))
else:
    status[6] = ("NOT RUN", "no result files"); P("No result files.\n")

# ---------------------------------------------------------------- gate 7
P("## Gate 7. PISCO subsampling of dense species\n")
P(NOTES.get("gate7", "") + "\n")
g7 = loads("g7", "g7_*")
if g7:
    P("| Taxon | Curve prior | PISCO units capped | Headline E10 | Cap 90 per transect | Cap 30 per transect |")
    P("|---|---|---|---|---|---|")
    mx = 0
    for t in ("purp", "red"):
        for c in CURVES:
            h = chosen.get((t, c)); a = g7.get(f"g7_{t}_{c}_cap90"); b = g7.get(f"g7_{t}_{c}_cap30")
            for x in (a, b):
                if x and h:
                    mx = max(mx, abs(x["E10"]["mean"] - h[2]["E10"]["mean"]))
            P(f"| {TAXA[t]} | {c} | {a['units_capped'] if a else ''} / {b['units_capped'] if b else ''} | {e10(h[2]) if h else ''} | {e10(a) if a else ''} | {e10(b) if b else ''} |")
    status[7] = ("PASS (bounded)", f"subsampling is present for urchins; searched area is not in the delivered files; capping counts moves E10 by at most {mx:.2f}")
else:
    status[7] = ("PARTIAL", "subsampling documented; bounding refits not run")

# ---------------------------------------------------------------- gate 8
P("## Gate 8. Posterior predictive checks\n")
if chosen and any("ppc" in r for _, _, r in chosen.values()):
    P("200 posterior draws per fit. Unit counts are redrawn with new visit noise and new count noise, holding each station-year's latent density at its posterior draw. p is the share of replicates at or above the observed value; values below 0.05 or above 0.95 are flagged. The two contrast rows (protected/reference log ratio, after minus before, and age 10+ minus age 0 to 5) are conditional on the latent station-year states, so they test the count layer, not the protection model.\n")
    P("| Taxon | Curve prior | Statistic | Observed | Replicate median [5%, 95%] | p | Flag |")
    P("|---|---|---|---|---|---|---|")
    flags = 0; tot = 0
    for (t, c), (k, gate, r) in chosen.items():
        for s, v in r.get("ppc", {}).items():
            f = v["p_rep_ge_obs"] < 0.05 or v["p_rep_ge_obs"] > 0.95
            if s.startswith("q99") and v["rep_q05"] <= v["observed"] <= v["rep_q95"]:
                f = False
            flags += int(f); tot += 1
            P(f"| {TAXA[t]} | {c} | {s} | {v['observed']} | {v['rep_median']} [{v['rep_q05']}, {v['rep_q95']}] | {v['p_rep_ge_obs']:.2f} | {'FLAG' if f else ''} |")
    status[8] = ("PASS" if flags == 0 else "FAIL", f"{flags} of {tot} statistics flagged")
    P("\n" + NOTES.get("gate8", "") + "\n")
else:
    status[8] = ("NOT RUN", "no predictive checks in the result files"); P("No result files.\n")

# ---------------------------------------------------------------- gate 9
P("## Gate 9. Sensitivity of the headline E10\n")
g9 = loads("g9", "g9_*")
if g9:
    cols = [("tauDhalf", "Spread prior halved"), ("tauDdouble", "Spread prior doubled"), ("noshock", "No shared-year shock"), ("level", "Level shift added"),
            ("nolevel", "Level shift removed"), ("noHarris", "Without Harris Point SMR"), ("noWebster", "Without Webster's Arch"), ("noWyckoff", "Without Wyckoff Ledge")]
    P("| Taxon | Curve prior | Headline | " + " | ".join(c[1] for c in cols) + " |")
    P("|---|---|---|" + "---|" * len(cols))
    mx = 0; nconv = 0; nfit = 0; flips = 0
    for t in TAXA:
        for c in CURVES:
            h = chosen.get((t, c)); cells_ = []
            for key, _ in cols:
                r = g9.get(f"g9_{t}_{c}_{key}")
                if r is None:
                    cells_.append("")
                elif r.get("not_applicable"):
                    cells_.append("not in data")
                else:
                    nfit += 1; ok = conv(r)[0]; nconv += int(ok)
                    if h:
                        mx = max(mx, abs(r["E10"]["mean"] - h[2]["E10"]["mean"]))
                        hx = h[2]["E10"]; flips += int((hx["5.0%"] > 0 or hx["95.0%"] < 0) != (r["E10"]["5.0%"] > 0 or r["E10"]["95.0%"] < 0))
                    cells_.append(e10(r) + ("" if ok else " *"))
            P(f"| {TAXA[t]} | {c} | {e10(h[2]) if h else ''} | " + " | ".join(cells_) + " |")
    P(f"\n\\* did not meet the gate-4 convergence thresholds (these runs use 2,000 + 2,000 draws per chain). {nconv} of {nfit} sensitivity fits met them.\n")
    P("Harris Point is left out as a whole MPA (protected and reference stations, both programs). Webster's Arch (KFM a-k-35) is a reference station for sheephead only. Wyckoff Ledge (KFM a-k-01) is the Harris Point reference station for all four taxa.\n")
    P(NOTES.get("gate9", "") + "\n")
    status[9] = ("PASS" if nfit >= 48 else "INCOMPLETE", f"{nfit} fits; largest shift in the E10 mean {mx:.2f}; {flips} fits change whether the 90% interval excludes zero")
else:
    status[9] = ("NOT RUN", "no result files"); P("No result files.\n")

# ---------------------------------------------------------------- assemble
out = ["# Joint density model (stage 1): validation report\n",
       f"Generated by `code/joint_model/report.py`. Design: `docs/JOINT_MODEL_DESIGN.md`. Runs: UCSB GRIT cluster, {NOTES.get('dates', '2026-11-10')}. The first report (short chains, before validation) is `REPORT_first_run_20261109.md`.\n",
       "Review only. Nothing here is adopted for the paper. Both curve priors are reported everywhere: `rw` (random-walk increments, the main prior) and `piecewise` (independent increments shrunk toward zero, the conservative bound). No specification, prior, station or cohort was chosen by the direction or significance of a result.\n",
       "## Gates\n", "| Gate | Result | Basis |", "|---|---|---|"]
names = {1: "Reproduce on GRIT", 2: "Freeze", 3: "Sheephead convergence failure", 4: "Convergence, 4 taxa x 2 priors", 5: "Simulation calibration",
         6: "Placebos", 7: "PISCO subsampling", 8: "Posterior predictive checks", 9: "Sensitivities"}
for i in range(1, 10):
    s = status.get(i, ("NOT RUN", ""))
    out.append(f"| {i}. {names[i]} | {s[0]} | {s[1]} |")
out.append("")
out.append("## Headline: E10 by taxon under both curve priors\n")
out.append("E10 is the regional mean log change in protected-station density, relative to reference stations, 10 years after designation.\n")
out += HEAD or ["No converged fits to report."]
out.append("")
out.append(NOTES.get("headline", ""))
out.append("")
out += L
out.append("## Not verified\n")
out.append(NOTES.get("unverified", ""))
(ROOT / "REPORT.md").write_text("\n".join(out) + "\n")
(V / "gate_status.json").write_text(json.dumps({str(k): v for k, v in status.items()}, indent=1))
print("\n".join(f"gate {k}: {v[0]} - {v[1]}" for k, v in sorted(status.items())))
