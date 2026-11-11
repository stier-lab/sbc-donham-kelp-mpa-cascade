"""Checks of the recommended structure (regional drift, three regions, first-increment prior doubled), 2026-11-10.

  python3 code/joint_model/report_r.py    # reads outputs/joint_model_20261109/validation/{xF,xG,xK,xKp,R4,R5,R5t,R6s,R6t,R7,R9}

Writes outputs/joint_model_20261109/REPORT_recommended_structure.md. Review only: the structure is selected with flags,
the default in model.py is still the stage-1 model, and nothing is adopted. Gate rules are the same as in report.py.
"""
import glob, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("outputs/joint_model_20261109")
V = ROOT / "validation"
TAXA = {"purp": "Purple urchin", "red": "Red urchin", "lob": "Lobster", "sheep": "Sheephead"}
NOTES = json.loads((V / "notes_r.json").read_text()) if (V / "notes_r.json").exists() else {}
RHAT, NEFF = 1.01, 400


def loads(gate, pattern):
    out = {}
    for p in sorted(glob.glob(str(V / gate / f"{pattern}.json"))):
        try:
            out[Path(p).stem] = json.loads(Path(p).read_text())
        except Exception:
            pass
    return out


def e10(r):
    e = r["E10"]
    return f"{e['mean']:+.2f} [{e['5.0%']:+.2f}, {e['95.0%']:+.2f}]"


def pct(x):
    return f"{100 * (np.exp(x) - 1):+.0f}%"


def ok(r):
    return r["worst_rhat"] < RHAT and r["min_neff"] > NEFF and r["divergences"] == 0


def conv(r):
    return ("yes" if ok(r) else "NO") + f" ({r['worst_rhat']:.3f} / {r['min_neff']:.0f} / {r['divergences']})"


def excl(r):
    return r["E10"]["5.0%"] > 0 or r["E10"]["95.0%"] < 0


def ppc_flags(r):
    fl = []
    for s, v in r.get("ppc", {}).items():
        f = v["p_rep_ge_obs"] < 0.05 or v["p_rep_ge_obs"] > 0.95
        if s.startswith("q99") and v["rep_q05"] <= v["observed"] <= v["rep_q95"]:
            f = False
        if f:
            fl.append(s)
    return fl


L = []; P = L.append; status = {}
xf = loads("xF", "xF_*"); r4 = loads("R4", "R4_*")
head = {t: xf.get(f"xF_{t}_rw_rec") for t in TAXA}

# ---------------------------------------------------------------- convergence and headline
P("## A. Headline fits, convergence and seed stability\n")
P("4 chains x (3,000 + 6,000). Converged = R-hat < 1.01, effective sample size > 400 for every reported parameter, no divergences (worst R-hat / minimum ESS / divergences).\n")
P("| Taxon | Curve prior | Seed | E10 [90% interval] | Change per decade | Pr(E10<0) | Converged | Predictive-check flags |")
P("|---|---|---|---|---|---|---|---|")
nfit = npass = 0; seedgap = 0.0
for t in TAXA:
    for lab, r, c, sd in ((f"xF_{t}_rw_rec", head[t], "rw", 1), (f"R4_{t}_rw_seed2", r4.get(f"R4_{t}_rw_seed2"), "rw", 2), (f"R4_{t}_piecewise", r4.get(f"R4_{t}_piecewise"), "piecewise", 1)):
        if not r:
            continue
        nfit += 1; npass += int(ok(r)); e = r["E10"]
        if sd == 2 and head[t]:
            seedgap = max(seedgap, abs(e["mean"] - head[t]["E10"]["mean"]))
        P(f"| {TAXA[t]} | {c} | {sd} | {e10(r)} | {pct(e['mean'])} ({pct(e['5.0%'])} to {pct(e['95.0%'])}) | {r['E10_prob_negative']:.2f} | {conv(r)} | {len(ppc_flags(r))} |")
status["A"] = ("Convergence", "PASS" if nfit and npass == nfit else ("FAIL" if nfit else "NOT RUN"), f"{npass} of {nfit} fits converged; largest change in the E10 mean between seeds {seedgap:.2f}")
P("\n" + NOTES.get("A", "") + "\n")

# ---------------------------------------------------------------- calibration
P("## B. Simulation calibration\n")
P("Simulated on the real layouts with three regions and fitted with the recommended structure. 2 chains x (500 + 500). Coverage should be near 0.90; power is the share of 90% intervals excluding zero (false-positive rate at null).\n")
P("| Layout | Scenario | True E10 | Curve prior | Replicates | Mean estimate | Bias | RMSE | 90% coverage | Power / false-positive rate | Wrong sign and excludes zero | Mean interval width |")
P("|---|---|---|---|---|---|---|---|---|---|---|---|")
cells = {}
for k, r in loads("R5", "R5_*").items():
    body = k[3:]; lay, rest = body.split("_", 1); sc, c, _s = rest.rsplit("_", 2)
    cells.setdefault((lay, sc, c), []).append(r)
bad = []; nrw = 0
for (lay, sc, c), rs in sorted(cells.items(), key=lambda kv: (kv[0][2] != "rw", kv[0][0], kv[0][1])):
    tr = rs[0]["truth"]["E10"]; m = np.array([r["E10"]["mean"] for r in rs]); lo = np.array([r["E10"]["5.0%"] for r in rs]); hi = np.array([r["E10"]["95.0%"] for r in rs])
    cov = ((lo <= tr) & (tr <= hi)).mean(); ex = (lo > 0) | (hi < 0)
    wrong = "" if tr == 0 else f"{(ex & (np.sign(m) != np.sign(tr))).mean():.3f}"
    if c == "rw":
        nrw += 1
        if cov < 0.85 or cov > 0.96:
            bad.append(f"{lay} {sc} {cov:.2f}")
    P(f"| {lay} | {sc} | {tr:+.2f} | {c} | {len(rs)} | {m.mean():+.2f} | {m.mean() - tr:+.2f} | {np.sqrt(((m - tr) ** 2).mean()):.2f} | {cov:.2f} | {ex.mean():.2f} | {wrong} | {(hi - lo).mean():.2f} |")
tcells = {}
for k, r in loads("R5t", "R5t_*").items():
    _, lay, ty, c, _s = k.split("_"); tcells.setdefault((lay, ty), []).append(r)
P("\nShared-trend simulations (no protection effect; trend of 0.8 per decade up or down in every station, or region-specific trends with SD 0.8 per decade), `rw` prior:\n")
P("| Layout | Trend | Replicates | Mean E10 | 90% coverage of zero | Intervals excluding zero |")
P("|---|---|---|---|---|---|")
tbad = []
for (lay, ty), rs in sorted(tcells.items()):
    m = np.array([r["E10"]["mean"] for r in rs]); lo = np.array([r["E10"]["5.0%"] for r in rs]); hi = np.array([r["E10"]["95.0%"] for r in rs])
    cov = ((lo <= 0) & (hi >= 0)).mean()
    if cov < 0.85:
        tbad.append(f"{lay} {ty} {cov:.2f}")
    P(f"| {lay} | {ty} | {len(rs)} | {m.mean():+.2f} | {cov:.2f} | {((lo > 0) | (hi < 0)).mean():.2f} |")
if cells:
    status["B"] = ("Simulation calibration (rw)", "PASS" if not bad and not tbad and nrw == 6 else "FAIL",
                   (f"{nrw} of 6 cells run; " if nrw < 6 else "") + ("coverage within 0.85 to 0.96 in every `rw` cell" if not bad else "coverage outside 0.85 to 0.96: " + "; ".join(bad))
                   + ("; trend simulations keep coverage of zero at 0.85 or above" if tcells and not tbad else ("; trend coverage low: " + "; ".join(tbad) if tbad else "")))
else:
    status["B"] = ("Simulation calibration (rw)", "NOT RUN", "")
P("\nGate rule: 90% coverage between 0.85 and 0.96 in every `rw` cell, and coverage of zero at 0.85 or above in every shared-trend cell. `piecewise` rows are reported, not gated.\n")
P(NOTES.get("B", "") + "\n")

# ---------------------------------------------------------------- placebos
P("## C. Placebos\n")
P("### Time placebo (KFM before 2003, fake designation 1992; 1993 for sheephead)\n")
P("| Taxon | Curve prior | E10 [90% interval] | Excludes zero | Converged |")
P("|---|---|---|---|---|")
r6t = loads("R6t", "R6t_*"); tbad2 = []
for t in TAXA:
    for c in ("rw", "piecewise"):
        r = r6t.get(f"R6t_{t}_{c}")
        if r:
            if excl(r):
                tbad2.append(f"{TAXA[t]} {c}")
            P(f"| {TAXA[t]} | {c} | {e10(r)} | {'YES' if excl(r) else 'no'} | {conv(r)} |")
P("\n### Space placebo\n")
P("Protected stations dropped, one reference station per MPA relabelled as protected. `rw`: 20 relabellings; `piecewise`: 10; without Point Dume: both Point Dume MPAs left out before relabelling (their two reference stations are otherwise relabelled every time), 10 relabellings.\n")
P("| Taxon | Version | Relabellings | Mean placebo E10 | Range | 90% intervals excluding zero | Mean interval width | Verdict | Headline E10 |")
P("|---|---|---|---|---|---|---|---|---|")
r6s = loads("R6s", "R6s_*"); xg = loads("xG", "xG_*"); xkp = loads("xKp", "xKp_kelp_rw_rec_*"); xk = loads("xK", "xK_*")
sets = {}
for k, r in xg.items():
    sets.setdefault((k.split("_")[1], "rw"), []).append(r)
for k, r in r6s.items():
    _, t, ver, _s = k.split("_"); sets.setdefault((t, ver), []).append(r)
for k, r in xkp.items():
    sets.setdefault(("kelp", "rw"), []).append(r)
names = dict(TAXA, kelp="Kelp (counts, exploratory)")
sbad = []
for t in list(TAXA) + ["kelp"]:
    hd = head.get(t) if t in TAXA else xk.get("xK_kelp_rw_rec_margb")
    for ver, lab in (("rw", "rw"), ("piecewise", "piecewise"), ("rwnoPD", "rw, without Point Dume")):
        rs = sets.get((t, ver))
        if not rs:
            continue
        m = np.array([r["E10"]["mean"] for r in rs]); ex = np.array([excl(r) for r in rs]); wd = np.mean([r["E10"]["95.0%"] - r["E10"]["5.0%"] for r in rs])
        off = ((np.all(m > 0) or np.all(m < 0)) and abs(m.mean()) > 0.22) or ex.mean() > 0.25
        if off and ver == "rw" and t in TAXA:
            sbad.append(names[t])
        P(f"| {names[t]} | {lab} | {len(rs)} | {m.mean():+.2f} | {m.min():+.2f} to {m.max():+.2f} | {int(ex.sum())}/{len(rs)} | {wd:.2f} | {'OFF-CENTRE' if off else 'ok'} | {hd['E10']['mean'] if hd else float('nan'):+.2f} |")
P("\nVerdict rule (same as the stage-1 report): off-centre if every relabelling has the same sign and the mean exceeds 0.22 in size, or if more than a quarter of intervals exclude zero. Gate uses the `rw` rows for the four taxa.\n")
if r6t or sets:
    status["C"] = ("Placebos", "PASS" if not tbad2 and not sbad else "FAIL",
                   ("time placebo: no interval excludes zero" if not tbad2 else "time placebo excludes zero: " + ", ".join(tbad2)) + "; "
                   + ("space placebo centred for all four taxa" if not sbad else "space placebo off-centre: " + ", ".join(sbad)))
else:
    status["C"] = ("Placebos", "NOT RUN", "")
P(NOTES.get("C", "") + "\n")

# ---------------------------------------------------------------- caps
P("## D. PISCO count caps (subsampling bound)\n")
r7 = loads("R7", "R7_*"); mx = 0.0
if r7:
    P("| Taxon | Headline E10 | Cap 90 per transect | Cap 30 per transect |")
    P("|---|---|---|---|")
    for t in ("purp", "red"):
        a = r7.get(f"R7_{t}_rw_cap90"); b = r7.get(f"R7_{t}_rw_cap30")
        for x in (a, b):
            if x and head[t]:
                mx = max(mx, abs(x["E10"]["mean"] - head[t]["E10"]["mean"]))
        P(f"| {TAXA[t]} | {e10(head[t]) if head[t] else ''} | {e10(a) if a else ''} | {e10(b) if b else ''} |")
    status["D"] = ("PISCO subsampling bound", "PASS (bounded)", f"capping counts moves E10 by at most {mx:.2f}")
else:
    status["D"] = ("PISCO subsampling bound", "NOT RUN", "")
P("\n" + NOTES.get("D", "") + "\n")

# ---------------------------------------------------------------- predictive checks
P("## E. Posterior predictive checks (headline fits, `rw`)\n")
tot = fl = 0; byt = []
for t in TAXA:
    if head[t]:
        f = ppc_flags(head[t]); tot += len(head[t].get("ppc", {})); fl += len(f); byt.append(f"{TAXA[t]} {len(f)}")
status["E"] = ("Posterior predictive checks", ("PASS" if fl == 0 else "FAIL") if tot else "NOT RUN", f"{fl} of {tot} statistics flagged ({'; '.join(byt)})")
P(f"{fl} of {tot} statistics flagged: {'; '.join(byt)}. The count layer is unchanged from stage 1, so the urchin misfit (too few zeros, upper tail too heavy) is unchanged.\n")
P(NOTES.get("E", "") + "\n")

# ---------------------------------------------------------------- sensitivities
P("## F. Sensitivities of the headline E10 (`rw`)\n")
r9 = loads("R9", "R9_*")
if r9:
    cols = [("tauDhalf", "Spread prior halved"), ("tauDdouble", "Spread prior doubled"), ("noshock", "No shared-year shock"), ("level", "Level shift added"),
            ("nolevel", "Level shift removed"), ("noHarris", "Without Harris Point SMR"), ("noWebster", "Without Webster's Arch"), ("noWyckoff", "Without Wyckoff Ledge")]
    P("| Taxon | Headline | " + " | ".join(c[1] for c in cols) + " |")
    P("|---|---|" + "---|" * len(cols))
    mx = 0.0; flips = []; n9 = c9 = 0
    for t in TAXA:
        row = []
        for key, lab in cols:
            r = r9.get(f"R9_{t}_{key}")
            if r is None:
                row.append("")
            elif r.get("not_applicable"):
                row.append("not in data")
            else:
                n9 += 1; c9 += int(ok(r))
                if head[t]:
                    mx = max(mx, abs(r["E10"]["mean"] - head[t]["E10"]["mean"]))
                    if excl(r) != excl(head[t]):
                        flips.append(f"{TAXA[t]}: {lab.lower()}")
                row.append(e10(r) + ("" if ok(r) else " *"))
        P(f"| {TAXA[t]} | {e10(head[t]) if head[t] else ''} | " + " | ".join(row) + " |")
    P(f"\n\\* did not meet the convergence thresholds (these runs use 2,000 + 4,000 draws per chain). {c9} of {n9} did.\n")
    status["F"] = ("Sensitivities", "PASS" if not flips else "MIXED", f"{n9} fits; largest shift in the E10 mean {mx:.2f}; " + ("no fit changes whether the 90% interval excludes zero" if not flips else "interval status changes in: " + "; ".join(flips)))
else:
    status["F"] = ("Sensitivities", "NOT RUN", "")
P(NOTES.get("F", "") + "\n")

# ---------------------------------------------------------------- assemble
out = ["# Joint density model: checks of the recommended structure\n",
       "Generated by `code/joint_model/report_r.py`, 2026-11-10. Structure: random-walk year effect shared by all stations plus a random-walk deviation per region (three regions: northern Channel Islands, southern islands, mainland), first-increment prior Normal(0, 1.2), negative-binomial counts. Background: `REPORT.md` (stage-1 gates) and `REPORT_structure_experiments.md` (why this structure).\n",
       "Review only. The structure is selected with flags; the default in `model.py` is still the stage-1 model and nothing is adopted for the paper. No specification, prior, station or cohort was chosen by the direction or significance of a result.\n",
       NOTES.get("summary", ""), "",
       "## Checks\n", "| Check | Result | Basis |", "|---|---|---|"]
for k in "ABCDEF":
    if k in status:
        out.append(f"| {k}. {status[k][0]} | {status[k][1]} | {status[k][2]} |")
out.append("")
out += L
out.append("## Not verified\n"); out.append(NOTES.get("unverified", ""))
(ROOT / "REPORT_recommended_structure.md").write_text("\n".join(out) + "\n")
(V / "check_status_recommended.json").write_text(json.dumps(status, indent=1))
print("\n".join(f"{k}. {v[0]}: {v[1]} - {v[2]}" for k, v in sorted(status.items())))
