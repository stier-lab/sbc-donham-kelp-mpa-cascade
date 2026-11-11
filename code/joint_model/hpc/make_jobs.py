"""Write the task lists for the validation gates (one line per fit: result name, TAB, run.py arguments).

  python3 code/joint_model/hpc/make_jobs.py            # writes code/joint_model/jobs/*.tsv

The sampler settings per taxon and curve prior for gates 6, 7 and 9 are read from jobs/sampler.json, which
records the gate-4 choice. That choice is made on convergence diagnostics only (R-hat, effective sample size,
divergences), never on the estimate.
"""
import json
from pathlib import Path

J = Path(__file__).resolve().parents[1] / "jobs"
J.mkdir(exist_ok=True)
TAXA = {"purp": '--taxon "S. purpuratus"', "red": '--taxon "M. franciscanus"',
        "lob": '--taxon "P. interruptus" --noncentered --level-dev', "sheep": '--taxon "S. pulcher" --noncentered --level-dev'}
CURVES = ["rw", "piecewise"]
LONG = "--warmup 1000 --samples 1000 --chains 4 --cpus 4"
MID = "--warmup 2000 --samples 2000 --chains 4 --cpus 4"
SPACE = "--warmup 1000 --samples 1000 --chains 2"
SHORT = "--warmup 500 --samples 500 --chains 2"


def write(name, rows):
    (J / f"{name}.tsv").write_text("".join(f"{a}.json\t{b}\n" for a, b in rows))
    print(name, len(rows))


# Gate 3: sheephead with the original start (legacy) and the new start (tight), three seeds each
rows = []
for init in ("legacy", "tight"):
    for seed in (1, 2, 3):
        rows.append((f"g3_sheep_rw_{init}_s{seed}", f'fit {TAXA["sheep"]} --curve rw --init {init} --depth 9 --seed {seed} --warmup 600 --samples 600 --chains 4 --cpus 4'))
    rows.append((f"g3_sheep_linear_{init}_s1", f'fit {TAXA["sheep"]} --curve linear --init {init} --depth 9 --seed 1 --warmup 600 --samples 600 --chains 4 --cpus 4'))
write("g3", rows)

# Gate 4, attempt 1: three sampler variants of the same model
VAR = {"base": "", "cb": "--center-b", "dense": "--dense", "cbdense": "--center-b --dense"}
write("g4a", [(f"g4_{t}_{c}_{v}", f"fit {TAXA[t]} --curve {c} {VAR[v]} --seed 1 {LONG} --ppc 200") for t in TAXA for c in CURVES for v in VAR])

# Gate 3, attempt 2: does a small initial step size stop chains flying to the Poisson plateau in early warm-up?
rows = []
for c in CURVES:
    for step0 in ("1.0", "0.01"):
        for seed in range(1, 7):
            rows.append((f"g3b_sheep_{c}_step{step0}_s{seed}", f'fit {TAXA["sheep"]} --curve {c} --step0 {step0} --seed {seed} {LONG}'))
write("g3b", rows)
rows_g3b = list(rows)

# Gate 4, attempt 2: station baselines integrated out (urchins), longer chains, higher target acceptance, small initial step
XL = "--warmup 2000 --samples 2000 --chains 4 --cpus 4 --step0 0.01"
rows = []
for t in ("purp", "red"):
    for c in CURVES:
        rows.append((f"g4_{t}_{c}_margb", f"fit {TAXA[t]} --curve {c} --marg-b --seed 1 {XL} --ppc 200"))
        rows.append((f"g4_{t}_{c}_margb98", f"fit {TAXA[t]} --curve {c} --marg-b --accept 0.98 --seed 1 {XL} --ppc 200"))
        rows.append((f"g4_{t}_{c}_margbdense", f"fit {TAXA[t]} --curve {c} --marg-b --dense --seed 1 {XL} --ppc 200"))
for t in ("lob", "sheep"):
    for c in CURVES:
        rows.append((f"g4_{t}_{c}_base98", f"fit {TAXA[t]} --curve {c} --accept 0.98 --seed 1 {XL} --ppc 200"))
        rows.append((f"g4_{t}_{c}_cb98", f"fit {TAXA[t]} --curve {c} --center-b --accept 0.98 --seed 1 {XL} --ppc 200"))
        rows.append((f"g4_{t}_{c}_margb98", f"fit {TAXA[t].replace(' --noncentered', '')} --curve {c} --marg-b --accept 0.98 --seed 1 {XL} --ppc 200"))
write("g4b", rows)

# Gate 3, attempt 2 continued: same design after the count likelihood was made numerically stable near the Poisson limit
write("g3c", [(a.replace("g3b_", "g3c_"), b) for a, b in rows_g3b])

# Gate 4, final attempt: longer chains and target acceptance 0.99, with the stable likelihood
XXL = "--warmup 3000 --samples 6000 --chains 4 --cpus 4 --step0 0.01 --accept 0.99"
rows = []
for t in ("purp", "red"):
    for c in CURVES:
        rows.append((f"g4_{t}_{c}_margb99", f"fit {TAXA[t]} --curve {c} --marg-b --seed 1 {XXL} --ppc 200"))
        rows.append((f"g4_{t}_{c}_margbdense99", f"fit {TAXA[t]} --curve {c} --marg-b --dense --seed 1 {XXL} --ppc 200"))
for t in ("lob", "sheep"):
    for c in CURVES:
        rows.append((f"g4_{t}_{c}_base99", f"fit {TAXA[t]} --curve {c} --seed 1 {XXL} --ppc 200"))
        rows.append((f"g4_{t}_{c}_cb99", f"fit {TAXA[t]} --curve {c} --center-b --seed 1 {XXL} --ppc 200"))
write("g4c", rows)

# Gate 5: simulation calibration on the real survey layouts
rows = []
for c in CURVES:
    for sc in ("null", "effect", "small", "large"):
        rows += [(f"sim_urchin_{sc}_{c}_s{s}", f'simulate --taxon "S. purpuratus" --truth urchin --scenario {sc} --curve {c} --seed {s} {SHORT}') for s in range(1, 201)]
    for sc in ("null", "lobster_effect"):
        rows += [(f"sim_lobster_{sc}_{c}_s{s}", f'simulate --taxon "P. interruptus" --truth lobster --scenario {sc} --curve {c} --noncentered --level-dev --seed {s} {SHORT}') for s in range(1, 201)]
write("g5", rows)

# Gate 6 follow-up: does a secular regional trend (which the model has no term for) leak into E10 when there is no effect?
rows = []
for c in CURVES:
    for tr_ in ("0.08", "-0.08"):
        tag = "up" if tr_ == "0.08" else "down"
        rows += [(f"simtrend_urchin_{tag}_{c}_s{s}", f'simulate --taxon "S. purpuratus" --truth urchin --scenario null --trend {tr_} --curve {c} --seed {s} {SHORT}') for s in range(1, 41)]
        rows += [(f"simtrend_lobster_{tag}_{c}_s{s}", f'simulate --taxon "P. interruptus" --truth lobster --scenario null --trend {tr_} --curve {c} --noncentered --level-dev --seed {s} {SHORT}') for s in range(1, 41)]
write("g6trend", rows)

sampler = J / "sampler.json"
if sampler.exists():
    S = json.loads(sampler.read_text())   # {"purp": {"rw": "--center-b", "piecewise": ""}, ...}

    def base(t, c):
        return f"{TAXA[t]} --curve {c} {S[t][c]}".strip()

    # Gate 6: placebos. Time: KFM before 2003, fake designation 1992 (1993 sheephead, as in the first report).
    rows = []
    for t in TAXA:
        yr = 1993 if t == "sheep" else 1992
        for c in CURVES:
            nc = "--noncentered" if t in ("lob", "sheep") else ""
            rows.append((f"g6_time_{t}_{c}", f'fit {TAXA[t].split(" --")[0]} {nc} --curve {c} {S[t][c]} --placebo-year {yr} --seed 1 {MID}'))
    # sheephead: which MPA and which fake year drive the off-centre placebo
    sh = '--taxon "S. pulcher" --noncentered ' + S["sheep"]["rw"]
    for c in CURVES:
        for yr in (1990, 1991, 1992, 1994, 1995):
            rows.append((f"g6_time_sheep_{c}_y{yr}", f"fit {sh} --curve {c} --placebo-year {yr} --seed 1 {MID}"))
        for m, tag in (("Harris Point SMR", "noHarris"), ("Santa Barbara Island SMR", "noSBI"), ("Scorpion SMR", "noScorpion")):
            rows.append((f"g6_time_sheep_{c}_{tag}", f'fit {sh} --curve {c} --placebo-year 1993 --omit-mpa "{m}" --seed 1 {MID}'))
        rows.append((f"g6_time_sheep_{c}_leveldev", f"fit {sh} --level-dev --curve {c} --placebo-year 1993 --seed 1 {MID}"))
    write("g6time", rows)
    # Space: protected stations dropped, one reference station per MPA relabelled as protected (20 random relabellings)
    write("g6space", [(f"g6_space_{t}_{c}_s{s}", f"fit {base(t, c)} --space-placebo-seed {s} --seed 1 {SPACE}")
                      for t in TAXA for c in CURVES for s in range(1, 21)])

    # Gate 7: PISCO urchin counts capped at the subsampling threshold (bounding check)
    write("g7", [(f"g7_{t}_{c}_cap{cap}", f"fit {base(t, c)} --pisco-cap {cap} --seed 1 {MID}")
                 for t in ("purp", "red") for c in CURVES for cap in (90, 30)])

    # Gate 9: sensitivities of the headline E10
    rows = []
    for t in TAXA:
        for c in CURVES:
            b = base(t, c)
            rows.append((f"g9_{t}_{c}_tauDhalf", f"fit {b} --tau-d-scale 0.5 --seed 1 {MID}"))
            rows.append((f"g9_{t}_{c}_tauDdouble", f"fit {b} --tau-d-scale 2.0 --seed 1 {MID}"))
            rows.append((f"g9_{t}_{c}_noshock", f"fit {b} --no-shock --seed 1 {MID}"))
            flip = b.replace(" --level-dev", "") if "--level-dev" in b else b + " --level-dev"
            rows.append((f"g9_{t}_{c}_{'nolevel' if '--level-dev' in b else 'level'}", f"fit {flip} --seed 1 {MID}"))
            rows.append((f"g9_{t}_{c}_noHarris", f'fit {b} --omit-mpa "Harris Point SMR" --seed 1 {MID}'))
            rows.append((f"g9_{t}_{c}_noWebster", f'fit {b} --omit-station "KFM|a-k-35" --seed 1 {MID}'))
            rows.append((f"g9_{t}_{c}_noWyckoff", f'fit {b} --omit-station "KFM|a-k-01" --seed 1 {MID}'))
    write("g9", rows)

    # ------------------------------------------------------------------------------------------------------------
    # Structure experiments, 2026-11-10 (review only; nothing adopted). x = experimental options in model.py.
    DR = ("global", "regional")
    SIMX = {"urchin": '--taxon "S. purpuratus" --truth urchin --marg-b --step0 0.01', "lobster": '--taxon "P. interruptus" --truth lobster --noncentered --level-dev --step0 0.01'}
    TREND = {"up": "--trend 0.08", "down": "--trend -0.08", "regvar": "--trend 0.0 --trend-sd 0.08"}
    rows = []
    for lay in SIMX:
        for d in DR:
            for ty in TREND:
                rows += [(f"xA_{lay}_{ty}_rw_{d}_s{s}", f"simulate {SIMX[lay]} --scenario null {TREND[ty]} --curve rw --drift {d} --seed {s} {SHORT}") for s in range(1, 41)]
            rows += [(f"xA_{lay}_up_piecewise_{d}_s{s}", f"simulate {SIMX[lay]} --scenario null {TREND['up']} --curve piecewise --drift {d} --seed {s} {SHORT}") for s in range(1, 41)]
    write("xA", rows)
    rows = []
    for lay, scs in (("urchin", ("null", "effect", "large")), ("lobster", ("null", "lobster_effect"))):
        for sc in scs:
            for d, k in (("global", 1), ("regional", 1), ("regional", 2)):
                rows += [(f"xB_{lay}_{sc}_{d}_inc{k}_s{s}", f"simulate {SIMX[lay]} --scenario {sc} --curve rw --drift {d} --inc-scale {k} --seed {s} {SHORT}") for s in range(1, 101)]
    write("xB", rows)
    write("xC", [(f"xC_{t}_{c}_{d}_s{s}", f"fit {base(t, c)} --drift {d} --space-placebo-seed {s} --seed 1 {SPACE}")
                 for t in TAXA for c in CURVES for d in DR for s in range(1, 21)])
    rows = []
    for t in TAXA:
        yr = 1993 if t == "sheep" else 1992
        nc = "--noncentered" if t in ("lob", "sheep") else ""
        for c in CURVES:
            for d in DR:
                rows.append((f"xD_{t}_{c}_{d}", f'fit {TAXA[t].split(" --")[0]} {nc} --curve {c} {S[t][c]} --drift {d} --placebo-year {yr} --seed 1 {MID}'))
    write("xD", rows)
    XLX = "--warmup 2000 --samples 4000 --chains 4 --cpus 4 --ppc 200 --seed 1"
    rows = []
    for t in TAXA:
        for c in CURVES:
            for d in DR:
                rows.append((f"xE_{t}_{c}_{d}", f"fit {base(t, c)} --drift {d} {XLX}"))
        rows.append((f"xE_{t}_rw_regional_inc2", f"fit {base(t, 'rw')} --drift regional --inc-scale 2 {XLX}"))
        rows.append((f"xE_{t}_rw_regional_three", f"fit {base(t, 'rw')} --drift regional --regions three {XLX}"))
        rows.append((f"xE_{t}_rw_regional_one", f"fit {base(t, 'rw')} --drift regional --regions one {XLX}"))
        rows.append((f"xE_{t}_rw_regional_zi", f"fit {base(t, 'rw')} --drift regional --zi {XLX}"))
        rows.append((f"xE_{t}_rw_none_zi", f"fit {base(t, 'rw')} --zi {XLX}"))
    for t in ("purp", "red"):
        rows.append((f"xE_{t}_piecewise_regional_zi", f"fit {base(t, 'piecewise')} --drift regional --zi {XLX}"))
    write("xE", rows)
