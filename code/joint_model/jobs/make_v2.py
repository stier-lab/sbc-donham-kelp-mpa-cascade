"""Task lists for validation round 2 (VALIDATION_PLAN_v2.md). One line per task: result file name <TAB> run_multi.py arguments (the layout array.slurm reads)."""
from pathlib import Path
here = Path(__file__).parent
TAXA = ["lob", "shp", "pur", "red", "kelp"]; CURVES = {"rw": "rw", "pw": "piecewise"}
S = "--seed 1 --cpus 4 --chains 4 --warmup 1500 --samples 3000 --thin 2"; L = "--seed 1 --cpus 4 --chains 4 --warmup 3000 --samples 6000 --thin 2"
perm = [f"p2_{g}_{c}_k{k:02d}.json\tfit --taxa {g} --curve {cv} --permute-seed {k} {S}" for g in TAXA for c, cv in CURVES.items() for k in range(1, 41)]
split = [f"e2_{g}_{c}_{e}.json\tfit --taxa {g} --curve {cv} --evidence {e} --mass {L}" for g in TAXA for c, cv in CURVES.items() for e in ("ba", "post")]
pre = [f"t2_{g}_{c}_pre.json\tfit --taxa {g} --curve {cv} --pretrend --mass {L}" for g in TAXA for c, cv in CURVES.items()]
pre += [f"t2_{g}_rw_y{y}{'_pre' if f else ''}.json\tfit --taxa {g} --curve rw --placebo-year {y} {f} {S}".replace("  ", " ")
        for g in ("shp", "pur", "red") for y in (1990, 1993, 1996, 1999) for f in ("", "--pretrend")]
sim = [f"s2_{g}_s{k:02d}{'_pre' if f else ''}.json\tsimulate --taxa {g} --scenario nopaths --curve rw --sim-pretrend 0.5 0.3 {f} --seed {k} {S[9:]}".replace("  ", " ")
       for g in ("shp", "pur") for k in range(1, 21) for f in ("", "--pretrend")]
for name, rows in (("v2_perm", perm), ("v2_split", split), ("v2_pre", pre), ("v2_sim", sim)):
    (here / f"{name}.tsv").write_text("\n".join(rows) + "\n"); print(name, len(rows))
