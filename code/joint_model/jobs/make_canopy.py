"""Task lists for the Landsat canopy response (CANOPY_PLAN.md). One line per task: result file name <TAB> run_multi.py arguments."""
from pathlib import Path
here = Path(__file__).parent
PRI = {"rw": "rw", "pw": "piecewise"}
S = "--seed 1 --cpus 4 --chains 4 --warmup 1500 --samples 3000 --thin 2"; L = "--seed 1 --cpus 4 --chains 4 --warmup 3000 --samples 6000 --thin 2"
B = "--taxa can"
fit = [f"can_{s}_{c}.json\tfit {B} --canopy-season {s} --curve {cv} {L}" for s in ("spring", "max", "autumn", "all") for c, cv in PRI.items()]
fit += [f"can_spring_{c}_next.json\tfit {B} --canopy-after next --curve {cv} {L}" for c, cv in PRI.items()]
fit += [f"can_spring_{c}_{d}.json\tfit {B} --canopy-depth {d} --curve {cv} {L}" for c, cv in PRI.items() for d in ("to10m", "over10m")]
fit += [f"can_spring_{c}_pre.json\tfit {B} --pretrend --curve {cv} {L}" for c, cv in PRI.items()]
fit += [f"can_spring_rw_s2.json\tfit {B} --curve rw {L.replace('--seed 1', '--seed 2')}"]
fit += [f"can_time_{c}_y{y}.json\tfit {B} --placebo-year {y} --curve {cv} {S}" for c, cv in PRI.items() for y in (1990, 1993, 1996, 1999)]
perm = [f"can_perm_{c}_k{k:02d}.json\tfit {B} --permute-seed {k} --curve {cv} {S}" for c, cv in PRI.items() for k in range(1, 41)]
sim = [f"can_sim_{n}_s{k}.json\tsimulate {B} {f} --curve rw --seed {k} {S[9:]}".replace("  ", " ") for n, f in (("std", ""), ("null", "--sim-null")) for k in range(201, 221)]
MPAS = ["Abalone Cove SMCA", "Cabrillo SMR", "Campus Point SMCA", "Carrington Pt SMR", "Farnsworth Onshore SMCA", "Gull Island SMR", "Harris Point SMR",
        "Judith Rk SMR", "Naples SMCA", "Point Conception SMR", "Point Dume SMCA", "Point Dume SMR", "Point Vicente SMCA", "Santa Barbara Island SMR",
        "Skunk Pt SMR", "South La Jolla SMR", "South Point SMR", "Swamis SMCA"]
loo = [f"can_loo_{i:02d}.json\tfit {B} --omit-mpa '{m}' --curve rw {S}" for i, m in enumerate(MPAS)]
for name, rows in (("c_fit", fit), ("c_perm", perm), ("c_sim", sim), ("c_loo", loo)):
    (here / f"{name}.tsv").write_text("\n".join(rows) + "\n"); print(name, len(rows))
