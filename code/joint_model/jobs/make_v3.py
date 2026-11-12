"""Task lists added 2026-11-12 after the first round-2 results (VALIDATION_PLAN_v2.md, addendum).
One line per task: result file name <TAB> run_multi.py arguments."""
from pathlib import Path
here = Path(__file__).parent
TAXA = ["lob", "shp", "pur", "red", "kelp"]
S = "--seed 1 --cpus 4 --chains 4 --warmup 1500 --samples 3000 --thin 2"; L = "--seed 1 --cpus 4 --chains 4 --warmup 3000 --samples 6000 --thin 2"
LN = "--mass --mass-obs lognormal"
logn = [f"l3_{g}_rw_ppc.json\tfit --taxa {g} {LN} --curve rw --mass-ppc 200 {L}" for g in TAXA]
logn += [f"l3_{g}_pw.json\tfit --taxa {g} {LN} --curve piecewise {L}" for g in TAXA]
logn += [f"l3_{g}_hold_s{k}.json\tfit --taxa {g} {LN} --curve rw --mass-holdout {k} --mass-ppc 200 {L}" for g in TAXA for k in (1, 2)]
lsim = [f"l3_sim_std_{g}_s{k}.json\tsimulate --taxa {g} {LN} --mass-scenario std --scenario nopaths --seed {k} --curve rw --cpus 2 --chains 2 --warmup 1500 --samples 1500"
        for g in TAXA for k in range(101, 111)]
clus = [f"c3_shp_rw_k{k:02d}.json\tfit --taxa shp --curve rw --permute-seed {k} --cluster-drift {S}" for k in range(1, 41)]
clus += [f"c3_shp_rw_real_cluster.json\tfit --taxa shp --curve rw --cluster-drift {S}", f"c3_shp_rw_real_base.json\tfit --taxa shp --curve rw {S}"]
for name, rows in (("v3_logn", logn), ("v3_lsim", lsim), ("v3_clus", clus)):
    (here / f"{name}.tsv").write_text("\n".join(rows) + "\n"); print(name, len(rows))
