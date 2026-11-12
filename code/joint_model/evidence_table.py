"""How much before-and-after information each MPA and taxon actually has (counts behind the density model).

  python3 code/joint_model/evidence_table.py [uses parquet] [output folder]

Writes validation/evidence_by_mpa.csv (one row per taxon and MPA) and evidence_summary.json. Totals only: animals
counted, area searched and survey years, by protected / reference and before / after designation. No survey records.
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
from data import build

TAXA = {"lob": "P. interruptus", "shp": "S. pulcher", "pur": "S. purpuratus", "red": "M. franciscanus", "kelp": "M. pyrifera"}
uses = sys.argv[1] if len(sys.argv) > 1 else "outputs/integrated_review_20261103/private/supported_measurement_uses.parquet"
out = Path(sys.argv[2] if len(sys.argv) > 2 else "outputs/joint_model_20261109/validation"); out.mkdir(parents=True, exist_ok=True)
rows, summ = [], {}
for g, tx in TAXA.items():
    D = build(uses, tx, regions="three"); c = D["cells"].copy(); c["A"] = c.area * c.n; c["post"] = c.year >= c["T"]
    for m, x in c.groupby("MPA"):
        o = dict(taxon=g, MPA=m, designated=int(x["T"].iloc[0]), programs="/".join(sorted(x.Source.unique())),
                 before_after=int(D["mpa_ba"][D["mpas"].index(m)]))
        for st, k in (("mpa", "prot"), ("reference", "ref")):
            o[f"{k}_stations"] = int(x[x.status == st].station.nunique())
            for p, lab in ((False, "before"), (True, "after")):
                y = x[(x.status == st) & (x.post == p)]
                o[f"{k}_{lab}_years"] = int(y.year.nunique()); o[f"{k}_{lab}_counted"] = int(y.y.sum())
                o[f"{k}_{lab}_per_1000m2"] = round(1000 * y.y.sum() / y.A.sum(), 3) if len(y) else np.nan
        rows.append(o)
    t = pd.DataFrame([r for r in rows if r["taxon"] == g]); ba = t[t.before_after == 1]
    summ[g] = dict(mpas=int(len(t)), before_after_mpas=int(len(ba)), post_only_mpas=int((t.before_after == 0).sum()),
                   protected_counted_before_at_before_after_mpas=int(ba.prot_before_counted.sum()),
                   largest_single_mpa_share_of_that=round(float(ba.prot_before_counted.max() / max(ba.prot_before_counted.sum(), 1)), 2),
                   before_after_mpas_with_under_20_counted_before=int((ba.prot_before_counted < 20).sum()),
                   mpas_with_no_reference_station=int((t.ref_stations == 0).sum()))
pd.DataFrame(rows).to_csv(out / "evidence_by_mpa.csv", index=False)
(out / "evidence_summary.json").write_text(json.dumps(summ, indent=1))
print(json.dumps(summ, indent=1))
