#!/usr/bin/env python3
"""Compare a sheephead-only correction with its preserved input/model checkpoint."""

import argparse
import hashlib
import json
import math
from pathlib import Path

from audit_monthly_kelp_refit import bh, keyed, rows


ROOT = Path(__file__).resolve().parents[1]


def target(row, taxon):
    return row[taxon] == "Semicossyphus pulcher" and row["source"] == "PISCO" and row["resp"] == "Bio"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--report", type=Path, default=ROOT / "docs/SHEEPHEAD_ALLOMETRY_REFIT_2026-10-12.json")
    args = parser.parse_args()
    report = {"scope": "PISCO sheephead biomass only", "checks": [], "input_changes": {}, "primary_changes": []}
    for name, taxon, numbers, extra in (
        ("harmonized_raw_responses.csv", "taxon_name", ["value"], ["status"]),
        ("harmonized_response_ratios.csv", "y", ["mpa", "reference", "Diff", "lnDiff"], [])):
        old_path = args.before / "data/harmonized" / name
        new_path = ROOT / "data/harmonized" / name
        old, new = rows(old_path), rows(new_path)
        keys = ["CA_MPA_Name_Short", "year", taxon, "source", "resp"] + extra
        a, b = keyed(old, keys), keyed(new, keys)
        assert a.keys() == b.keys() and len(old) == len(new)
        old_lines, new_lines = old_path.read_bytes().splitlines(), new_path.read_bytes().splitlines()
        changed = 0
        for i, r in enumerate(old, 1):
            k = tuple(r[f] for f in keys)
            assert all(r[f] == b[k][f] for f in r if f not in numbers)
            if not target(r, taxon):
                assert old_lines[i] == new_lines[i] and r == b[k]
            differs = any(not math.isclose(float(r[f]), float(b[k][f]), rel_tol=1e-12, abs_tol=1e-12) for f in numbers)
            if differs:
                assert target(r, taxon)
                changed += 1
        report["input_changes"][name] = {"rows": len(old), "changed_rows": changed,
            "target_rows": sum(target(r, taxon) for r in old),
            "non_target_rows_byte_preserved": sum(not target(r, taxon) for r in old),
            "after_md5": hashlib.md5(new_path.read_bytes()).hexdigest()}
    for name in ("harmonized_landsat_rr.csv", "harmonized_landsat_annual_index.csv", "harmonized_site_metadata.csv"):
        assert (args.before / "data/harmonized" / name).read_bytes() == (ROOT / "data/harmonized" / name).read_bytes()
    report["checks"].append("All record identities and non-target bytes preserved; satellite/site files unchanged")
    keys = ["MPA", "Taxa", "Source", "Resp"]
    a = keyed(rows(args.before / "data/sumstats_final.csv"), keys)
    b = keyed(rows(ROOT / "data/sumstats_final.csv"), keys)
    assert a.keys() == b.keys() and len(a) == 144
    for k in a:
        assert all(a[k][f] == b[k][f] for f in ("N", "MPA_Start"))
        affected = a[k]["Taxa"] == "S. pulcher" and a[k]["Source"] == "PISCO" and a[k]["Resp"] == "Bio"
        changed = any(not math.isclose(float(a[k][f]), float(b[k][f]), rel_tol=1e-10, abs_tol=1e-10)
                      for f in ("Mean", "SE")) or a[k]["Model"] != b[k]["Model"] or a[k]["AnalysisType"] != b[k]["AnalysisType"]
        if not affected:
            assert not changed, k
        if changed:
            report["primary_changes"].append({"key": k, "before": {f: a[k][f] for f in ("Mean", "SE", "Model", "AnalysisType")},
                                               "after": {f: b[k][f] for f in ("Mean", "SE", "Model", "AnalysisType")}})
    report["checks"].append("All 144 primary keys and sample counts retained; only PISCO sheephead biomass histories can change")
    report["model_or_design_changes"] = [k for k in a if any(a[k][f] != b[k][f] for f in ("Model", "AnalysisType"))]
    for name, ignored in (("reference_primary_keys.csv", set()), ("reference_sites.csv", {"exclusion_definitions_md5"})):
        strip = lambda path: [{k: v for k, v in r.items() if k not in ignored} for r in rows(path)]
        assert strip(args.before / "data/harmonized" / name) == strip(ROOT / "data/harmonized" / name)
    report["checks"].append("Reference identities, exclusions and survey footprints retained after explicit fingerprint refresh")
    old_meta = keyed(rows(args.before / "tables/table_02_meta_analysis.csv"), ["Taxa", "Response"])
    meta = rows(ROOT / "tables/table_02_meta_analysis.csv")
    assert len(meta) == len(old_meta) == 9
    report["pooled_changes"] = []
    for r, q in zip(meta, bh([float(r["pval"]) for r in meta])):
        assert math.isclose(float(r["pval_fdr"]), q, abs_tol=1e-12)
        old = old_meta[(r["Taxa"], r["Response"])]
        report["pooled_changes"].append({"taxon": r["Taxa"], "response": r["Response"],
          "old_estimate": float(old["Estimate"]), "new_estimate": float(r["Estimate"]),
          "new_ratio": math.exp(float(r["Estimate"])), "new_ci": [math.exp(float(r[f])) for f in ("CI_lower", "CI_upper")],
          "new_q": q, "direction_changed": float(r["Estimate"]) * float(old["Estimate"]) < 0,
          "fdr_classification_changed": (q < .05) != (float(old["pval_fdr"]) < .05)})
    report["checks"].append("Nine-test BH adjustment independently reconstructed; inference changes reported explicitly")
    old_pathways = keyed(rows(args.before / "tables/table_03_cross_taxa_meta_regression.csv"), ["Relationship"])
    pathways = keyed(rows(ROOT / "tables/table_03_cross_taxa_meta_regression.csv"), ["Relationship"])
    assert pathways.keys() == old_pathways.keys() and len(pathways) == 6
    report["pathway_changes"] = []
    for key, r in pathways.items():
        old = old_pathways[key]
        assert r["k"] == old["k"]
        if key != ("S. pulcher biomass -> S. purpuratus biomass",):
            assert r == old, key
        if r != old:
            report["pathway_changes"].append({"relationship": key[0], "k": int(r["k"]),
              "before": {f: float(old[f]) for f in ("Slope_Est", "Slope_pval")},
              "after": {f: float(r[f]) for f in ("Slope_Est", "Slope_pval")},
              "nominal_classification_changed": (float(r["Slope_pval"]) < .05) != (float(old["Slope_pval"]) < .05)})
    report["checks"].append("Six pathway memberships retained; the other five regressions are unchanged")
    for name in ("table_resistance_recovery.csv", "table_resistance_recovery_sensitivity.csv",
                 "table_s_historical_heatwave_summary.csv", "table_s_historical_heatwave_pairs.csv"):
        assert rows(args.before / "tables" / name) == rows(ROOT / "tables" / name)
    report["checks"].append("Modern kelp resistance/recovery and historical Landsat results unchanged")
    report["status"] = "PASS"
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
