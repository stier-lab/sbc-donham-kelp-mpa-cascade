#!/usr/bin/env python3
"""Audit one authorised membership addition against the pre-inclusion checkpoint."""

import argparse
import hashlib
import json
import math
from pathlib import Path

from audit_monthly_kelp_refit import bh, keyed, rows

ROOT = Path(__file__).resolve().parents[1]
KEY = ("Naples SMCA", "S. pulcher", "LTER", "Bio")


def linear_contrast(annual):
    annual = sorted(annual, key=lambda r: int(r["year"]))
    n_after = sum(r["BA"] == "After" for r in annual)
    x = [0.0] * (len(annual) - n_after) + list(range(n_after))
    y = [float(r["lnDiff"]) for r in annual]
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    ssx = sum((v - mx) ** 2 for v in x)
    slope = sum((a - mx) * (b - my) for a, b in zip(x, y)) / ssx
    intercept = my - slope * mx
    variance = sum((b - intercept - slope * a) ** 2 for a, b in zip(x, y)) / (n - 2)
    return 11 * slope, 11 * math.sqrt(variance / ssx)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--report", type=Path, default=ROOT / "docs/NAPLES_BIOMASS_INCLUSION_2026-10-12.json")
    args = parser.parse_args()
    report = {"scope": "Add Naples LTER sheephead biomass only", "checks": []}
    inputs = ("harmonized_raw_responses.csv", "harmonized_response_ratios.csv",
              "harmonized_landsat_rr.csv", "harmonized_landsat_annual_index.csv", "harmonized_site_metadata.csv")
    for name in inputs:
        assert (args.before / "data/harmonized" / name).read_bytes() == (ROOT / "data/harmonized" / name).read_bytes(), name
    report["checks"].append("All annual, satellite and site inputs byte-preserved")
    cols = ["MPA", "Taxa", "Source", "Resp"]
    old = keyed(rows(args.before / "data/sumstats_final.csv"), cols)
    current = keyed(rows(ROOT / "data/sumstats_final.csv"), cols)
    assert len(old) == 144 and len(current) == 145 and current.keys() - old.keys() == {KEY} and not old.keys() - current.keys()
    for key in old:
        for field in ("Model", "AnalysisType", "N", "MPA_Start", "Mean", "SE", "CI", "SD", "DW_stat", "DW_pval"):
            assert old[key][field] == current[key][field], (key, field)
    candidate = current[KEY]
    assert candidate["Model"] == "Linear" and candidate["AnalysisType"] == "pBACIPS" and int(candidate["N"]) == 23
    annual = [r for r in rows(ROOT / "data/harmonized/harmonized_response_ratios.csv")
              if r["CA_MPA_Name_Short"] == KEY[0] and r["y"] == "Semicossyphus pulcher" and r["source"] == "LTER" and r["resp"] == "Bio"]
    mean, se = linear_contrast(annual)
    assert math.isclose(float(candidate["Mean"]), mean, abs_tol=1e-10)
    assert math.isclose(float(candidate["SE"]), se, abs_tol=1e-10)
    report["new_history"] = {field: candidate[field] for field in (*cols, "Model", "AnalysisType", "N", "Mean", "SE", "DW_stat", "DW_pval")}
    report["checks"].append("Exactly one eligible history added; all 144 existing estimates, designs and diagnostics preserved; new OLS contrast independently reconstructed")
    for name, keys in (("reference_primary_keys.csv", ["series_id"]), ("reference_sites.csv", ["series_id", "reference_id"])):
        a = keyed(rows(args.before / "data/harmonized" / name), keys)
        b = keyed(rows(ROOT / "data/harmonized" / name), keys)
        assert not a.keys() - b.keys()
        assert all(a[k] == b[k] for k in a)
        assert all(k[0] == "|".join(KEY) for k in b.keys() - a.keys())
    report["checks"].append("Original reference identities and footprints retained; new registry records only for the added history")
    previous = keyed(rows(args.before / "tables/table_02_meta_analysis.csv"), ["Taxa", "Response"])
    pooled = rows(ROOT / "tables/table_02_meta_analysis.csv")
    omitted = keyed([r for r in rows(ROOT / "tables/table_s_naples_biomass_inclusion.csv")
                     if r["scenario"] == "omit_naples_lter_sheephead_biomass"], ["Taxa", "Resp"])
    assert len(pooled) == len(previous) == len(omitted) == 9
    report["pooled"] = []
    for r, q in zip(pooled, bh([float(r["pval"]) for r in pooled])):
        a = previous[(r["Taxa"], r["Response"])]
        o = omitted[(r["Taxa"], "Bio" if r["Response"] == "Biomass" else "Den")]
        assert math.isclose(float(r["pval_fdr"]), q, abs_tol=1e-12)
        assert math.isclose(float(o["estimate"]), float(a["Estimate"]), abs_tol=1e-7)
        assert math.isclose(float(o["q"]), float(a["pval_fdr"]), abs_tol=1e-7)
        assert (q < .05) == (float(a["pval_fdr"]) < .05) and float(a["Estimate"]) * float(r["Estimate"]) > 0
        report["pooled"].append({"taxon": r["Taxa"], "response": r["Response"],
            "excluded": math.exp(float(a["Estimate"])), "included": math.exp(float(r["Estimate"])),
            "ci": [math.exp(float(r[k])) for k in ("CI_lower", "CI_upper")], "q": q})
    report["checks"].append("Omission reproduces the preserved primary model; nine-test BH independently checked; all pooled signs/detection classifications retained")
    a = keyed(rows(args.before / "tables/table_03_cross_taxa_meta_regression.csv"), ["Relationship"])
    b = keyed(rows(ROOT / "tables/table_03_cross_taxa_meta_regression.csv"), ["Relationship"])
    assert a.keys() == b.keys() and len(a) == 6
    report["pathways_changed"] = []
    for key in a:
        assert a[key]["k"] == b[key]["k"]
        if a[key] != b[key]:
            assert key == ("S. pulcher biomass -> S. purpuratus biomass",)
            report["pathways_changed"].append({"relationship": key[0], "slope_before": a[key]["Slope_Est"],
                "slope_after": b[key]["Slope_Est"], "p_after": b[key]["Slope_pval"]})
    for name in ("table_resistance_recovery.csv", "table_resistance_recovery_sensitivity.csv",
                 "table_s_historical_heatwave_summary.csv", "table_s_historical_heatwave_pairs.csv"):
        assert rows(args.before / "tables" / name) == rows(ROOT / "tables" / name), name
    report["checks"].append("All pathway MPA memberships, five other regressions and modern/historical kelp resilience results retained")
    report["status"] = "PASS"
    report["input_md5"] = {name: hashlib.md5((ROOT / "data/harmonized" / name).read_bytes()).hexdigest() for name in inputs}
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print("PASS: one history added; 144 preserved; inputs/references preserved; pooled omission, BH, pathways and kelp results checked.")


if __name__ == "__main__":
    main()
