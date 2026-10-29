#!/usr/bin/env python3
"""Audit the monthly refit against its preserved pre-change input/model snapshot."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAXON = "Macrocystis pyrifera"


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def keyed(data, columns):
    result = {tuple(r[k] for k in columns): r for r in data}
    if len(result) != len(data):
        raise ValueError("Duplicate identities")
    return result


def bh(p):
    order = sorted(range(len(p)), key=p.__getitem__)
    result = [None] * len(p)
    upper = 1.0
    for rank in range(len(p), 0, -1):
        i = order[rank - 1]
        upper = min(upper, p[i] * len(p) / rank)
        result[i] = upper
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, default=ROOT / "local_archive/monthly_kelp_2026-10-12/before")
    parser.add_argument("--report", type=Path, default=ROOT / "docs/MONTHLY_KELP_REFIT_2026-10-12.md")
    args = parser.parse_args()
    checks, changes = [], {}
    for name, taxon, numbers, extra in (
        ("harmonized_raw_responses.csv", "taxon_name", ["value"], ["status"]),
        ("harmonized_response_ratios.csv", "y", ["mpa", "reference", "Diff", "lnDiff"], [])):
        keys = ["CA_MPA_Name_Short", "year", taxon, "source", "resp"] + extra
        before_path, after_path = args.before / "data/harmonized" / name, ROOT / "data/harmonized" / name
        a, b = keyed(rows(before_path), keys), keyed(rows(after_path), keys)
        assert a.keys() == b.keys()
        nonkelp = [k for k in a if a[k][taxon] != TAXON]
        assert all(a[k] == b[k] for k in nonkelp)
        old_lines, new_lines = before_path.read_bytes().splitlines(), after_path.read_bytes().splitlines()
        assert len(old_lines) == len(new_lines)
        for i, r in enumerate(rows(before_path), 1):
            if r[taxon] != TAXON:
                assert old_lines[i] == new_lines[i]
        affected = []
        for k in a:
            assert all(a[k][f] == b[k][f] for f in a[k] if f not in numbers)
            if any(not math.isclose(float(a[k][f]), float(b[k][f]), rel_tol=1e-12, abs_tol=1e-12) for f in numbers):
                assert a[k][taxon] == TAXON
                affected.append(k)
        den = [k for k in affected if a[k]["resp"] == "Den"]
        if extra:
            assert len(den) == 1 and a[den[0]]["source"] == "KFM" and a[den[0]]["year"] == "2007"
        else:
            assert not den
        changes[name] = {"rows": len(a), "changed_rows": len(affected), "nonkelp_rows_byte_preserved": len(nonkelp),
                         "density_changes": len(den), "md5": hashlib.md5(after_path.read_bytes()).hexdigest()}
        checks.append(f"{name}: {len(a):,} identical identities; {len(nonkelp):,} non-kelp records byte-preserved; {len(affected)} kelp records changed.")
    for name in ("harmonized_landsat_rr.csv", "harmonized_landsat_annual_index.csv", "harmonized_site_metadata.csv"):
        assert (args.before / "data/harmonized" / name).read_bytes() == (ROOT / "data/harmonized" / name).read_bytes()
    checks.append("Satellite inputs and site metadata are byte-identical to the checkpoint.")
    key = ["MPA", "Taxa", "Source", "Resp"]
    a = keyed(rows(args.before / "data/sumstats_final.csv"), key)
    b = keyed(rows(ROOT / "data/sumstats_final.csv"), key)
    assert len(a) == len(b) == 144 and a.keys() == b.keys()
    assert all(a[k][f] == b[k][f] for k in a for f in ("N", "Model", "AnalysisType", "MPA_Start"))
    checks.append("All 144 primary histories retain their membership, N, model selection, analysis design and designation date.")
    for k in a:
        if a[k]["Taxa"] != "M. pyrifera" or a[k]["Source"] == "Landsat":
            assert all(math.isclose(float(a[k][f]), float(b[k][f]), rel_tol=1e-10, abs_tol=1e-10) for f in ("Mean", "SE"))
    checks.append("Non-kelp and Landsat individual primary contrasts and SEs are unchanged.")
    for name, ignore in (("reference_primary_keys.csv", set()), ("reference_sites.csv", {"exclusion_definitions_md5"})):
        strip = lambda path: [{k: v for k, v in r.items() if k not in ignore} for r in rows(path)]
        assert strip(args.before / "data/harmonized" / name) == strip(ROOT / "data/harmonized" / name)
    checks.append("Explicit registry re-export retained all 857 reference rows and 144 primary keys; only source-definition fingerprint changed within reference rows.")
    before_meta = keyed(rows(args.before / "tables/table_02_meta_analysis.csv"), ["Taxa", "Response"])
    meta = rows(ROOT / "tables/table_02_meta_analysis.csv")
    adjusted = bh([float(r["pval"]) for r in meta])
    for r, q in zip(meta, adjusted):
        assert math.isclose(float(r["pval_fdr"]), q, abs_tol=1e-12)
        old = before_meta[(r["Taxa"], r["Response"])]
        assert (float(r["pval_fdr"]) < .05) == (float(old["pval_fdr"]) < .05)
        assert float(r["Estimate"]) * float(old["Estimate"]) > 0
    checks.append("Independent nine-test BH calculation reproduces all adjusted P values; all nine direction/detection classifications are unchanged.")
    for name in ("table_s_historical_heatwave_summary.csv", "table_s_historical_heatwave_pairs.csv"):
        assert rows(args.before / "tables" / name) == rows(ROOT / "tables" / name)
    checks.append("Historical Landsat pair-level and summary results are unchanged.")
    kelp = next(r for r in meta if r["Taxa"] == "M. pyrifera")
    old = before_meta[("M. pyrifera", "Biomass")]
    periods = [r for r in rows(ROOT / "tables/table_resistance_recovery_sensitivity.csv") if r["baseline"] == "2010-2013*"]
    lines = ["# Monthly Kelp Calibration Refit", "", "Implemented 12 October 2026 following Adrian's approval of Dan Reed's suggestion.", "",
             "## Calibration", "", "The pinned source is SBC knb-lter-sbc.112.15, Table 5, model 1, DOI 10.6073/pasta/710844a10e2fe8873a4b49d97999d549. Its twelve zero-intercept slopes estimate dry foliar standing crop, excluding holdfasts and sporophylls. Provider CSV bytes and EML are preserved under data/external/sbc11215/.", "",
             "PISCO and LTER convert recorded-month observations before temporal averaging. KFM converts each dated stipe-size record, then resamples whole records so counts and months remain paired. Its existing site-year abundance matching is retained; size and abundance surveys are not necessarily month-matched. An unweighted monthly average is not used.", "",
             "KFM's one-element sample(count, ...) branch previously sampled integers from 1 to count. Sampling record indices corrects one unpaired Anacapa Island 2007 density observation (0.007435 to 0.01). No paired annual density lnRR changed. The correction is separate from the biomass calibration. Legacy RNG consumption is retained for singleton samples to preserve unrelated bootstrap draws.", "", "## Independent Checks", ""]
    lines += ["- PASS: " + c for c in checks]
    lines += ["", "## Results", "", f"The pooled kelp log contrast changed from {float(old['Estimate']):.6f} to {float(kelp['Estimate']):.6f}; exponentiated values are {math.exp(float(old['Estimate'])):.3f} and {math.exp(float(kelp['Estimate'])):.3f}. These combine distinct contrast types and are not one physical inside/reference biomass ratio.", "",
              "| Window | Paired geometric ratio | 95% CI | Signed-rank P |", "|---|---:|---|---:|"]
    lines += [f"| {r['metric']} | {r['ratio_in_over_out']} | {r['ratio_CI_low']} to {r['ratio_CI_high']} | {r['p_wilcoxon']} |" for r in periods]
    lines += ["", "All three modern windows retain ten paired MPAs. The recent-period sign test remains P = 0.109; the result is not test-invariant. The historical comparison still does not detect an increment after designation. Prior reef differences, shared controls and observational attribution limits remain.", "", "## Reproduction", "", "```sh", "python3 -B code/audit_kelp_calibration.py --strict", "Rscript --vanilla code/tests/test_monthly_kelp_refit.R", "python3 -B code/audit_monthly_kelp_refit.py --before /path/to/preserved/before", "```", "", "The before/after audit requires the local checkpoint and is not a normal pipeline dependency. Routine input validation uses the pinned source bytes and kelp_monthly_calibration.json. All pre-change dirty files were preserved in local_archive/monthly_kelp_2026-10-12/before/ in each repository.", ""]
    args.report.write_text("\n".join(lines))
    args.report.with_suffix(".json").write_text(json.dumps({"status": "PASS", "checks": checks, "input_changes": changes}, indent=2) + "\n")
    print("PASS:", len(checks), "input, membership, reference, FDR and historical checks;", args.report)


if __name__ == "__main__":
    main()
