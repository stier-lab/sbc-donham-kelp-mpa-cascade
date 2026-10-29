"""Compare the corrected run with the dated pre-correction checkpoint."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    args = parser.parse_args()
    before = args.before
    old = pd.read_csv(before / "sumstats_final.csv")
    new = pd.read_csv(ROOT / "data/sumstats_final.csv")
    keys = ["MPA", "Taxa", "Source", "Resp"]
    paired = old.merge(new, on=keys, suffixes=("_before", "_after"), validate="one_to_one")
    assert len(paired) == len(old) == len(new) == 145
    diver = paired.Source != "Landsat"
    for field in ["Mean", "SE", "SD", "CI", "N"]:
        assert np.allclose(paired.loc[diver, field + "_before"], paired.loc[diver, field + "_after"],
                           equal_nan=True, rtol=0, atol=1e-8), field
    changed = paired[(abs(paired.Mean_before - paired.Mean_after) > 1e-8) |
                     (abs(paired.SE_before - paired.SE_after) > 1e-8)]
    assert len(changed) == 6 and (changed.Source == "Landsat").all()
    old_pool = pd.read_csv(before / "tables/table_02_meta_analysis.csv")
    new_pool = pd.read_csv(ROOT / "tables/table_02_meta_analysis.csv")
    pooled = old_pool.merge(new_pool, on=["Taxa", "Response"], suffixes=("_before", "_after"), validate="one_to_one")
    assert (np.sign(pooled.Estimate_before) == np.sign(pooled.Estimate_after)).all()
    assert ((pooled.pval_fdr_before < .05) == (pooled.pval_fdr_after < .05)).all()
    unchanged_inputs = {}
    for name in ["harmonized_raw_responses.csv", "harmonized_response_ratios.csv", "harmonized_site_metadata.csv"]:
        unchanged_inputs[name] = ((before / "harmonized" / name).read_bytes() ==
                                  (ROOT / "data/harmonized" / name).read_bytes())
    assert all(unchanged_inputs.values())
    historical = pd.read_csv(ROOT / "tables/table_s_historical_heatwave_summary.csv")
    historical = historical[(historical.specification == "primary") & (historical.Region == "All")]
    modern_old = pd.read_csv(before / "tables/table_resistance_recovery.csv")
    modern_new = pd.read_csv(ROOT / "tables/table_resistance_recovery.csv")
    pd.testing.assert_frame_equal(modern_old, modern_new, check_exact=False, atol=1e-10, rtol=0)
    kelp = new_pool[(new_pool.Taxa == "M. pyrifera") & (new_pool.Response == "Biomass")].iloc[0]
    marker = json.loads((ROOT / "data/harmonized/landsat_reference_correction.json").read_text())
    matched = pd.read_csv(ROOT / "tables/table_s_design_matched_summary.csv")
    matched = matched[(matched.correlation == "Annual means") & (matched.source_subset == "All programmes")
        & (matched.Region == "All") & (matched.Taxa == "M. pyrifera") & (matched.Resp == "Bio")]
    change = matched[matched.design == "BACIPS_step"].iloc[0]
    post = matched[matched.design == "CI_after"].iloc[0]
    assert change.p_fdr < .05 and post.p_fdr >= .05
    report = {"release": marker["release"], "primary_cohort": len(new),
        "old_satellite_pairs": marker["old_rr_rows"], "new_satellite_pairs": marker["new_rr_rows"],
        "diver_inputs_unchanged": unchanged_inputs, "diver_primary_estimates_unchanged": True,
        "all_nine_signs_and_BH_classifications_unchanged": True,
        "modern_diver_resilience_unchanged": True,
        "changed_primary_rows": changed[keys + ["Mean_before", "Mean_after", "SE_before", "SE_after", "N_before", "N_after"]].to_dict(orient="records"),
        "pooled_before_after": pooled.to_dict(orient="records"),
        "historical_primary": historical.to_dict(orient="records"),
        "matched_kelp": matched[["design", "estimate", "ci_low", "ci_high", "p_fdr", "n_mpa"]].to_dict(orient="records"),
        "kelp_exponentiated_contrast": {"estimate": float(np.exp(kelp.Estimate)),
            "ci_low": float(np.exp(kelp.CI_lower)), "ci_high": float(np.exp(kelp.CI_upper)),
            "p_fdr": float(kelp.pval_fdr)}}
    path = ROOT / "docs/LANDSAT_REFERENCE_REFIT_2026-10-16.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    lines = ["# Empty Landsat reference correction", "",
        "Six documented empty second-reference slots were removed before averaging.",
        "The source CSV is preserved. Observed zeros in nonempty slots remain.",
        "", f"Paired satellite MPA-years: {marker['old_rr_rows']} -> {marker['new_rr_rows']}.",
        "South La Jolla 1991 loses its artificial zero reference and is no longer a paired observation.",
        "The 145 primary histories and MPA membership are unchanged. All diver inputs",
        "and diver primary effect estimates are unchanged. Six Landsat estimates change.",
        "Shared variance components cause small changes in other pooled biomass estimates.",
        "All nine pooled signs and BH classifications persist.", "",
        f"Kelp's exponentiated pooled contrast is {np.exp(kelp.Estimate):.2f} "
        f"(95% CI {np.exp(kelp.CI_lower):.2f}-{np.exp(kelp.CI_upper):.2f}; adjusted P = {kelp.pval_fdr:.4f}).",
        "This is an exponentiated synthesis of design-specific contrasts, not a common physical abundance ratio.",
        "The modern diver resistance/recovery results are unchanged.", "",
        f"The secondary matched before/after kelp check now meets its adjusted threshold: "
        f"log change {change.estimate:.3f}, 95% CI {change.ci_low:.3f}-{change.ci_high:.3f}, "
        f"adjusted P = {change.p_fdr:.3f}, {change.n_mpa} MPAs (previous adjusted P = 0.054).",
        f"Its post-only contrast remains uncertain (adjusted P = {post.p_fdr:.3f}).",
        "This secondary change does not establish a predator-mediated cascade or erase placement confounding.", "",
        "No-take rules and inherited programme-specific inclusion are documented separately.",
        "Gull Island KFM remains primary; its omission sensitivity remains.",
        "The missing habitat-specific classification list and exact extraction masks/coverage rules",
        "are still unresolved. This correction does not reclassify other histories.", "",
        "Rebuild satellite inputs explicitly through `code/R/export_landsat_annual_index.R`",
        "with the source workspace path. After the primary refit, explicitly rebuild the",
        "reference registry with `code/R/export_reference_registry.R` and the same source",
        "workspace path, then rerun `Rscript code/R/run_all.R` to validate all dependencies.",
        "The correction marker is validated before analyses and cached figure regeneration.",
        "The before/after comparison and focused tests are reproducible from this checkpoint", 
        "with `python3 code/audit_landsat_reference_refit.py --before <checkpoint>`."]
    path.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print("PASS: six changed satellite fits; 145 histories; diver inputs/estimates unchanged; nine pooled classifications unchanged")
    print(lines[12])


if __name__ == "__main__":
    main()
