#!/usr/bin/env python3
"""Independently check uncertainty exports and preservation of primary results."""

import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import median

from audit_monthly_kelp_refit import bh, keyed, rows

ROOT = Path(__file__).resolve().parents[1]
KEY = ["MPA", "Taxa", "Source", "Resp"]


def close(a, b):
    return math.isclose(float(a), float(b), rel_tol=1e-7, abs_tol=1e-7)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path,
                        default=ROOT / "local_archive/uncertainty_checks_2026-10-12/before")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "docs/PRIMARY_UNCERTAINTY_2026-10-12.json")
    args = parser.parse_args()
    report = {"scope": "Exact primary-contrast HAC and retained-diver temporal uncertainty", "checks": []}
    inputs = ["harmonized_raw_responses.csv", "harmonized_response_ratios.csv",
              "harmonized_landsat_rr.csv", "harmonized_landsat_annual_index.csv",
              "harmonized_site_metadata.csv", "reference_primary_keys.csv", "reference_sites.csv"]
    for name in inputs:
        assert (args.before / "data/harmonized" / name).read_bytes() == (ROOT / "data/harmonized" / name).read_bytes(), name
    assert rows(args.before / "data/sumstats_final.csv") == rows(ROOT / "data/sumstats_final.csv")
    primary = keyed(rows(ROOT / "data/sumstats_final.csv"), KEY)
    assert len(primary) == 145
    preserved = ["table_02_meta_analysis.csv", "table_03_cross_taxa_meta_regression.csv",
                 "table_resistance_recovery.csv", "table_resistance_recovery_sensitivity.csv",
                 "table_s_historical_heatwave_summary.csv", "table_s_historical_heatwave_pairs.csv"]
    for name in preserved:
        assert rows(args.before / "tables" / name) == rows(ROOT / "tables" / name), name
    report["checks"].append("Annual/raw/satellite inputs, reference identities, all 145 primary effects, Tables 2/3 and modern/historical resilience results preserved")
    histories = rows(ROOT / "tables/table_s_primary_hac_histories.csv")
    assert len(histories) == 435
    scenarios = ["lag1", "lag3", "automatic"]
    for scenario in scenarios:
        x = keyed([r for r in histories if r["scenario"] == scenario], KEY)
        assert x.keys() == primary.keys()
        for k, r in x.items():
            p = primary[k]
            assert all(r[f] == p[f] for f in ["Model", "AnalysisType", "N"])
            assert close(r["primary_estimate"], r["reconstructed_estimate"])
            assert close(r["primary_se"], r["reconstructed_se"])
            assert close(r["synthesis_se"], max(float(r["primary_se"]), float(r["hac_se"])))
            assert float(r["synthesis_se"]) >= float(r["primary_se"])
            expected_lag = {"lag1": 1, "lag3": 3,
                            "automatic": max(1, math.floor(4 * (int(r["N"]) / 100) ** (2 / 9)))}[scenario]
            assert int(r["lag_years"]) == expected_lag
    report["checks"].append("435 unique history/bandwidth records cover the exact cohort; contrasts/SEs reconstruct; bandwidth rule and conservative variance floor checked")
    pooled = rows(ROOT / "tables/table_s_primary_hac_pooled.csv")
    assert len(pooled) == 36
    published = keyed(rows(ROOT / "tables/table_02_meta_analysis.csv"), ["Taxa", "Response"])
    baseline = keyed([r for r in pooled if r["scenario"] == "retained_primary"], ["Taxa", "Resp"])
    assert len(baseline) == 9
    for k, r in baseline.items():
        p = published[(k[0], "Biomass" if k[1] == "Bio" else "Density")]
        assert all(close(r[a], p[b]) for a, b in [("estimate", "Estimate"), ("se", "SE"),
                                                ("ci_low", "CI_lower"), ("ci_high", "CI_upper"), ("q", "pval_fdr")])
    for scenario in ["retained_primary"] + ["calendar_HAC_" + s for s in scenarios]:
        group = [r for r in pooled if r["scenario"] == scenario]
        assert len(keyed(group, ["Taxa", "Resp"])) == 9
        for r, q in zip(group, bh([float(r["p"]) for r in group])):
            assert close(r["q"], q)
            base = baseline[(r["Taxa"], r["Resp"])]
            assert float(r["estimate"]) * float(base["estimate"]) > 0
            assert (q < .05) == (float(base["q"]) < .05)
            assert r["direction_changed"] == r["bh_detection_changed"] == "FALSE"
            assert float(r["ci_low"]) < float(r["ci_high"])
    report["checks"].append("Original-variance joint model reproduces Table 2; independent nine-test BH check and unchanged directions/detections in every bandwidth")
    auto = [r for r in histories if r["scenario"] == "automatic"]
    report["automatic_hac"] = {"individual_se_ratio_median": median(float(r["se_ratio"]) for r in auto),
                               "individual_se_ratio_max": max(float(r["se_ratio"]) for r in auto),
                               "n_floored": sum(r["variance_floor_applied"] == "TRUE" for r in auto),
                               "pooled_ci_width_ratios": [float(r["ci_width_ratio"]) for r in pooled if r["scenario"] == "calendar_HAC_automatic"]}
    annual = rows(ROOT / "tables/table_s_temporal_primary_annual.csv")
    assert len(keyed(annual, ["history", "year"])) == 831
    assert len({r["history"] for r in annual}) == 63
    assert len({r["CA_MPA_Name_Short"] for r in annual}) == 15
    taxa = {"Panulirus interruptus": "P. interruptus", "Semicossyphus pulcher": "S. pulcher",
            "Strongylocentrotus purpuratus": "S. purpuratus", "Mesocentrotus franciscanus": "M. franciscanus",
            "Macrocystis pyrifera": "M. pyrifera"}
    for r in annual:
        assert r["source"] != "Landsat" and r["resp"] == "Bio"
        assert (r["CA_MPA_Name_Short"], taxa[r["y"]], r["source"], r["resp"]) in primary
        assert float(r["time"]) >= 0
    fit = keyed(rows(ROOT / "tables/table_s_lmer_model_fit.csv"), ["model"])
    assert fit[("bio",)]["n_obs"] == "831" and fit[("den",)]["n_obs"] == "919"
    for response in ("bio", "den"):
        assert fit[(response,)]["singular_fit"] == "FALSE"
        attempts = rows(ROOT / "tables" / f"table_s_temporal_candidates_{response}.csv")
        assert len(attempts) == 5 and attempts[0]["selected"] == "TRUE"
        assert attempts[0]["singular"] == "FALSE" and attempts[0]["converged"] == "TRUE"
    report["checks"].append("Figure 4 uses 831 rows, 63 retained diver histories, 15 MPAs; response models retain converged nonsingular full random-slope/source structure")
    predictions = rows(ROOT / "tables/table_s_temporal_predictions.csv")
    assert len(keyed(predictions, ["resp", "Species", "time"])) == 1359
    assert all(all(math.isfinite(float(r[f])) for f in ["estimate", "se", "df", "ci_low", "ci_high"]) and
               float(r["df"]) > 0 and float(r["se"]) > 0 and float(r["ci_low"]) < float(r["ci_high"])
               for r in predictions)
    assert len({r["df"] for r in predictions if r["resp"] == "Bio" and r["Species"] == "Macrocystis pyrifera"}) > 100
    attempts = rows(ROOT / "tables/table_s_temporal_ar1_attempts.csv")
    selected = [r for r in attempts if r["selected"] == "TRUE"]
    assert len(selected) == 1 and selected[0]["candidate"] == "independent_slopes_source"
    assert selected[0]["valid"] == selected[0]["pdHess"] == "TRUE" and selected[0]["convergence"] == "0"
    assert all(r["n_obs"] == "831" for r in attempts) and attempts[0]["valid"] == "FALSE"
    diagnostic = selected[0]
    assert diagnostic["n_history"] == "63" and diagnostic["n_mpa"] == "15"
    process_sd = float(diagnostic["ar1_process_sd"])
    annual_sd = float(diagnostic["annual_error_sd"])
    phi = float(diagnostic["ar1_process_phi"])
    conditional_correlation = float(diagnostic["conditional_lag1_correlation"])
    assert process_sd > 0 and annual_sd > 0 and 0 < phi < 1
    assert math.isfinite(float(diagnostic["max_abs_gradient"]))
    assert 0 <= float(diagnostic["source_sd"]) < .001
    assert close(conditional_correlation, phi * process_sd**2 / (process_sd**2 + annual_sd**2))
    assert 0 < conditional_correlation < phi
    slopes = keyed(rows(ROOT / "tables/table_s_temporal_ar1_slopes.csv"), ["Species"])
    assert len(slopes) == 5 and all(r["sign_changed"] == "FALSE" for r in slopes.values())
    purple = slopes[("Strongylocentrotus purpuratus",)]
    kelp = slopes[("Macrocystis pyrifera",)]
    assert float(purple["ci_low"]) < 0 < float(purple["ci_high"]) and float(purple["primary_ci_high"]) < 0
    assert float(kelp["ci_low"]) > 0 and float(kelp["primary_ci_low"]) > 0
    for r in slopes.values():
        assert close(r["ar1_phi"], phi)
        assert close(r["ci_low"], float(r["estimate"]) - 1.959963984540054 * float(r["se"]))
        assert close(r["ci_high"], float(r["estimate"]) + 1.959963984540054 * float(r["se"]))
    report["checks"].append("1,359 exact-contrast prediction rows have finite varying df; same-row AR1 rejects invalid fit and retains positive kelp but purple-urchin interval includes zero")
    report["checks"].append("AR1 process and independent annual-error SDs give the exported conditional lag-one correlation; near-zero source component is recorded separately")
    report["temporal_ar1"] = {"process_phi": phi, "process_sd": process_sd,
                              "annual_error_sd": annual_sd,
                              "conditional_lag1_correlation": conditional_correlation,
                              "source_sd": float(diagnostic["source_sd"]),
                              "max_abs_gradient": float(diagnostic["max_abs_gradient"]),
                              "slopes": list(slopes.values())}
    report["input_md5"] = {name: hashlib.md5((ROOT / "data/harmonized" / name).read_bytes()).hexdigest() for name in inputs}
    report["status"] = "PASS"
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print("PASS: ", len(report["checks"]), " independent preservation, HAC, membership, model and interval checks")


if __name__ == "__main__":
    main()
