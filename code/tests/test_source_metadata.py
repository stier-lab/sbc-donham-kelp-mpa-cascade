"""Regression checks for the recovered source evidence; no network or model fits."""
import csv
import hashlib
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EVIDENCE = ROOT / "data/external/source_metadata_2026-10-16"


def rows(name):
    with (EVIDENCE / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


class SourceEvidenceTest(unittest.TestCase):
    def test_legal_scope_is_separate_from_curation(self):
        data = {r["MPA"]: r for r in rows("no_take_review.csv")}
        self.assertEqual(len(data), 12)
        for name in ["Blue Cavern Onshore SMCA", "Long Point SMR", "Point Dume SMR", "Santa Barbara Island SMR"]:
            self.assertEqual(data[name]["current_rule_class"], "no-take")
        self.assertEqual(data["Swamis SMCA"]["ordinary_sheephead_take"], "not fully prohibited")
        self.assertTrue(all(r["source_url"].startswith("https://wildlife.ca.gov/") for r in data.values()))

    def test_pixel_coordinates_are_recovered_not_inferred(self):
        coverage = rows("satellite_pixel_coverage.csv")
        self.assertEqual(len(coverage), 72)
        self.assertEqual(sum(int(r["n_unmatched"]) for r in coverage), 0)
        overlaps = rows("satellite_reference_overlap.csv")
        focal = [r for r in overlaps if {r["site_a"], r["site_b"]} == {"15", "16"}]
        self.assertEqual(len(focal), 1)
        self.assertEqual(int(focal[0]["n_shared_pixels"]), 320)

    def test_numeric_matches_do_not_hide_missingness(self):
        receipt = json.loads((EVIDENCE / "workbook_comparison.json").read_text())
        raw = ROOT / "data/LANDSAT/MPA_Runs_new.csv"
        self.assertEqual(hashlib.sha256(raw.read_bytes()).hexdigest(), receipt["input_sha256"])
        matches = [r for r in rows("satellite_pixel_comparison.csv")
                   if r["method"] == "sum" and float(r["max_abs_difference"]) < 1e-6]
        self.assertEqual(len(matches), 38)
        self.assertEqual(sum(r["same_missing"] == "TRUE" for r in matches), 6)
        metadata = json.loads((EVIDENCE / "satellite_netcdf_metadata.json").read_text())
        self.assertEqual(metadata["attributes"]["biomass"]["Units"], "wet Kg/900m^2 pixel")
        self.assertEqual(metadata["n_quarters"], 152)

    def test_empty_reference_slots_are_not_observed_zero_biomass(self):
        data = rows("satellite_empty_reference_rows.csv")
        flagged = [r for r in data if r["all_delivered_years_zero"] == "TRUE"]
        self.assertEqual(len(flagged), 6)
        self.assertTrue(all(int(r["n_pixels"]) == 0 for r in flagged))
        arithmetic = rows("satellite_empty_slot_arithmetic.csv")
        ratios = {float(r["candidate_over_delivered"]) for r in arithmetic if r["candidate_over_delivered"]}
        self.assertEqual(ratios, {2.0})


if __name__ == "__main__":
    unittest.main()
