"""Synthetic tests only: no live inputs, pipeline models, or output writes."""

import importlib.util
import csv
import io
import math
from pathlib import Path
import sys
import tempfile
import unittest


MODULE = Path(__file__).resolve().parents[1] / "audit_lter_count_refit.py"
spec = importlib.util.spec_from_file_location("lter_audit", MODULE)
audit = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = audit
spec.loader.exec_module(audit)


def snapshot(text):
    return audit.CsvSnapshot.parse(Path("synthetic.csv"), text.encode("utf-8"))


def failures(report):
    return [item["name"] for item in report.checks if item["status"] == "FAIL"]


class CsvByteTests(unittest.TestCase):
    def test_multiline_crlf_and_quoted_records_preserved(self):
        data = b'key,text,value\r\n1,"two\r\nlines, ""quoted""",2\r\n2,last,3\r\n'
        parsed = audit.CsvSnapshot.parse(Path("synthetic.csv"), data)
        self.assertEqual(parsed.rows[0]["text"], 'two\r\nlines, "quoted"')
        self.assertEqual(parsed.header + b"".join(parsed.records), data)
        self.assertEqual(parsed.records[1], b"2,last,3\r\n")

    def test_duplicate_key_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicate key"):
            snapshot("key,value\n1,a\n1,b\n").indexed(("key",))

    def test_missing_key_column_rejected(self):
        with self.assertRaisesRegex(ValueError, "missing key"):
            snapshot("key,value\n1,a\n").indexed(("absent",))

    def test_malformed_record_rejected(self):
        with self.assertRaisesRegex(ValueError, "malformed record"):
            snapshot("key,value\n1,a,b\n")

    def test_duplicate_column_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicate column"):
            snapshot("key,key\n1,2\n")


class InputTests(unittest.TestCase):
    HEADER = "CA_MPA_Name_Short,year,taxon_name,source,resp,status,value,BA\n"
    CAMPUS = "Campus Point SMCA,2020,Semicossyphus pulcher,LTER,Bio,reference,3,After\n"
    OTHER = "Campus Point SMCA,2020,Semicossyphus pulcher,LTER,Den,reference,1,After\n"
    KEYS = ("CA_MPA_Name_Short", "year", "taxon_name", "source", "resp", "status")

    def compare(self, new):
        report = audit.Audit()
        audit.compare_inputs(report, snapshot(self.HEADER+self.CAMPUS+self.OTHER),
                             snapshot(new), self.KEYS, {"value"}, "raw")
        return report

    def test_only_target_value_can_change(self):
        new = self.HEADER+self.CAMPUS.replace(",3,After", ",4,After")+self.OTHER
        self.assertEqual(failures(self.compare(new)), [])

    def test_non_target_formatting_change_fails_byte_check(self):
        new = (self.HEADER+self.CAMPUS.replace(",3,After", ",4,After")
               +self.OTHER.replace(",1,After", ",1.0,After"))
        self.assertIn("raw: all non-target records byte-identical", failures(self.compare(new)))

    def test_target_metadata_change_fails(self):
        new = self.HEADER+self.CAMPUS.replace(",3,After", ",4,Before")+self.OTHER
        self.assertIn("raw: target metadata preserved", failures(self.compare(new)))

    def test_missing_non_target_record_fails(self):
        new = self.HEADER+self.CAMPUS.replace(",3,After", ",4,After")
        self.assertIn("raw: identical record membership", failures(self.compare(new)))

    def test_no_correction_is_flagged_as_stale(self):
        self.assertIn("raw: corrected target values present",
                      failures(self.compare(self.HEADER+self.CAMPUS+self.OTHER)))


class PrimaryTests(unittest.TestCase):
    HEADER = "MPA,Taxa,Source,Resp,Model,AnalysisType,Mean,SE,CI,N\n"
    CAMPUS = "Campus Point SMCA,S. pulcher,LTER,Bio,Linear,pBACIPS,1,0.1,0.2,22\n"
    OTHER = "Elsewhere,S. pulcher,PISCO,Bio,Mean,CI,2,0.2,0.4,10\n"

    def compare(self, new):
        report = audit.Audit()
        audit.compare_primary(report, snapshot(self.HEADER+self.CAMPUS+self.OTHER), snapshot(new))
        return report

    def test_campus_only_numeric_correction_passes(self):
        new = self.HEADER+self.CAMPUS.replace(",1,0.1,0.2,", ",1.1,0.11,0.22,")+self.OTHER
        report = self.compare(new)
        self.assertEqual(failures(report), [])
        self.assertEqual(report.details["primary"]["counts"]["after"]["total"], 2)

    def test_model_change_does_not_pass_as_numeric_correction(self):
        new = self.HEADER+self.CAMPUS.replace("Linear", "Mean")+self.OTHER
        self.assertIn("Primary: only Campus LTER sheephead Bio numerics changed",
                      failures(self.compare(new)))

    def test_non_target_numeric_change_fails(self):
        new = self.HEADER+self.CAMPUS.replace(",1,0.1,", ",1.1,0.1,")+self.OTHER.replace(",2,0.2,", ",2.1,0.2,")
        self.assertIn("Primary: only Campus LTER sheephead Bio numerics changed",
                      failures(self.compare(new)))

    def test_naples_admission_fails_membership_and_routing(self):
        naples = self.CAMPUS.replace("Campus Point", "Naples")
        report = self.compare(self.HEADER+self.CAMPUS+self.OTHER+naples)
        self.assertIn("Primary: identical membership", failures(report))
        self.assertIn("Primary: Campus retained; Naples still excluded", failures(report))


class PhaseTwoTests(unittest.TestCase):
    def table(self, step_fixed=False, campus_changed=False, bad_se=False, missing_flip=False):
        text = io.StringIO(newline="")
        writer = csv.writer(text)
        writer.writerow(("MPA", "Taxa", "Source", "Resp", "Model", "AnalysisType", "Mean", "SE", "CI", "N"))
        writer.writerow((*audit.CAMPUS_KEY, "Linear", "pBACIPS", 1.1 if campus_changed else 1, .1, .2, 22))
        for i, key in enumerate(sorted(audit.STEP_KEYS)):
            mean = i + 2
            if step_fixed and not (missing_flip and i == 0):
                mean = -mean
            writer.writerow((*key, "Step", "pBACIPS", mean, .3 if bad_se and i == 0 else .2, .4, 23))
        return snapshot(text.getvalue())

    def test_phase_two_requires_exact_six_sign_reversals_plus_campus(self):
        report = audit.Audit()
        audit.compare_primary(report, self.table(), self.table(True, True), step_direction_fixed=True)
        self.assertEqual(failures(report), [])

    def test_phase_two_uncertainty_change_fails(self):
        report = audit.Audit()
        audit.compare_primary(report, self.table(), self.table(True, True, bad_se=True), step_direction_fixed=True)
        self.assertIn("Primary: only authorized Campus/six-Step changes", failures(report))

    def test_phase_two_missing_flip_fails(self):
        report = audit.Audit()
        audit.compare_primary(report, self.table(), self.table(True, True, missing_flip=True), step_direction_fixed=True)
        self.assertIn("Primary: all and only six stored Step means reversed", failures(report))

    def test_count_only_to_combined_campus_must_stay_fixed(self):
        report = audit.Audit()
        audit.compare_count_only_primary(report, self.table(campus_changed=True), self.table(True, True))
        self.assertEqual(failures(report), [])
        bad = audit.Audit()
        audit.compare_count_only_primary(bad, self.table(campus_changed=True), self.table(True, False))
        self.assertIn("Phase 2 vs COUNT-only: six Mean sign reversals, all other fields unchanged", failures(bad))


class NumericTests(unittest.TestCase):
    def test_bh_unsorted_ties_and_bound(self):
        actual = audit.bh_adjust([.04, .001, .04, .9, 1])
        expected = [1/15, .005, 1/15, 1, 1]
        for left, right in zip(actual, expected):
            self.assertAlmostEqual(left, right)

    def test_bh_empty(self):
        self.assertEqual(audit.bh_adjust([]), [])

    def test_invalid_pvalues_rejected(self):
        for value in (-.01, 1.01, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                audit.bh_adjust([value])

    def test_tolerance_and_missing_values(self):
        self.assertTrue(audit.same_value("1", "1.000000000001"))
        self.assertFalse(audit.same_value("1", "1.001"))
        self.assertTrue(audit.same_value("NA", "NaN"))
        self.assertFalse(audit.same_value("NA", "0"))

    def test_bad_fdr_is_rejected(self):
        table = snapshot("Taxa,Response,Estimate,SE,CI_lower,CI_upper,pval,pval_fdr,RR,RR_lower,RR_upper\n"
                         "species,Biomass,0,1,-1,1,0.04,0.9,1,0.3678794411714423,2.718281828459045\n")
        report = audit.Audit()
        audit.validate_pooled(report, table, "test")
        self.assertIn("test: independent BH and exponentiation checks", failures(report))


class PooledImpactTests(unittest.TestCase):
    def table(self, altered_p=False, altered_density=False):
        pvalues = [.001, .002, .003, .08, .2, .04, .05, .1, .6]
        if altered_p:
            pvalues[2] = .07
        qvalues = audit.bh_adjust(pvalues)
        text = io.StringIO(newline="")
        writer = csv.writer(text)
        writer.writerow(("Taxa", "Response", "k", "n_MPAs", "model_type", "Estimate", "SE",
                         "CI_lower", "CI_upper", "pval", "pval_fdr", "RR", "RR_lower", "RR_upper"))
        for i, (p, q) in enumerate(zip(pvalues, qvalues)):
            estimate = .2 if altered_density and i == 5 else 0
            writer.writerow((f"taxon{i}", "Biomass" if i < 5 else "Density", 10, 9, "joint",
                             estimate, 1, -1, 1, p, q, math.exp(estimate), math.exp(-1), math.exp(1)))
        return snapshot(text.getvalue())

    def test_density_fdr_can_change_without_density_refit(self):
        old, new = self.table(), self.table(altered_p=True)
        report = audit.Audit()
        audit.compare_pooled(report, old, new)
        self.assertEqual(failures(report), [])
        density = next(row for row in report.details["pooled"]
                       if row["key"] == ("taxon5", "Density"))
        self.assertEqual(set(density["changes"]), {"pval_fdr"})

    def test_changed_density_coefficient_fails(self):
        report = audit.Audit()
        audit.compare_pooled(report, self.table(), self.table(altered_density=True))
        self.assertIn("Table 2: density fits and membership unchanged; q may change", failures(report))

    def test_combined_phase_allows_density_refit(self):
        report = audit.Audit()
        audit.compare_pooled(report, self.table(), self.table(altered_density=True), step_direction_fixed=True)
        self.assertEqual(failures(report), [])


class RawSummaryTests(unittest.TestCase):
    def test_only_target_summary_and_unchanged_n(self):
        header = "taxon_name,source,resp,N,value,sd,se,ci\n"
        target = "Semicossyphus pulcher,LTER,Bio,90,8,16,1.7,3.4\n"
        other = "Semicossyphus pulcher,LTER,Den,90,1,2,0.2,0.4\n"
        before = snapshot(header+target+other)
        report = audit.Audit()
        audit.compare_raw_summary(report, before, snapshot(header+target.replace(",8,16,", ",9,19,")+other))
        self.assertEqual(failures(report), [])
        bad = audit.Audit()
        audit.compare_raw_summary(bad, before, snapshot(header+target.replace(",90,", ",91,")+other))
        self.assertIn("Table 1: only LTER sheephead Bio summaries changed; N unchanged", failures(bad))


class OffsetFidelityTests(unittest.TestCase):
    RAW_HEADER = "CA_MPA_Name_Short,year,taxon_name,source,resp,status,value\n"
    RR_HEADER = "CA_MPA_Name_Short,year,y,source,resp,lnDiff\n"

    def test_adaptive_ratio_matches_normalized_offset_without_scale_dependence(self):
        raw = snapshot(self.RAW_HEADER + "A,2000,taxon,source,Bio,mpa,0\n"
                       "A,2000,taxon,source,Bio,reference,2\n"
                       "A,2001,taxon,source,Bio,mpa,4\n"
                       "A,2001,taxon,source,Bio,reference,2\n")
        ratios = snapshot(self.RR_HEADER + f"A,2000,taxon,source,Bio,{math.log(1/3)}\n"
                          f"A,2001,taxon,source,Bio,{math.log(5/3)}\n")
        self.assertEqual(audit.offset_fidelity(raw, ratios), (2, 2))
        scaled = snapshot(self.RAW_HEADER + "A,2000,taxon,source,Bio,mpa,0\n"
                          "A,2000,taxon,source,Bio,reference,20\n"
                          "A,2001,taxon,source,Bio,mpa,40\n"
                          "A,2001,taxon,source,Bio,reference,20\n")
        self.assertEqual(audit.offset_fidelity(scaled, ratios), (2, 2))

    def test_unmatched_and_all_zero_groups_are_not_compared(self):
        raw = snapshot(self.RAW_HEADER + "A,2000,taxon,source,Bio,mpa,0\n"
                       "A,2000,taxon,source,Bio,reference,0\n"
                       "B,2000,taxon,source,Bio,mpa,1\n"
                       "B,2000,taxon,source,Bio,reference,2\n")
        ratios = snapshot(self.RR_HEADER + "A,2000,taxon,source,Bio,0\n")
        self.assertEqual(audit.offset_fidelity(raw, ratios), (0, 0))

    def test_duplicate_status_fails_instead_of_overwriting(self):
        raw = snapshot(self.RAW_HEADER + "A,2000,taxon,source,Bio,mpa,1\n"
                       "A,2000,taxon,source,Bio,mpa,2\n")
        with self.assertRaisesRegex(ValueError, "Duplicate raw status"):
            audit.offset_fidelity(raw, snapshot(self.RR_HEADER))


class StabilityTests(unittest.TestCase):
    def test_frozen_result_paths_do_not_replace_inputs_or_original_baseline(self):
        root = Path("/tmp/audit-fixture").resolve()
        frozen = Path("/tmp/frozen-fixture").resolve()
        report = audit.Audit(root, frozen)
        self.assertEqual(report.resolve(root / "data/sumstats_final.csv"), frozen / "sumstats_final.csv")
        self.assertEqual(report.resolve(root / "tables/table_02_meta_analysis.csv"),
                         frozen / "tables/table_02_meta_analysis.csv")
        for path in (root / "code/R/08_effect_sizes.R", root / "data/harmonized/example.csv",
                     root / "local_archive/before/tables/table_02_meta_analysis.csv"):
            self.assertEqual(report.resolve(path), path)

    def test_changed_input_fails_final_stability_check(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.csv"
            path.write_bytes(b"before")
            report = audit.Audit()
            report.read_bytes(path)
            path.write_bytes(b"after")
            report.verify_stable()
            self.assertEqual(failures(report), ["All inspected files stable throughout audit"])

    def test_repeated_read_detects_race_immediately(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.csv"
            path.write_bytes(b"before")
            report = audit.Audit()
            report.read_bytes(path)
            path.write_bytes(b"after")
            with self.assertRaisesRegex(ValueError, "changed during audit"):
                report.read_bytes(path)


if __name__ == "__main__":
    unittest.main()
