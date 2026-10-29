"""Synthetic calibration guards; no raw data, model fits or output writes."""

import importlib.util
from pathlib import Path
import tempfile
import unittest


spec = importlib.util.spec_from_file_location(
    "kelp_calibration", Path(__file__).resolve().parents[1] / "audit_kelp_calibration.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def rows():
    return [{"model": "1", "MONTH": str(i),
             "Independent_variable": f"Frond density (num m-2) in {name}",
             "Dependent_variable": f"FSC dry (kg dry mass m-2) in {name}",
             "slope": str(.07 + i / 1000), "intercept": "0",
             "pvalue": "<.0001", "rsquare": ".9", "N": "40"}
            for i, name in enumerate(audit.MONTH_NAMES, 1)]


class CalibrationTests(unittest.TestCase):
    def test_months_are_sorted_and_other_models_are_not_selected(self):
        source = list(reversed(rows()))
        source.append({**source[0], "model": "2", "Dependent_variable": "Canopy FSC dry"})
        self.assertEqual(list(audit.biomass_calibration(source)), list(range(1, 13)))

    def test_missing_and_duplicate_months_fail(self):
        for source in (rows()[:-1], rows() + [rows()[0]],
                       [{**rows()[0], "MONTH": "13"}] + rows()[1:]):
            with self.subTest(source=source), self.assertRaises(ValueError):
                audit.biomass_calibration(source)

    def test_wrong_units_canopy_npp_or_month_label_fail(self):
        for variable in ("FSC wet (kg wet mass m-2) in January", "Canopy FSC dry",
                         "Seasonal NPP dry", "FSC dry (kg dry mass m-2) in February"):
            source = rows()
            source[0]["Dependent_variable"] = variable
            with self.subTest(variable=variable), self.assertRaises(ValueError):
                audit.biomass_calibration(source)

    def test_nonfinite_negative_and_nonzero_intercept_fail(self):
        for field, value in (("slope", "nan"), ("slope", "inf"), ("slope", "-1"),
                             ("slope", "0"), ("intercept", ".01")):
            source = rows()
            source[0][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                audit.biomass_calibration(source)

    def test_wrong_predictor_or_missing_column_fails(self):
        for variable in ("Plant density", "Frond density (num m-2) in February"):
            source = rows()
            source[0]["Independent_variable"] = variable
            with self.subTest(variable=variable), self.assertRaises(ValueError):
                audit.biomass_calibration(source)
        source = rows()
        del source[0]["slope"]
        with self.assertRaises(ValueError):
            audit.biomass_calibration(source)

    def test_source_tampering_fails_before_using_coefficients(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "metadata.xml").write_bytes(b"<eml/>")
            (path / "table5.csv").write_bytes(b"model,MONTH\n1,1\n")
            with self.assertRaisesRegex(ValueError, "Metadata differs"):
                audit.read_calibration(path)


if __name__ == "__main__":
    unittest.main()
