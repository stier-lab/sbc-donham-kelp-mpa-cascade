#!/usr/bin/env python3
"""Verify the current SBC 127.5 calibration and the historical release chain; never write analysis outputs."""

import argparse
import csv
import hashlib
import io
import json
import math
import statistics
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = "knb-lter-sbc.112.15"
ENTITY = "6bc968141edf7b9e0a6bac923f5f78b2"
TABLE_MD5 = "50eda5d3664bf70780b3bead9674fac1"
METADATA_SHA256 = "b2369950e40c82d0672baef1aa2c0d34afb99c52c8694ae66873295713e02f32"
FIELDS = {"model", "MONTH", "Independent_variable", "Dependent_variable",
          "slope", "intercept", "pvalue", "rsquare", "N"}
MONTH_NAMES = ("January", "February", "March", "April", "May", "June",
               "July", "August", "September", "October", "November", "December")


def biomass_calibration(rows):
    """Select total dry foliar standing crop, not canopy biomass or NPP."""
    result = {}
    for row in rows:
        if set(row) != FIELDS or any(value is None for value in row.values()):
            raise ValueError("Unexpected calibration columns or missing values")
        if row["model"] != "1":
            continue
        month = int(row["MONTH"])
        if not 1 <= month <= 12 or month in result:
            raise ValueError("Invalid or duplicate biomass month")
        name = MONTH_NAMES[month - 1]
        if row["Independent_variable"] != f"Frond density (num m-2) in {name}":
            raise ValueError("Wrong independent variable or month label")
        if row["Dependent_variable"] != f"FSC dry (kg dry mass m-2) in {name}":
            raise ValueError("Wrong biomass quantity, units or month label")
        slope, intercept = float(row["slope"]), float(row["intercept"])
        if not math.isfinite(slope) or slope <= 0 or intercept != 0:
            raise ValueError("Expected finite positive slopes and zero intercepts")
        result[month] = slope
    if set(result) != set(range(1, 13)):
        raise ValueError("Expected exactly twelve monthly biomass calibrations")
    return dict(sorted(result.items()))


def read_calibration(directory):
    metadata = (directory / "metadata.xml").read_bytes()
    table = (directory / "table5.csv").read_bytes()
    if hashlib.sha256(metadata).hexdigest() != METADATA_SHA256:
        raise ValueError("Metadata differs from the verified release snapshot")
    root = ET.fromstring(metadata)
    if root.attrib.get("packageId") != PACKAGE:
        raise ValueError("Wrong package revision")
    tables = root.findall("./dataset/dataTable")
    if len(tables) != 6 or tables[4].attrib.get("id") != ENTITY:
        raise ValueError("Wrong fifth data entity")
    expected_md5 = tables[4].findtext("./physical/authentication")
    if expected_md5 != TABLE_MD5 or hashlib.md5(table).hexdigest() != expected_md5:
        raise ValueError("Table bytes do not match the provider MD5")
    rows = list(csv.DictReader(io.StringIO(table.decode("utf-8"), newline="")))
    if len(rows) != int(tables[4].findtext("numberOfRecords")) or len(rows) != 108:
        raise ValueError("Table record count does not match the release")
    return biomass_calibration(rows)


def read_fixed_calibration(directory):
    metadata = (directory / "sbc1275_eml.xml").read_bytes()
    if hashlib.sha256(metadata).hexdigest() != "a128209f4fd5acc6c359dc1e7d0227ea3b9183c77dcb9dc766b36f5c2ddfcbc4":
        raise ValueError("SBC127.5 metadata differs from the verified snapshot")
    root = ET.fromstring(metadata)
    if root.attrib.get("packageId") != "knb-lter-sbc.127.5":
        raise ValueError("Wrong fixed-coefficient package")
    name = "Macrocystis_Biomass_Relationships_20261016.csv"
    entities = [x for x in root.findall("./dataset/dataTable") if x.findtext("physical/objectName") == name]
    if len(entities) != 1:
        raise ValueError("Missing or ambiguous Macrocystis table")
    table = (directory / name).read_bytes()
    if entities[0].findtext("physical/authentication") != "8ddaf77458d21adfb6d0398cc1c51fda" or hashlib.md5(table).hexdigest() != "8ddaf77458d21adfb6d0398cc1c51fda":
        raise ValueError("Fixed table does not match provider MD5")
    rows = list(csv.DictReader(io.StringIO(table.decode("utf-8-sig"))))
    if len(rows) != 108 or int(entities[0].findtext("numberOfRecords")) != 108:
        raise ValueError("Wrong fixed-table record count")
    for row in rows:
        if row.pop("SP_CODE") != "MAPY" or row.pop("SCIENTIFIC_NAME") != "Macrocystis pyrifera":
            raise ValueError("Wrong kelp species")
        row.pop("Equation")
        row["Independent_variable"] = row.pop("x_variable")
        row["Dependent_variable"] = row.pop("y_variable")
    return biomass_calibration(rows)


def production_calibration(data_repo):
    expression = (
        'e <- new.env(); sys.source(commandArgs(TRUE)[1], envir=e); '
        'cat(jsonlite::toJSON(list(slopes=e$bio_macro(rep(1,12),1:12)/1000),'
        'auto_unbox=TRUE,digits=16),"\\n")')
    result = subprocess.run(["Rscript", "--vanilla", "-e", expression,
                             str(data_repo / "code/R/01_utils.R")],
                            check=True, capture_output=True, text=True)
    return json.loads(result.stdout.splitlines()[-1])["slopes"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path,
                        default=ROOT / "data/external/sbc1275")
    parser.add_argument("--data-repo", type=Path,
                        default=ROOT.parent / "Donham-Stier-CA-MPA-Data-2026")
    parser.add_argument("--strict", action="store_true",
                        help="Require the adopted monthly-calibration input release")
    args = parser.parse_args()
    calibration = read_fixed_calibration(args.source_dir)
    production = production_calibration(args.data_repo)
    if len(production) != 12 or any(not math.isclose(a, b, rel_tol=0, abs_tol=1e-12)
            for a, b in zip(production, calibration.values())):
        raise ValueError("Production conversion does not reproduce all monthly slopes")
    fixed = json.loads((ROOT / "data/harmonized/kelp_fixed_calibration.json").read_text())
    if (fixed.get("package") != "knb-lter-sbc.127.5" or fixed.get("model") != 1
            or fixed.get("calibration_md5") != "8ddaf77458d21adfb6d0398cc1c51fda"
            or fixed.get("applied") is not True):
        raise ValueError("Current fixed-coefficient release is not applied")
    marker = json.loads((ROOT / "data/harmonized/kelp_monthly_calibration.json").read_text())
    if not marker.get("applied") or marker.get("calibration_md5") != TABLE_MD5:
        raise ValueError("Monthly-calibration harmonized release was not applied")
    followup_path = ROOT / "data/harmonized/sheephead_allometry.json"
    followup = json.loads(followup_path.read_text()) if followup_path.exists() else None
    for name, info in marker["files"].items():
        expected_md5 = info["after_md5"]
        if followup is not None:
            expected_source = {"reference_doi": "10.1098/rspb.2014.1817", "a": 0.0144, "b": 3.04,
                               "length_type": "TL", "length_unit": "cm", "weight_unit": "g"}
            if (followup.get("release") != "hamilton_caselle_2015_sheephead" or not followup.get("applied")
                    or any(followup.get("source", {}).get(k) != v for k, v in expected_source.items())
                    or followup.get("scope") != {"programme": "PISCO", "taxon": "Semicossyphus pulcher", "response": "Bio"}
                    or followup.get("files", {}).get(name, {}).get("before_md5") != expected_md5):
                raise ValueError("Sheephead correction does not descend from the checked monthly release")
            expected_md5 = followup["files"][name]["after_md5"]
        for marker, release in (("lobster_bootstrap_correction.json", "lobster_singleton_sampling_2026_10_18"),
                                ("urchin_bootstrap_correction.json", "urchin_singleton_sampling_2026_10_18"),
                                ("kelp_fixed_calibration.json", "kelp_fixed_sbc1275_2026_10_29")):
            child_path = ROOT / "data/harmonized" / marker
            if child_path.exists():
                child = json.loads(child_path.read_text())
                if (child.get("release") != release or not child.get("applied") or
                        child["files"][name]["before_md5"] != expected_md5):
                    raise ValueError("Broken singleton correction lineage: " + marker)
                expected_md5 = child["files"][name]["after_md5"]
        if hashlib.md5((ROOT / "data/harmonized" / name).read_bytes()).hexdigest() != expected_md5:
            raise ValueError("Harmonized release fingerprint changed: " + name)
    report = {
        "source_status": "PASS",
        "package": "knb-lter-sbc.127.5",
        "doi": "10.6073/pasta/3a6272c195d7de89eb412cf1f0c4d63e",
        "table_md5": "8ddaf77458d21adfb6d0398cc1c51fda",
        "biomass_quantity": "dry foliar standing crop, excluding holdfasts and sporophylls",
        "biomass_units": "kg dry mass m-2",
        "monthly_slopes": calibration,
        "production_monthly_slopes": production,
        "production_status": "PASS",
        "reason": "Functional R conversion and adopted input release match the source",
        "unweighted_monthly_mean_not_a_replacement": statistics.mean(calibration.values()),
        "method": "Apply the survey-month coefficient before temporal aggregation; refit and check results",
        "inputs_changed": True,
    }
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
