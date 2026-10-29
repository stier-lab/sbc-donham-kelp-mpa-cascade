#!/usr/bin/env python3
"""Read-only audit of the authorized LTER sheephead COUNT correction.

Run after run_all.R AND 25_prereview_sensitivities.R finish. This script never
sources pipeline modules or writes outputs. --pretrend-only evaluates just the
two existing before-trend helpers in isolated R environments; it is safe while
the parent refit runs. --figures checks metadata, not rendered visual quality.
"""

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import sys


PRIMARY_KEY = ("MPA", "Taxa", "Source", "Resp")
CAMPUS_KEY = ("Campus Point SMCA", "S. pulcher", "LTER", "Bio")
NAPLES_KEY = ("Naples SMCA", "S. pulcher", "LTER", "Bio")
STEP_KEYS = {
    ("Gull Island SMR", "M. franciscanus", "KFM", "Den"),
    ("Naples SMCA", "M. franciscanus", "LTER", "Den"),
    ("Santa Barbara Island SMR", "S. purpuratus", "KFM", "Den"),
    ("Scorpion SMR", "S. pulcher", "KFM", "Den"),
    ("Scorpion SMR", "M. franciscanus", "KFM", "Den"),
    ("Scorpion SMR", "M. franciscanus", "KFM", "Bio"),
}
TARGET_MPAS = {CAMPUS_KEY[0], NAPLES_KEY[0]}
MISSING = {"", "NA", "NaN", "nan"}
PRIMARY_NUMBERS = {"Mean", "SE", "SD", "CI", "DW_stat", "DW_pval"}
PREREVIEW_TABLES = (
    "ar1_inflation", "crosstaxa_permutation", "dw_distribution",
    "offset_coverage", "offset_sensitivity", "point_conception",
)


@dataclass
class CsvSnapshot:
    path: Path
    fields: tuple
    rows: list
    records: list
    header: bytes
    data: bytes

    @classmethod
    def parse(cls, path, data):
        # csv.reader.line_num counts physical lines, including quoted newlines.
        # Keep their exact bytes so numeric equality cannot mask other edits.
        lines = data.splitlines(keepends=True)
        reader = csv.reader(io.StringIO(data.decode("utf-8-sig"), newline=""),
                            strict=True)
        fields = tuple(next(reader))
        if not fields or len(set(fields)) != len(fields):
            raise ValueError(f"{path}: empty or duplicate column names")
        previous = reader.line_num
        header = b"".join(lines[:previous])
        rows, records = [], []
        for values in reader:
            end = reader.line_num
            record = b"".join(lines[previous:end])
            previous = end
            if len(values) != len(fields):
                raise ValueError(f"{path}: malformed record ending at line {end}")
            rows.append(dict(zip(fields, values)))
            records.append(record)
        return cls(Path(path), fields, rows, records, header, data)

    def indexed(self, columns):
        absent = set(columns) - set(self.fields)
        if absent:
            raise ValueError(f"{self.path}: missing key columns {sorted(absent)}")
        result = {}
        for row, record in zip(self.rows, self.records):
            key = tuple(row[column] for column in columns)
            if key in result:
                raise ValueError(f"{self.path}: duplicate key {key}")
            result[key] = (row, record)
        return result


class Audit:
    def __init__(self, root=None, results_snapshot=None):
        self.checks = []
        self.details = {}
        self.snapshots = {}
        self.root = Path(root).resolve() if root else None
        self.results_snapshot = Path(results_snapshot).resolve() if results_snapshot else None

    def resolve(self, path):
        path = Path(path).resolve()
        if not self.results_snapshot or not path.is_relative_to(self.root):
            return path
        relative = path.relative_to(self.root)
        if str(relative) == "data/sumstats_final.csv":
            return self.results_snapshot / "sumstats_final.csv"
        if str(relative) == "data/cache/figures_snapshot.rds":
            return self.results_snapshot / "figures_snapshot.rds"
        if relative.parts[0] in {"tables", "plots", "outputs"}:
            return self.results_snapshot / relative
        return path

    def check(self, name, passed, detail=None, warning=False):
        self.checks.append({"name": name,
                            "status": "PASS" if passed else
                            ("REVIEW" if warning else "FAIL"),
                            "detail": detail})

    def read_bytes(self, path):
        path = self.resolve(path)
        data = path.read_bytes()
        if path in self.snapshots and self.snapshots[path] != data:
            raise ValueError(f"Input changed during audit: {path}")
        self.snapshots[path] = data
        return data

    def read_csv(self, path):
        return CsvSnapshot.parse(self.resolve(path), self.read_bytes(path))

    def verify_stable(self):
        changed = [str(path) for path, data in self.snapshots.items()
                   if not path.exists() or path.read_bytes() != data]
        self.check("All inspected files stable throughout audit", not changed, changed)

    def result(self):
        return {"checks": self.checks, "details": self.details,
                "input_sha256": {str(path): hashlib.sha256(data).hexdigest()
                                 for path, data in self.snapshots.items()},
                "limitations": [
                    "No refits, pipeline modules, or manuscript edits performed.",
                    "Byte preservation covers exported harmonized raw/RR records, "
                    "not a reconstruction from fish observations.",
                    "Figure metadata does not establish visual correctness.",
                    "The archived model-selection table may omit primary metadata.",
                ]}


def same_value(left, right, tolerance=1e-10):
    if left == right:
        return True
    if left in MISSING or right in MISSING:
        return left in MISSING and right in MISSING
    try:
        a, b = float(left), float(right)
    except ValueError:
        return False
    return (math.isfinite(a) and math.isfinite(b)
            and math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance))


def changes(left, right, fields):
    return {column: {"before": left[column], "after": right[column]}
            for column in sorted(fields)
            if not same_value(left[column], right[column])}


def target_input(row):
    return (row["CA_MPA_Name_Short"] in TARGET_MPAS
            and row["source"] == "LTER" and row["resp"] == "Bio"
            and row.get("taxon_name", row.get("y")) == "Semicossyphus pulcher")


def compare_inputs(audit, before, after, keys, allowed_fields, name):
    old, new = before.indexed(keys), after.indexed(keys)
    audit.check(f"{name}: exact header preservation", before.header == after.header)
    audit.check(f"{name}: identical record membership", old.keys() == new.keys(),
                {"removed": sorted(old.keys() - new.keys()),
                 "added": sorted(new.keys() - old.keys())})
    byte_changes, forbidden, numerical_changes = [], [], []
    old_non_target = [key for key in old if not target_input(old[key][0])]
    new_non_target = [key for key in new if not target_input(new[key][0])]
    audit.check(f"{name}: non-target record order preserved",
                old_non_target == new_non_target)
    for key in old.keys() & new.keys():
        row_old, record_old = old[key]
        row_new, record_new = new[key]
        delta = changes(row_old, row_new, set(before.fields) & set(after.fields))
        if not target_input(row_old):
            if record_old != record_new:
                byte_changes.append(key)
        else:
            unexpected = set(delta) - set(allowed_fields)
            if unexpected:
                forbidden.append({"key": key, "fields": sorted(unexpected)})
            if delta:
                numerical_changes.append({"key": key, "changes": delta})
    audit.check(f"{name}: all non-target records byte-identical",
                not byte_changes, byte_changes)
    audit.check(f"{name}: target metadata preserved", not forbidden, forbidden)
    audit.check(f"{name}: corrected target values present", bool(numerical_changes))
    audit.details[name] = {"before_rows": len(old), "after_rows": len(new),
                           "non_target_rows_checked": len(old_non_target),
                           "changed_target_records": sorted(numerical_changes,
                                                            key=lambda x: x["key"])}


def compare_primary(audit, before, after, step_direction_fixed=False):
    old, new = before.indexed(PRIMARY_KEY), after.indexed(PRIMARY_KEY)
    audit.check("Primary: identical membership", old.keys() == new.keys(),
                {"removed": sorted(old.keys() - new.keys()),
                 "added": sorted(new.keys() - old.keys())})
    shared = set(before.fields) & set(after.fields)
    required = set(PRIMARY_KEY) | {"Mean", "SE", "CI", "N", "Model", "AnalysisType"}
    audit.check("Primary: full-precision required columns available", required <= shared)
    unexpected, all_changes = [], []
    for key in sorted(old.keys() & new.keys()):
        delta = changes(old[key][0], new[key][0], shared)
        if not delta:
            continue
        item = {"key": key, "changes": delta}
        all_changes.append(item)
        permitted = key == CAMPUS_KEY and not set(delta) - PRIMARY_NUMBERS
        if step_direction_fixed and key in STEP_KEYS:
            permitted = (set(delta) == {"Mean"}
                         and same_value(str(-float(old[key][0]["Mean"])), new[key][0]["Mean"]))
        if not permitted:
            unexpected.append(item)
    audit.check("Primary: only authorized Campus/six-Step changes" if step_direction_fixed else
                "Primary: only Campus LTER sheephead Bio numerics changed",
                not unexpected, unexpected)
    if step_direction_fixed:
        actual_steps = {key for key, row in old.items() if row[0]["Model"] == "Step"}
        valid = {item["key"] for item in all_changes if item["key"] in STEP_KEYS
                 and set(item["changes"]) == {"Mean"}
                 and same_value(str(-float(old[item["key"]][0]["Mean"])),
                                new[item["key"]][0]["Mean"])}
        audit.check("Primary: all and only six stored Step means reversed",
                    actual_steps == valid == STEP_KEYS)
    audit.check("Primary: Campus estimate/uncertainty refreshed",
                any(item["key"] == CAMPUS_KEY and
                    set(item["changes"]) & {"Mean", "SE", "CI"}
                    for item in all_changes))
    audit.check("Primary: Campus retained; Naples still excluded",
                CAMPUS_KEY in new and NAPLES_KEY not in new)
    counts = {stage: {"total": len(index),
                      "Model": dict(Counter(row[0]["Model"] for row in index.values())),
                      "AnalysisType": dict(Counter(row[0]["AnalysisType"]
                                                   for row in index.values()))}
              for stage, index in (("before", old), ("after", new))}
    audit.check("Primary: model/design counts unchanged", counts["before"] == counts["after"])
    audit.details["primary"] = {"baseline_file": str(before.path),
                                 "compared_fields": sorted(shared),
                                 "unavailable_baseline_fields": sorted(set(after.fields)-shared),
                                 "counts": counts, "changes": all_changes}
    return new


def compare_count_only_primary(audit, count_only, after):
    old, new = count_only.indexed(PRIMARY_KEY), after.indexed(PRIMARY_KEY)
    audit.check("Phase 2 vs COUNT-only: exact primary membership", old.keys() == new.keys())
    shared = set(count_only.fields) & set(after.fields)
    unexpected, altered = [], []
    for key in sorted(old.keys() & new.keys()):
        delta = changes(old[key][0], new[key][0], shared)
        if delta:
            altered.append({"key": key, "changes": delta})
            if (key not in STEP_KEYS or set(delta) != {"Mean"}
                    or not same_value(str(-float(old[key][0]["Mean"])), new[key][0]["Mean"])):
                unexpected.append(key)
    audit.check("Phase 2 vs COUNT-only: six Mean sign reversals, all other fields unchanged",
                not unexpected and {row["key"] for row in altered} == STEP_KEYS,
                sorted(unexpected))
    audit.details["primary_vs_count_only"] = altered


def compare_raw_summary(audit, before, after):
    keys = ("taxon_name", "source", "resp")
    old, new = before.indexed(keys), after.indexed(keys)
    audit.check("Table 1: identical taxon/source/response membership", old.keys() == new.keys())
    impacts, unexpected = [], []
    for key in sorted(old.keys() & new.keys()):
        delta = changes(old[key][0], new[key][0], set(before.fields) & set(after.fields))
        if delta:
            impacts.append({"key": key, "changes": delta})
        if delta and (key != ("Semicossyphus pulcher", "LTER", "Bio") or "N" in delta):
            unexpected.append({"key": key, "changes": delta})
    audit.check("Table 1: only LTER sheephead Bio summaries changed; N unchanged",
                not unexpected, unexpected)
    audit.details["raw_summary"] = impacts


def bh_adjust(pvalues):
    if any(not math.isfinite(p) or not 0 <= p <= 1 for p in pvalues):
        raise ValueError("BH adjustment requires finite P values in [0,1]")
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    adjusted, minimum = [0.0] * len(pvalues), 1.0
    for rank in range(len(order), 0, -1):
        index = order[rank - 1]
        minimum = min(minimum, pvalues[index] * len(order) / rank)
        adjusted[index] = minimum
    return adjusted


def validate_pooled(audit, table, label):
    pvalues = [float(row["pval"]) for row in table.rows]
    adjusted = bh_adjust(pvalues)
    audit.check(f"{label}: nine-test FDR family", len(pvalues) == 9)
    errors = []
    for row, q in zip(table.rows, adjusted):
        key = (row["Taxa"], row["Response"])
        if not same_value(str(q), row["pval_fdr"]):
            errors.append({"key": key, "issue": "BH q", "expected": q})
        estimate, lower, upper, se = [float(row[field]) for field in
                                      ("Estimate", "CI_lower", "CI_upper", "SE")]
        if (not all(math.isfinite(x) for x in (estimate, lower, upper, se))
                or not lower <= estimate <= upper or se <= 0):
            errors.append({"key": key, "issue": "nonfinite/invalid log estimate or CI"})
        for source, target in (("Estimate", "RR"), ("CI_lower", "RR_lower"),
                               ("CI_upper", "RR_upper")):
            if not same_value(str(math.exp(float(row[source]))), row[target]):
                errors.append({"key": key, "issue": f"exp({source}) != {target}"})
    audit.check(f"{label}: independent BH and exponentiation checks", not errors, errors)


def compare_pooled(audit, before, after, step_direction_fixed=False):
    old = before.indexed(("Taxa", "Response"))
    new = after.indexed(("Taxa", "Response"))
    audit.check("Table 2: identical taxa/response membership", old.keys() == new.keys())
    validate_pooled(audit, before, "Baseline Table 2")
    validate_pooled(audit, after, "Current Table 2")
    unexpected, impacts = [], []
    for key in sorted(old.keys() & new.keys()):
        left, right = old[key][0], new[key][0]
        delta = changes(left, right, set(before.fields) & set(after.fields))
        forbidden = {"Taxa", "Response", "k", "n_MPAs", "model_type"} & set(delta)
        if key[1] != "Biomass" and not step_direction_fixed:
            forbidden |= set(delta) - {"pval_fdr"}
        if forbidden:
            unexpected.append({"key": key, "fields": sorted(forbidden)})
        impacts.append({"key": key, "changes": delta,
                        "nominal_support_before": float(left["pval"]) < .05,
                        "nominal_support_after": float(right["pval"]) < .05,
                        "FDR_support_before": float(left["pval_fdr"]) < .05,
                        "FDR_support_after": float(right["pval_fdr"]) < .05})
    audit.check("Table 2: combined refit membership/design unchanged" if step_direction_fixed else
                "Table 2: density fits and membership unchanged; q may change",
                not unexpected, unexpected)
    audit.details["pooled"] = impacts
    audit.details["pooled_qualitative_counts"] = {
        stage: {"nominal_supported": sum(float(row[0]["pval"]) < .05 for row in index.values()),
                "FDR_supported": sum(float(row[0]["pval_fdr"]) < .05 for row in index.values()),
                "FDR_positive_biomass": sum(row[0]["Response"] == "Biomass"
                                             and float(row[0]["Estimate"]) > 0
                                             and float(row[0]["pval_fdr"]) < .05
                                             for row in index.values()),
                "FDR_positive_density": sum(row[0]["Response"] == "Density"
                                             and float(row[0]["Estimate"]) > 0
                                             and float(row[0]["pval_fdr"]) < .05
                                             for row in index.values()),
                "FDR_negative_density": sum(row[0]["Response"] == "Density"
                                             and float(row[0]["Estimate"]) < 0
                                             and float(row[0]["pval_fdr"]) < .05
                                             for row in index.values())}
        for stage, index in (("before", old), ("after", new))}


def compare_pathways(audit, before, after, step_direction_fixed=False):
    keys = ("Predictor", "Response_Var")
    old, new = before.indexed(keys), after.indexed(keys)
    audit.check("Table 3: identical pathway membership", old.keys() == new.keys())
    unexpected, impacts = [], []
    affected = ("S. pulcher biomass", "S. purpuratus biomass")
    for key in sorted(old.keys() & new.keys()):
        delta = changes(old[key][0], new[key][0], set(before.fields) & set(after.fields))
        if delta:
            impacts.append({"key": key, "changes": delta})
        if delta and ((key != affected and not step_direction_fixed)
                      or set(delta) & {"k", "Relationship", "Predictor", "Response_Var"}):
            unexpected.append({"key": key, "changes": delta})
    audit.check("Table 3: combined refit pathway membership/k unchanged" if step_direction_fixed else
                "Table 3: only direct sheephead-biomass pathway may change",
                not unexpected, unexpected)
    audit.details["pathways"] = impacts


def check_coverage(audit, baseline, root, primary):
    atlas_name = "table_s_pair_history_atlas_index.csv"
    before = audit.read_csv(baseline / "tables" / atlas_name)
    after = audit.read_csv(root / "tables" / atlas_name)
    old, new = before.indexed(PRIMARY_KEY), after.indexed(PRIMARY_KEY)
    audit.check("Atlas: identical primary keys", old.keys() == new.keys() == primary.keys())
    coordinates = [(row["page"], row["panel"]) for row in after.rows]
    expected = [(str(page), panel) for page in range(1, math.ceil(len(primary)/4)+1)
                for panel in "ABCD"][:len(primary)]
    audit.check("Atlas: unique consecutive page/panel addresses",
                len(set(coordinates)) == len(coordinates)
                and sorted(coordinates, key=lambda x: (int(x[0]), x[1])) == expected)
    conflicts = []
    for key in primary.keys() & new.keys():
        row, effect = new[key][0], primary[key][0]
        if (row["coverage_status"] != "matched"
                or not same_value(row["final_N"], effect["N"])
                or not same_value(row["n_pairs"], effect["N"])
                or row["Model"] != effect["Model"]
                or row["AnalysisType"] != effect["AnalysisType"]):
            conflicts.append(key)
    audit.check("Atlas: histories match primary N/model/design", not conflicts, conflicts)
    footprint_fields = ("n_before", "n_after", "year_min", "year_max", "missing_years",
                        "n_pairs", "Model", "AnalysisType", "page", "panel")
    footprint_changes = [{"key": key, "changes": changes(old[key][0], new[key][0], footprint_fields)}
                         for key in sorted(old.keys() & new.keys())
                         if changes(old[key][0], new[key][0], footprint_fields)]
    audit.check("Atlas: annual footprints and addresses unchanged",
                not footprint_changes, footprint_changes)
    mean_before = sum(row["Model"] == "Mean" and row["AnalysisType"] == "CI"
                      and int(row["n_before"]) > 0 for row in after.rows)
    cov_name = "table_s_design_matched_coverage.csv"
    cov_old = audit.read_csv(baseline / "tables" / cov_name)
    cov_new = audit.read_csv(root / "tables" / cov_name)
    oi, ni = cov_old.indexed(("series_id",)), cov_new.indexed(("series_id",))
    cov_changes = [key for key in oi.keys() & ni.keys()
                   if changes(oi[key][0], ni[key][0], set(cov_old.fields) & set(cov_new.fields))]
    audit.check("Matched-design eligibility and footprints unchanged",
                oi.keys() == ni.keys() and not cov_changes, sorted(cov_changes))
    matched = sum(row["matched_eligible"] == "TRUE" for row in cov_new.rows)
    audit.details["coverage"] = {"selected_histories": len(new),
                                  "original_figure_pages": len(set(row["page"] for row in after.rows)),
                                  "Mean_CI_with_before": mean_before,
                                  "matched_eligible": matched}
    registry = audit.read_csv(root / "data/harmonized/reference_sites.csv")
    registry_keys = audit.read_csv(root / "data/harmonized/reference_primary_keys.csv")
    ri = registry_keys.indexed(PRIMARY_KEY)
    audit.check("Registry: exact primary key membership", ri.keys() == primary.keys())
    series = {row["series_id"] for row in registry.rows}
    audit.check("Registry: all primary series have reference footprints",
                series == {row["series_id"] for row in registry_keys.rows})
    audit.details["coverage"].update(reference_rows=len(registry.rows), reference_series=len(series))


def check_exports(audit, baseline, root, pending_25=False, step_direction_fixed=False):
    identical = []
    for name in ("harmonized_landsat_rr.csv", "harmonized_landsat_annual_index.csv",
                 "harmonized_site_metadata.csv"):
        old, new = baseline / "harmonized" / name, root / "data/harmonized" / name
        if old.exists():
            audit.check(f"Unaffected export byte identity: {name}",
                        audit.read_bytes(old) == audit.read_bytes(new))
            identical.append(name)
    audit.details["unchanged_export_files_checked"] = identical
    # Inventory all archived tables without assuming row order equals model identity.
    inventory, unaffected_checked, unaffected_changed = [], [], []
    for old_path in sorted((baseline / "tables").glob("*.csv")):
        if pending_25 and old_path.name.startswith("table_s_prereview_"):
            continue
        new_path = root / "tables" / old_path.name
        if not audit.resolve(new_path).exists():
            audit.check(f"Archived table still present: {old_path.name}", False)
            continue
        old, new = audit.read_csv(old_path), audit.read_csv(new_path)
        if old_path.name.startswith(("table_heatwave", "table_compound", "table_resilience",
                                     "table_resistance", "table_ecological_memory",
                                     "table_s_historical_heatwave")):
            unaffected_checked.append(old_path.name)
            if old.data != new.data:
                unaffected_changed.append(old_path.name)
        if old.data != new.data:
            exact_rows_equal = (set(old.fields) == set(new.fields)
                                and Counter(tuple(row.get(f) for f in sorted(old.fields))
                                            for row in old.rows)
                                == Counter(tuple(row.get(f) for f in sorted(old.fields))
                                           for row in new.rows))
            inventory.append({"file": old_path.name, "before_rows": len(old.rows),
                              "after_rows": len(new.rows),
                              "only_order_or_serialization_changed": exact_rows_equal})
    audit.details["table_file_changes"] = inventory
    audit.details["unaffected_disturbance_tables_checked"] = unaffected_checked
    if step_direction_fixed:
        historical_changes = [name for name in unaffected_changed
                              if name.startswith("table_s_historical_heatwave")]
        audit.check("Historical-heatwave tables byte-identical", not historical_changes, historical_changes)
        audit.details["disturbance_tables_changed_in_combined_refit"] = unaffected_changed
    else:
        audit.check("Disturbance/resilience/historical tables byte-identical",
                    not unaffected_changed, unaffected_changed)
    if pending_25:
        audit.check("Module 25: excluded while parent runner is active", False,
                    "Incomplete audit: rerun without --pending-25 after module 25 finishes.",
                    warning=True)
        return
    if audit.results_snapshot:
        missing = [suffix for suffix in PREREVIEW_TABLES
                   if not audit.resolve(root / "tables" / f"table_s_prereview_{suffix}.csv").exists()]
        audit.check("Module 25: frozen six-table set present (copy timestamps not used)", not missing, missing)
        audit.details["frozen_module25_limit"] = ("Completion is recorded by the parent; copied-file "
                                                 "timestamps are not a freshness proof.")
        return
    primary_mtime = audit.resolve(root / "data/sumstats_final.csv").stat().st_mtime
    stale = []
    for suffix in PREREVIEW_TABLES:
        path = audit.resolve(root / "tables" / f"table_s_prereview_{suffix}.csv")
        if not path.exists() or path.stat().st_mtime < primary_mtime:
            stale.append(str(path))
    audit.check("Module 25: all six dependent tables refreshed after primary input",
                not stale, stale)


def offset_fidelity(raw, ratios):
    """Reconstruct only the adaptive annual ratios, never effect-size models."""
    minimum, pairs = {}, {}
    for row in raw.rows:
        group = tuple(row[c] for c in ("CA_MPA_Name_Short", "taxon_name", "source", "resp"))
        key = (group[0], row["year"], *group[1:])
        value = float(row["value"])
        if value > 0:
            minimum[group] = min(value, minimum.get(group, value))
        if row["status"] in pairs.setdefault(key, {}):
            raise ValueError(f"Duplicate raw status in offset reconstruction: {key}")
        pairs[key][row["status"]] = value
    index = ratios.indexed(("CA_MPA_Name_Short", "year", "y", "source", "resp"))
    exact, compared = 0, 0
    for key, values in pairs.items():
        group = (key[0], *key[2:])
        if key not in index or group not in minimum or set(values) != {"mpa", "reference"}:
            continue
        offset = minimum[group] / 2
        reconstructed = math.log(values["mpa"] + offset) - math.log(values["reference"] + offset)
        compared += 1
        exact += abs(reconstructed - float(index[key][0]["lnDiff"])) < 1e-8
    return exact, compared


def check_prereview(audit, root):
    table = lambda suffix: audit.read_csv(root / "tables" / f"table_s_prereview_{suffix}.csv")
    primary = audit.read_csv(root / "tables/table_model_selection.csv")
    pathways = audit.read_csv(root / "tables/table_03_cross_taxa_meta_regression.csv").indexed(("Relationship",))
    permutation = table("crosstaxa_permutation")
    perm_index = permutation.indexed(("Relationship",))
    stale = []
    for key, (row, _) in perm_index.items():
        if key not in pathways:
            stale.append(key)
            continue
        published = pathways[key][0]
        if (int(row["k"]) != int(published["k"])
                or abs(float(row["slope_published"]) - float(published["Slope_Est"])) > .000050001
                or abs(float(row["p_parametric_published"]) - float(published["Slope_pval"])) > .000050001):
            stale.append(key)
    audit.check("Module 25: all six published references match current Table 3",
                len(perm_index) == 6 and perm_index.keys() == pathways.keys() and not stale, stale)
    invalid_modes = [row["Relationship"] for row in permutation.rows
                     if (row["permutation_mode"], int(row["n_permutations"])) !=
                     (("exact", math.factorial(int(row["k"]))) if int(row["k"]) <= 7
                      else ("monte-carlo(1e4)", 10000))]
    audit.check("Module 25: permutation counts/modes match six retained k values", not invalid_modes, invalid_modes)
    concordance = [row["Relationship"] for row in permutation.rows if row["reproduces_published"] == "TRUE"]
    audit.details["module25_permutation"] = {
        "rows": len(permutation.rows), "primary_concordant_relations": concordance,
        "nonprimary_p_lt_05": [row for row in permutation.rows
                               if row["reproduces_published"] == "FALSE" and float(row["p_permutation"]) < .05],
        "interpretation": "Grouped reconstructions are not the primary regressions; their P values do not replace Table 3.",
    }
    inflation = table("ar1_inflation")
    pooled = audit.read_csv(root / "tables/table_02_meta_analysis.csv").indexed(("Taxa", "Response"))
    ar1 = inflation.indexed(("Taxa", "Response"))
    audit.check("Module 25: nine AR1 comparisons retain primary taxon/response keys",
                len(ar1) == 9 and ar1.keys() == pooled.keys())
    invalid_flags = [key for key, (row, _) in ar1.items()
                     if ((float(row["est_reported_se"]) * float(row["est_ar1_se"]) < 0) != (row["sign_flip"] == "TRUE")
                         or ((float(row["p_reported"]) < .05) != (float(row["p_ar1"]) < .05)) !=
                         (row["significance_flip"] == "TRUE"))]
    audit.check("Module 25: AR1 sign/nominal flags agree with tabulated estimates/P", not invalid_flags, invalid_flags)
    rho = [max(0, min(.95, 1 - float(row["DW_stat"]) / 2))
           if row["DW_stat"] not in MISSING else 0
           for row in primary.rows if row["Mean"] not in MISSING and row["SE"] not in MISSING and float(row["SE"]) > 0]
    audit.details["module25_ar1"] = {
        "rows": len(ar1), "median_ci_width_ratio": statistics.median(float(row["ci_width_ratio"]) for row in inflation.rows),
        "sign_flips": sum(row["sign_flip"] == "TRUE" for row in inflation.rows),
        "nominal_flips": sum(row["significance_flip"] == "TRUE" for row in inflation.rows),
        "median_rho": statistics.median(rho), "rho_histories": len(rho),
        "interpretation": "Within the separate legacy-grouped fits, not a reproduction of the primary joint model.",
    }
    dw = table("dw_distribution").indexed(("Model",))
    tested = [row for row in primary.rows if row["DW_pval"] not in MISSING]
    invalid_dw = []
    for key, (row, _) in dw.items():
        subset = tested if key == ("ALL",) else [r for r in tested if r["Model"] == key[0]]
        n, flagged = len(subset), sum(float(r["DW_pval"]) < .05 for r in subset)
        if (int(row["n_tested"]) != n or int(row["n_p_lt_05"]) != flagged
                or abs(float(row["pct_p_lt_05"]) - 100 * flagged / n) > .050001
                or abs(float(row["median_DW"]) - statistics.median(float(r["DW_stat"]) for r in subset)) > .000500001):
            invalid_dw.append(key)
    audit.check("Module 25: DW groups/counts/medians match current primary diagnostics",
                dw.keys() == {("ALL",), ("Mean",), ("Linear",), ("Step",)} and not invalid_dw, invalid_dw)
    audit.details["module25_dw"] = {key[0]: row for key, (row, _) in dw.items()}
    offsets = table("offset_sensitivity")
    offsets.indexed(("taxon_name", "resp"))
    exact, compared = offset_fidelity(audit.read_csv(root / "data/harmonized/harmonized_raw_responses.csv"),
                                    audit.read_csv(root / "data/harmonized/harmonized_response_ratios.csv"))
    coverage = dict(csv.reader(io.StringIO(audit.read_bytes(root / "tables/table_s_prereview_offset_coverage.csv").decode())))
    total = sum(int(row["n"]) for row in offsets.rows)
    audit.check("Module 25: ten offset summaries and independent adaptive fidelity agree",
                len(offsets.rows) == 10 and total == exact == int(coverage["rows_compared"])
                and abs(float(coverage["reconstruction_fidelity_exact_fraction"]) - exact / compared) <= .000050001)
    audit.details["module25_offset"] = {"summary_rows": len(offsets.rows), "exact": exact, "compared": compared,
                                        "fraction": exact / compared, "scope": "Annual-ratio shifts, not refitted primary models."}
    pc = table("point_conception")
    kelp = pooled[("M. pyrifera", "Biomass")][0]
    audit.check("Module 25: Point Conception row references current pooled kelp fit",
                len(pc.rows) == 1 and pc.rows[0]["MPA"] == "Point Conception SMR"
                and int(pc.rows[0]["published_kelp_k"]) == int(kelp["k"])
                and abs(float(pc.rows[0]["published_kelp_estimate"]) - float(kelp["Estimate"])) <= .000050001)
    audit.details["module25_point_conception"] = {"rows": pc.rows, "interpretation": "Post-only summary plus main-fit reference; not a joint-model inclusion refit."}


def pretrend_check(audit, baseline, root, rscript):
    """Evaluate only parsed helper definitions, never source the module."""
    code_paths = (baseline / "R/08_effect_sizes.R", root / "code/R/08_effect_sizes.R")
    rr_paths = (baseline / "harmonized/harmonized_response_ratios.csv",
                root / "data/harmonized/harmonized_response_ratios.csv")
    for path in (*code_paths, *rr_paths):
        audit.read_bytes(path)
    rcode = r'''
args <- commandArgs(trailingOnly=TRUE)
for (stage in 1:2) {
  ast <- parse(args[[stage]])
  env <- new.env(parent=globalenv())
  for (expr in ast) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in%
        c("add_time_columns", "test_linear_before")) eval(expr, env)
  }
  stopifnot(is.function(env$add_time_columns), is.function(env$test_linear_before))
  x <- read.csv(args[[stage+2]], stringsAsFactors=FALSE)
  x <- subset(x, source=="LTER" & y=="Semicossyphus pulcher" & resp=="Bio")
  x <- env$add_time_columns(x)
  screen <- env$test_linear_before(x)
  for (name in c("Campus Point SMCA", "Naples SMCA")) {
    d <- x[x$CA_MPA_Name_Short == name, ]
    b <- subset(d, BA=="Before")
    fit <- lm(lnDiff ~ time.true, b)
    z <- summary(fit)$coefficients["time.true", ]
    p <- screen$pSlope[screen$site == name]
    stopifnot(length(p)==1, isTRUE(all.equal(p, z[[4]])))
    cat(paste(c(c("before","after")[[stage]], name, nrow(b), sum(d$BA=="After"),
                paste(b$year,collapse=","), paste(b$time.true,collapse=","),
                sprintf("%.17g",z[[1]]),sprintf("%.17g",z[[2]]),sprintf("%.17g",p)),
              collapse="\t"), "\n", sep="")
  }
}
'''
    result = subprocess.run([rscript, "--vanilla", "-e", rcode,
                             *(str(path) for path in (*code_paths, *rr_paths))],
                            check=True, capture_output=True, text=True, timeout=60)
    rows = []
    for line in result.stdout.splitlines():
        stage, mpa, before_n, after_n, years, indices, slope, se, p = line.split("\t")
        rows.append({"stage": stage, "MPA": mpa, "n_before": int(before_n),
                     "n_after": int(after_n), "before_years": years,
                     "before_time_true": indices, "slope": float(slope),
                     "SE": float(se), "p": float(p)})
    audit.details["before_trend"] = rows
    audit.check("Before-trend helper: both stages and MPAs inspected", len(rows) == 4)
    audit.check("Naples exclusion reproduced by current biomass P<0.05 screen",
                any(row["stage"] == "after" and row["MPA"] == NAPLES_KEY[0]
                    and row["p"] < .05 for row in rows),
                "REVIEW only: hardcoded exclusion predates correction; do not change routing here.",
                warning=True)


def check_figures(audit, baseline, root):
    from PIL import Image
    from pypdf import PdfReader

    names = ("fig_01_mpa_map", "fig_02_cascade_case_studies", "fig_03_mean_effects",
             "fig_04_recovery_curves", "fig_kelp_resilience_paired", "fig_s_primary_pathway")
    metadata = []
    for name in names:
        current = {}
        for suffix in ("png", "pdf"):
            path = root / "plots" / f"{name}.{suffix}"
            data = audit.read_bytes(path)
            old_path = baseline / "plots" / path.name
            current[suffix] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(),
                               "same_as_baseline": old_path.exists()
                               and audit.read_bytes(old_path) == data}
            if suffix == "png":
                with Image.open(io.BytesIO(data)) as im:
                    current[suffix].update(size=list(im.size), dpi=im.info.get("dpi"))
                audit.check(f"Figure metadata {name}: nonempty PNG",
                            min(current[suffix]["size"]) > 0)
            else:
                pdf = PdfReader(io.BytesIO(data))
                current[suffix].update(pages=len(pdf.pages),
                                      page_sizes_pt=[[float(p.mediabox.width),
                                                      float(p.mediabox.height)] for p in pdf.pages])
                audit.check(f"Figure metadata {name}: one-page PDF", len(pdf.pages) == 1)
        metadata.append({"figure": name, **current})
    atlas = audit.read_csv(root / "tables/table_s_pair_history_atlas_index.csv")
    combined = PdfReader(io.BytesIO(audit.read_bytes(root / "plots/fig_s_pair_history_atlas.pdf")))
    audit.check("Atlas metadata: original vector PDF page count matches index",
                len(combined.pages) == len(set(row["page"] for row in atlas.rows)))
    for name in sorted(set(row["pdf"] for row in atlas.rows)):
        pdf = PdfReader(io.BytesIO(audit.read_bytes(root / name)))
        audit.check(f"Atlas metadata: {name} has one PDF page", len(pdf.pages) == 1)
    membership = audit.read_csv(root / "tables/table_s_primary_pathway_membership.csv")
    md5s = {"primary_input_md5": hashlib.md5(audit.read_bytes(root / "data/sumstats_final.csv")).hexdigest(),
            "table3_input_md5": hashlib.md5(audit.read_bytes(root / "tables/table_03_cross_taxa_meta_regression.csv")).hexdigest()}
    audit.check("Pathway figure: input hashes match refreshed primary/Table 3",
                bool(membership.rows) and all(row[k] == value for row in membership.rows
                                             for k, value in md5s.items()), md5s)
    if not audit.results_snapshot:
        legend = audit.read_bytes(root / "docs/primary_pathway_legend.md").decode("utf-8")
        audit.check("Pathway legend: current input hashes recorded",
                    all(value in legend for value in md5s.values()))
    else:
        audit.details["snapshot_legend_limit"] = "Legends were not archived; live legend excluded from frozen replay."
    audit.details["figure_metadata"] = metadata


def check_shape_display(audit, root):
    from PIL import Image
    from pypdf import PdfReader

    primary = audit.read_csv(root / "data/sumstats_final.csv").indexed(PRIMARY_KEY)
    selected = {key: row[0] for key, row in primary.items() if row[0]["AnalysisType"] == "pBACIPS"}
    table_dir = root / "tables"
    fit = audit.read_csv(table_dir / "table_s_pbacips_shape_fit_audit.csv").indexed(PRIMARY_KEY)
    audit.check("Module 32: fit audit covers exactly retained pBACIPS keys", fit.keys() == selected.keys())
    issues = []
    for key in selected.keys() & fit.keys():
        row, effect = fit[key][0], selected[key]
        if (row["status"] != "verified" or row["Model"] != effect["Model"]
                or not same_value(row["N"], effect["N"])
                or not same_value(row["primary_effect"], effect["Mean"])
                or not same_value(row["reconstructed_effect"], effect["Mean"], 1e-6)):
            issues.append(key)
    audit.check("Module 32: verified display effects agree with refreshed primary", not issues, sorted(issues))
    provenance = audit.read_csv(table_dir / "table_s_pbacips_shape_provenance.csv")
    hash_issues = []
    for row in provenance.rows:
        path = (root / row["input"]).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ValueError(f"Module-32 provenance points outside repository: {path}")
        if hashlib.md5(audit.read_bytes(path)).hexdigest() != row["md5"]:
            hash_issues.append(row["input"])
    audit.check("Module 32: eight input/code fingerprints match", len(provenance.rows) == 8 and not hash_issues,
                hash_issues)
    histories = PdfReader(io.BytesIO(audit.read_bytes(root / "plots/fig_s_pbacips_selected_histories.pdf")))
    audit.check("Module 32: companion page count matches retained cohort",
                len(histories.pages) == math.ceil(len(selected)/4))
    nonlinear_n = sum(row["Model"] in {"Sigmoid", "Asymptotic"} for row in selected.values())
    stem = root / "plots/fig_s_pbacips_nonlinear_comparison"
    figure_data = audit.read_bytes(stem.with_suffix(".png"))
    with Image.open(io.BytesIO(figure_data)) as im:
        metadata = {"size": list(im.size), "dpi": im.info.get("dpi")}
    pdf = PdfReader(io.BytesIO(audit.read_bytes(stem.with_suffix(".pdf"))))
    audit.check("Module 32: nonlinear comparison has one PDF page", len(pdf.pages) == 1)
    metadata["pages"] = len(pdf.pages)
    metadata["page_sizes_pt"] = [[float(p.mediabox.width), float(p.mediabox.height)] for p in pdf.pages]
    audit.details["shape_display"] = {"retained_histories": len(selected), "nonlinear_histories": nonlinear_n,
                                      "companion_pages": len(histories.pages), "nonlinear_figure": metadata}


def run_audit(baseline, root, figures=False, pretrend=False, rscript=None, pending_25=False,
              shape_display=False, results_snapshot=None, step_direction_fixed=False,
              count_only=None):
    audit = Audit(root, results_snapshot)
    audit.details["phase"] = ("COUNT + six-Step correction" if step_direction_fixed else "COUNT-only")
    audit.details["results_snapshot"] = str(results_snapshot) if results_snapshot else None
    for name, keys, allowed in (
        ("raw", ("CA_MPA_Name_Short", "year", "taxon_name", "source", "resp", "status"), {"value"}),
        ("ratios", ("CA_MPA_Name_Short", "year", "y", "source", "resp"), {"lnDiff", "mpa", "reference", "Diff"}),
    ):
        filename = f"harmonized_{'raw_responses' if name == 'raw' else 'response_ratios'}.csv"
        compare_inputs(audit, audit.read_csv(baseline / "harmonized" / filename),
                       audit.read_csv(root / "data/harmonized" / filename), keys, allowed, name)
    candidates = (baseline / "sumstats_final.csv", baseline / "data/sumstats_final.csv",
                  baseline / "tables/table_model_selection.csv")
    before_primary = next(path for path in candidates if path.exists())
    after_primary = audit.read_csv(root / "data/sumstats_final.csv")
    primary = compare_primary(audit, audit.read_csv(before_primary), after_primary, step_direction_fixed)
    if step_direction_fixed:
        count_only = count_only or root / "local_archive/lter_count_refit_2026-10-10/count_only"
        compare_count_only_primary(audit, audit.read_csv(count_only / "sumstats_final.csv"), after_primary)
    compare_raw_summary(audit, audit.read_csv(baseline / "tables/average_responses.csv"),
                        audit.read_csv(root / "tables/average_responses.csv"))
    for filename, compare in (("table_02_meta_analysis.csv", compare_pooled),
                              ("table_03_cross_taxa_meta_regression.csv", compare_pathways)):
        compare(audit, audit.read_csv(baseline / "tables" / filename),
                audit.read_csv(root / "tables" / filename), step_direction_fixed)
    check_coverage(audit, baseline, root, primary)
    check_exports(audit, baseline, root, pending_25, step_direction_fixed)
    if not pending_25:
        check_prereview(audit, root)
    old_code = audit.read_bytes(baseline / "R/08_effect_sizes.R").decode("utf-8")
    new_code = audit.read_bytes(root / "code/R/08_effect_sizes.R").decode("utf-8")
    if not results_snapshot:
        if step_direction_fixed:
            def outside_step(code):
                prefix, rest = code.replace("\r\n", "\n").split("add_step_effect_size <- function", 1)
                _, suffix = rest.split("#' Add a single linear-model", 1)
                return prefix, suffix
            if rscript:
                paths = (baseline / "R/08_effect_sizes.R", root / "code/R/08_effect_sizes.R",
                         baseline / "R/00c_analysis_constants.R", root / "code/R/00c_analysis_constants.R")
                for path in paths:
                    audit.read_bytes(path)
                rcode = r'''
args <- commandArgs(trailingOnly=TRUE)
outside_step <- function(path) {
  x <- as.list(parse(path, keep.source=FALSE))
  Filter(function(e) !(is.call(e) && identical(e[[1]], as.name("<-")) &&
                      identical(e[[2]], as.name("add_step_effect_size"))), x)
}
cat(identical(outside_step(args[[1]]), outside_step(args[[2]])), "\n")
cat(identical(parse(args[[3]], keep.source=FALSE), parse(args[[4]], keep.source=FALSE)), "\n")
'''
                parsed = subprocess.run([rscript, "--vanilla", "-e", rcode, *map(str, paths)],
                                        check=True, capture_output=True, text=True).stdout.split()
                audit.check("Effect-size executable AST unchanged outside authorized Step helper", parsed[0] == "TRUE")
                audit.check("Analysis constants executable AST unchanged", parsed[1] == "TRUE")
            else:
                audit.check("Effect-size module text unchanged outside authorized Step helper",
                            outside_step(old_code) == outside_step(new_code),
                            "Rscript can distinguish comment-only changes by parsed AST comparison.", warning=True)
        else:
            audit.check("Effect-size module byte-identical to archived rules", old_code == new_code,
                        "If changed, inspect the diff; this check does not infer a methodological change.",
                        warning=True)
    else:
        audit.details["frozen_source_rule"] = ("Frozen results were audited before the Step edit. "
                                                "Live Step code is not attributed to this earlier phase.")
    if figures:
        check_figures(audit, baseline, root)
    if shape_display:
        check_shape_display(audit, root)
    if pretrend:
        pretrend_check(audit, baseline, root, rscript)
    audit.verify_stable()
    return audit.result()


def main():
    root_default = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=root_default)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--figures", action="store_true")
    parser.add_argument("--shape-display", action="store_true",
                        help="Check finalized module-32 display provenance/metadata without fits")
    parser.add_argument("--pretrend", action="store_true")
    parser.add_argument("--pretrend-only", action="store_true")
    parser.add_argument("--pending-25", action="store_true",
                        help="Exclude module-25 tables while parent runs them; not a final pass")
    parser.add_argument("--results-snapshot", type=Path,
                        help="Replay frozen results/cache/plots with live corrected harmonized inputs")
    parser.add_argument("--step-direction-fixed", action="store_true",
                        help="Explicit phase-2 gate: six Step sign reversals plus Campus vs original")
    parser.add_argument("--count-only", type=Path, help="Preserved phase-1 snapshot for phase-2 comparison")
    parser.add_argument("--rscript", default=shutil.which("Rscript"))
    parser.add_argument("--json", action="store_true", help="Print full report to stdout, never a file")
    args = parser.parse_args()
    root = args.root.resolve()
    baseline = (args.baseline or root / "local_archive/lter_count_refit_2026-10-10/before").resolve()
    if (args.pretrend or args.pretrend_only) and not args.rscript:
        parser.error("Rscript is required only for the optional before-trend screen")
    try:
        if args.pretrend_only:
            audit = Audit()
            pretrend_check(audit, baseline, root, args.rscript)
            audit.verify_stable()
            report = audit.result()
        else:
            report = run_audit(baseline, root, args.figures, args.pretrend, args.rscript,
                               args.pending_25, args.shape_display, args.results_snapshot,
                               args.step_direction_fixed, args.count_only)
    except (OSError, ValueError, StopIteration, csv.Error, subprocess.SubprocessError,
            ImportError, KeyError) as error:
        print(f"Audit incomplete: {error}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(report, indent=2, allow_nan=False))
    else:
        print("LTER COUNT refit audit (read-only)")
        counts = Counter(item["status"] for item in report["checks"])
        print("Checks: " + ", ".join(f"{key}={value}" for key, value in sorted(counts.items())))
        for item in report["checks"]:
            if item["status"] != "PASS":
                print(f"{item['status']}: {item['name']}: {item['detail']}")
        print(json.dumps(report["details"], indent=2, allow_nan=False))
    return int(any(item["status"] == "FAIL" for item in report["checks"]))


if __name__ == "__main__":
    sys.exit(main())
