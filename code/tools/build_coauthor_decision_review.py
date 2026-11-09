"""Build KFM site-map data for the current analysis decisions page."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "outputs/dataset_inventory_20261101"


def records(name: str) -> list[dict]:
    frame = pd.read_csv(DATA / f"{name}.csv", low_memory=False).fillna("")
    return frame.to_dict("records")


comparisons = records("comparison_decisions")
year_roster = records("integrated_comparison_year_roster")
sites = records("site_response_coverage")
land_source = ROOT / "data/external/map_context/ne_10m_land_channel_islands.geojson"
land_features = json.loads(land_source.read_text(encoding="utf-8"))["features"]
map_layers_path = DATA / "map_layers.json"
if not map_layers_path.exists():
    raise FileNotFoundError(f"Build satellite and boundary layers first: {map_layers_path}")
map_layers = json.loads(map_layers_path.read_text(encoding="utf-8"))
land_rings = []
for feature in land_features:
    geometry = feature["geometry"]
    polygons = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"]
    land_rings.extend(polygon[0] for polygon in polygons)
assert len(land_rings) >= 5
issues = [r for r in records("unresolved_issues") if r["issue_id"] not in {"U12", "U16", "U17"}]
rules = records("current_decisions_and_reasons")
names = {r["site_id"]: r["site_name"] for r in sites if r["program"] == "KFM"}
coordinates = {}
for station in sites:
    if station["program"] != "KFM" or station["station_latitude"] == "" or station["station_longitude"] == "":
        continue
    point = {"lat": float(station["station_latitude"]), "lon": float(station["station_longitude"])}
    previous = coordinates.setdefault(station["site_id"], point)
    if abs(previous["lat"] - point["lat"]) > 1e-6 or abs(previous["lon"] - point["lon"]) > 1e-6:
        raise ValueError(f"Conflicting KFM coordinates for {station['site_id']}")
assert len(names) >= 30 and names["a-k-01"] == "Wyckoff Ledge"
roster_by_series: dict[str, list[dict]] = {}
for record in year_roster:
    if record["Source"] == "KFM":
        roster_by_series.setdefault(record["series_id"], []).append(record)

for row in comparisons:
    if row["Source"] != "KFM":
        continue
    for role in ("protected", "reference"):
        codes = [code for code in str(row[f"mapped_{role}_sites"]).split(";") if code]
        missing = [code for code in codes if code not in names]
        if missing:
            raise ValueError(f"Unmapped KFM {role} code(s) in {row['series_id']}: {missing}")
        row[f"{role}_named"] = [{"code": code, "name": names[code]} for code in codes]
    supported_years = {int(y) for y in str(row["integrated_years"]).split(";") if y}
    used = [r for r in roster_by_series.get(row["series_id"], [])
            if int(r["year"]) in supported_years]
    for role in ("protected", "reference"):
        field = f"{role}_stations"
        years_by_code: dict[str, set[int]] = {}
        for survey in used:
            for code in str(survey[field]).split(";"):
                if code:
                    years_by_code.setdefault(code, set()).add(int(survey["year"]))
        actual_codes = sorted(years_by_code)
        missing = [code for code in actual_codes if code not in names]
        if missing:
            raise ValueError(f"Unmapped contributing KFM {role} code(s): {missing}")
        row[f"{role}_used"] = [{"code": code, "name": names[code], "years": sorted(years_by_code[code])} for code in actual_codes]

kfm = [r for r in comparisons if r["Source"] == "KFM"]
groups: dict[str, list[dict]] = {}
for row in kfm:
    groups.setdefault(row["MPA"], []).append(row)
groups.setdefault("Carrington Point SMR", [])

site_notes = {
    "Gull Island SMR": {
        "assumption": "Gull Island and Yellow Banks form the mapped KFM pair. Yellow Banks is about 15 miles away. The KFM histories are excluded; the PISCO Gull comparison is separate.",
        "check": "Does this exclusion reflect the sites and habitat differences correctly?",
    },
    "Harris Point SMR": {
        "assumption": "Hare Rock is compared with Wyckoff Ledge for the retained benthic and kelp responses; sheephead also uses Miracle Mile. The protected site is north of the reference sites.",
        "check": "Are these locations reasonable matches for each response?",
    },
    "Santa Barbara Island SMR": {
        "assumption": "Webster’s Arch is an offshore pinnacle. It is already excluded from benthic and kelp responses, but contributes to the fish comparison in the current review.",
        "check": "Should Webster’s Arch also be left out of fish comparisons?",
    },
    "Carrington Point SMR": {
        "assumption": "Rodes Reef is mapped as a possible reference but contributes no measurements to the current fitted histories.",
        "check": "Is this treatment of Rodes Reef correct?",
    },
}

payload = {"groups": groups, "comparisons": comparisons, "names": names, "coordinates": coordinates,
           "land_rings": land_rings, "map_layers": map_layers, "issues": issues, "rules": rules, "site_notes": site_notes,
           "counts": {"candidate": len(comparisons), "kfm": len(kfm)}}
safe_json = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))

(DATA / "site_map_payload.json").write_text(safe_json, encoding="utf-8")
print(f"Wrote {DATA / 'site_map_payload.json'}: {len(kfm)} KFM histories, {len(names)} station codes")
