"""Compare public release 9 with local inputs without replacing analysis data."""
import hashlib
import io
import json
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd

ANALYSIS = Path(__file__).resolve().parents[2]
HERE = ANALYSIS / "data/external/pisco_release9"
HERE.mkdir(parents=True, exist_ok=True)
url = "https://opc.dataone.org/metacat/d1/mn/v2/object/doi:10.25494%2FP6%2FMLPA_kelpforest.9"
metadata = urllib.request.urlopen(url, timeout=60).read()
root = ET.fromstring(metadata)
assert root.get("packageId") == "doi:10.25494/P6/MLPA_kelpforest.9"
(HERE / "MLPA_kelpforest_release9.eml.xml").write_bytes(metadata)
checks = []
for table in root.findall("./dataset/dataTable"):
    name = table.findtext("./physical/objectName")
    if name not in ("MLPA_kelpforest_fish.6.csv", "MLPA_kelpforest_swath.7.csv"):
        continue
    local = ANALYSIS / "data/PISCO" / ("MLPA_kelpforest_fish_2024.csv" if "fish" in name else "MLPA_kelpforest_swath_2024.csv")
    remote_url = table.findtext("./physical/distribution/online/url")
    cache = Path("/tmp") / name
    public = cache.read_bytes() if cache.exists() else urllib.request.urlopen(remote_url, timeout=120).read()
    if not cache.exists():
        cache.write_bytes(public)
    expected = table.findtext("./physical/authentication")
    assert hashlib.md5(public).hexdigest() == expected, name
    private = local.read_bytes()
    def csv(z):
        try:
            return pd.read_csv(io.BytesIO(z), dtype=str, keep_default_na=False)
        except UnicodeDecodeError:
            return pd.read_csv(io.BytesIO(z), dtype=str, keep_default_na=False, encoding="cp1252")
    a, b = csv(private), csv(public)
    common = sorted(set(a.columns) & set(b.columns))
    record_match = False
    if set(a.columns) == set(b.columns) and len(a) == len(b):
        ha = pd.util.hash_pandas_object(a[common], index=False).sort_values().to_numpy()
        hb = pd.util.hash_pandas_object(b[common], index=False).sort_values().to_numpy()
        record_match = bool((ha == hb).all())
    checks.append({"local_file": str(local.relative_to(ANALYSIS)), "public_entity": name,
        "public_md5": expected, "local_md5": hashlib.md5(private).hexdigest(),
        "bytes_match": private == public, "csv_records_match": record_match,
        "local_rows": len(a), "public_rows": len(b), "local_columns": list(a.columns),
        "public_columns": list(b.columns), "public_entity_url": remote_url,
        "local_only_columns": sorted(set(a.columns) - set(b.columns)),
        "ordered_cell_differences": {c: int((a[c] != b[c]).sum()) for c in common} if len(a) == len(b) else None})
    print(name, "bytes match", private == public, "CSV records match", record_match, "rows", len(a), len(b))
(HERE / "pisco_release_check.json").write_text(json.dumps({"doi": "10.25494/P6/MLPA_kelpforest.9",
    "published": root.findtext("./dataset/pubDate"), "checks": checks}, indent=2) + "\n")
