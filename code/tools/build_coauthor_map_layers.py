"""Extract review-map layers from the original satellite and MPA source files.

These layers are context for checking site assignments, not analysis inputs.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import openpyxl
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/dataset_inventory_20261101/map_layers.json"
SOURCE = ROOT / "data/external/tom_bell_2026-10-20"
NETCDF = ROOT / "local_archive/metadata_audit_2026-10-16/CAkelpCanopyEnv_2021.nc"
areas = sorted(set(pd.read_csv(ROOT / "outputs/dataset_inventory_20261101/comparison_decisions.csv")
                   .query("Source == 'KFM'")["MPA"]))
areas.append("Carrington Point SMR")
stations = pd.read_csv(ROOT / "outputs/dataset_inventory_20261101/site_response_coverage.csv")
stations = stations[(stations.program == "KFM") & stations.station_latitude.notna()
                    & stations.station_longitude.notna()]
centers = stations.groupby("site_id")[['station_latitude', 'station_longitude']].first()
comparisons = pd.read_csv(ROOT / "outputs/dataset_inventory_20261101/comparison_decisions.csv")

mpas = gpd.read_file(SOURCE / "02_MPA_CA.shp")
mpas = mpas[mpas.NAME.isin([a.replace(" 2003", "") for a in areas])]
boundaries = {}
for _, row in mpas.iterrows():
    geom = row.geometry.simplify(0.00004, preserve_topology=True)
    coords = [geom.exterior.coords] if geom.geom_type == "Polygon" else [p.exterior.coords for p in geom.geoms]
    boundaries[row.NAME] = [[[round(x, 6), round(y, 6)] for x, y in ring] for ring in coords]

site_list = pd.read_csv(SOURCE / "03_Site_List.csv")
mask_names = dict(zip(site_list['Site_#'], site_list.Site.str.strip("'")))
book = openpyxl.load_workbook(SOURCE / "03_MPA_Ref_Locations.xlsx", read_only=True, data_only=True)
sheet = book["Pixel Locations"]
header = list(next(sheet.iter_rows(min_row=1, max_row=1, values_only=True)))
roles = list(next(sheet.iter_rows(min_row=2, max_row=2, values_only=True)))
mask_columns = {}
for col in range(1, len(header) - 1):
    if isinstance(header[col], (int, float)) and roles[col] in {"MPA", "Ref_1", "Ref_2"}:
        mask_columns[col] = (int(header[col]), roles[col])
masks = {name: {"MPA": [], "Ref_1": [], "Ref_2": []} for name in mask_names.values()}
for row in sheet.iter_rows(min_row=3, values_only=True):
    for col, (site_no, role) in mask_columns.items():
        lat, lon = row[col], row[col + 1]
        if isinstance(lat, (float, int)) and isinstance(lon, (float, int)):
            masks[mask_names[site_no]][role].append([round(lon, 5), round(lat, 5)])
masks = {name: value for name, value in masks.items() if name in boundaries}
assert set(boundaries) == {area.replace(" 2003", "") for area in areas}

with xr.open_dataset(NETCDF) as ds:
    lat = ds.lat.values
    lon = ds.lon.values
    periods = np.flatnonzero(ds.year.values == 2021)
    # Canopy presence is the most direct map layer: any positive canopy area in 2021.
    canopy = np.zeros(len(lat), dtype=bool)
    for t in periods:
        canopy |= ds.area.isel(time=int(t)).values > 0

# The delivered mask coordinates should fall on the same Landsat grid as the NetCDF.
mask_points = np.array([point for roles in masks.values() for group in roles.values() for point in group])
grid = cKDTree(np.column_stack((lon * 93, lat * 111.32)))
distance_km, _ = grid.query(np.column_stack((mask_points[:, 0] * 93, mask_points[:, 1] * 111.32)))
if distance_km.max() > 0.002:
    raise ValueError(f"Satellite mask is misaligned with the kelp grid: max {distance_km.max() * 1000:.1f} m")

kelp = {}
for area in areas:
    group = comparisons[(comparisons.Source == "KFM") & (comparisons.MPA == area)]
    codes = set()
    for field in ("mapped_protected_sites", "mapped_reference_sites"):
        for value in group[field].dropna():
            codes.update(str(value).split(";"))
    if area == "Carrington Point SMR":
        codes.add("a-k-05")
    located = centers.loc[centers.index.intersection(codes)]
    if located.empty:
        kelp[area] = []
        continue
    near = np.zeros(len(lat), dtype=bool)
    for _, site in located.iterrows():
        dlat = (lat - site.station_latitude) * 111.32
        dlon = (lon - site.station_longitude) * 111.32 * np.cos(np.deg2rad(site.station_latitude))
        near |= dlat * dlat + dlon * dlon <= 12 * 12
    ids = np.flatnonzero(near & canopy)
    kelp[area] = [[round(float(lon[i]), 5), round(float(lat[i]), 5)] for i in ids]

payload = {
    "provenance": {"canopy": str(NETCDF.relative_to(ROOT)), "period": "2021 Q1–Q4",
                   "definition": "30 m Landsat grid pixel with canopy area > 0 in at least one 2021 quarter",
                   "boundaries": str((SOURCE / "02_MPA_CA.shp").relative_to(ROOT)),
                   "boundary_metadata_date": "2018-12-11",
                   "masks": str((SOURCE / "03_MPA_Ref_Locations.xlsx").relative_to(ROOT)),
                   "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in (NETCDF, SOURCE / "02_MPA_CA.shp",
                                           SOURCE / "03_MPA_Ref_Locations.xlsx")}},
    "boundaries": boundaries, "canopy": kelp, "satellite_masks": masks,
}
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
print("Wrote", OUT, "bytes", OUT.stat().st_size,
      "canopy pixels", {k: len(v) for k, v in kelp.items()},
      "mask areas", sorted(masks))
