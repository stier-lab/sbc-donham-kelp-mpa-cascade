# Land context for the coauthor site map

The coauthor review page embeds `ne_10m_land_channel_islands.geojson` as generalized land context. It is **not** an MPA boundary, habitat map, or source of station coordinates. Station positions come from the KFM site crosswalk in the inventory tables.

Source: Natural Earth 1:10m land polygons, version 5.1.1, downloaded November 4, 2026 from `https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_land.zip`. SHA-256 of the downloaded ZIP: `e547d749445eaa0964aba76738090ec88f5e63c4585122170f98c67a7ea922dc`. Natural Earth data are public domain.

Processing: intersected the land polygons with WGS84 bounds 120.6°–119.0° W, 33.35°–34.25° N; transformed to UTM zone 11N (EPSG:32611); simplified at 30 m with topology preserved; transformed back to WGS84; wrote GeoJSON. The map draws these polygons in an aspect-correct local kilometer projection. The simplification and source scale are appropriate for orientation, not for measuring shoreline distance or site habitat.
