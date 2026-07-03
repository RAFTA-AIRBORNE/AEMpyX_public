# PROJECT_viz_kml_overlay_readme.md

## Script: `PROJECT_viz_kml_overlay.py`

**Purpose:** AEMpyX KML/KMZ ground-overlay generation. Reads flight-line
data, interpolates selected components onto a regular grid (using
`griddata` or RBF), renders each component as a borderless, transparent,
georeferenced PNG matching its data bounding box, reprojects the box to
WGS84 lon/lat, and writes the corresponding KML `GroundOverlay` element(s).
Individual `.kml`/`.png` pairs are written per component; optionally all
components are bundled into a single `.kmz` (with a `doc.kml` containing
one `GroundOverlay` per component) for direct loading in Google Earth.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `shutil`, `zipfile`, `getpass`, `inspect` |
| Third-party | `numpy`, `matplotlib`, `scipy` (interpolate, spatial), `pyproj` |
| AEMpyX modules | `version`, `util`, `aesys`, `inverse` |

---

## Key settings

| Variable | Purpose |
|---|---|
| `EPSGCode` | Source projected CRS of Easting/Northing in the data files (must match the survey's processing CRS) |
| `WGS84Code` | Target CRS for the KML `LatLonBox`, always `'epsg:4326'` |
| `CompList` | Components to render, with optional `[vmin, vmax, step]` colour limits |
| `InterpMethod` | `griddata`/`rbf` interpolation onto the regular mesh |
| `MaskNeg`, `MaskPoly`, `MaskDist` | Masking options carried over from `PROJECT_viz_data_area.py` |
| `KmzBundle` | If `True`, bundle all per-component overlays into one `.kmz` |

---

## Notes

- The overlay PNG canvas is rendered with `ax.add_axes([0,0,1,1])`, no
  ticks/frame/margin, and `transparent=True`, so it aligns exactly with
  the `LatLonBox` declared in the KML — any deviation here (e.g. default
  matplotlib margins) will misplace the raster in Google Earth.
- Coordinates are transformed from the survey's projected CRS to WGS84
  lon/lat via `pyproj.Transformer`; `EPSGCode` must be set correctly for
  the dataset in use (e.g. Irish Transverse Mercator `epsg:2157`, UTM
  zone 29N `epsg:32629`).
- Follows the same file-discovery, data-merging, and interpolation/mesh
  logic as `PROJECT_viz_data_area.py` for consistency across the
  visualisation scripts.

---

## Bugs Fixed / Design Notes (Jul 2026)

| Location | Issue | Fix |
|---|---|---|
| Overlay rendering | Default matplotlib axes/margins shift raster relative to declared bbox | Borderless full-canvas axes (`add_axes([0,0,1,1])`), `axis('off')`, `pad_inches=0` |
| Georeferencing | `LatLonBox` requires WGS84 lon/lat, data are in a projected CRS | `pyproj.Transformer` reprojection of the interpolated grid's bounding box |
