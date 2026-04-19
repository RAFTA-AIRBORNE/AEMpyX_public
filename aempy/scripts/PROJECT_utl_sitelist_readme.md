# PROJECT_utl_sitelist_readme.md

## Script: `PROJECT_utl_sitelist.py`

**Purpose:** MT site list utility for WALDIM and related analyses.
Reads a directory of EDI files, extracts station coordinates (lat/lon
or UTM), and writes a CSV site list with optional KML styling columns.

> **Note:** This script belongs to the **py4mt / mtpy** ecosystem,
> not AEMpyX. It uses `PY4MT_ROOT` / `PY4MT_DATA` environment
> variables and the `mtpy` library. It is included here as a companion
> utility in the same project tree.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `csv`, `inspect` |
| Third-party | `numpy` (as `np`), `mtpy` (`mtpy.core.mt.MT`), `utm` |
| py4mt modules | `util`, `version` |

---

## Workflow

```
Set EdiDir, CSVFile, to_utm, to_kml
os.listdir → .edi files
Loop over EDI files
  └─ MT().read(file) → lat, lon, elev
  └─ utm.from_latlon (if to_utm)
  └─ csv.writer.writerow → [E, N, name] + to_kml columns
```

---

## Output

`CSVFile` (e.g. `Sitelist.dat`) with columns:
`Easting, Northing, SiteName, [KML colour, scale, icon, ...]`

---

## Usage Notes

- Requires `PY4MT_ROOT` and `PY4MT_DATA` environment variables.
- `to_kml` list appended to each row controls Google Earth styling
  when the CSV is consumed by a KML generation script.
- No changes made to environment variable references or module paths —
  these are correct for the py4mt ecosystem.
