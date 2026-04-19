# PROJECT_pre_import_data_readme.md

## Script: `PROJECT_pre_import_data.py`

**Purpose:** AEMpyX data ingestion and flight-line extraction.
Reads raw XYZ survey files, optionally reprojects coordinates,
selects a spatial subset (rectangle, polygon, union, intersection,
or explicit line list), and writes individual flight-line NPZ files.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `inspect` |
| Third-party | `numpy` |
| AEMpyX modules | `version`, `util`, `aesys` |

---

## Workflow

```
Set AEM_system, survey, directories
Set DataSelect mode
Loop over XYZ files
  └─ aesys.read_survey_data → Data
  └─ optional ITM→UTM reprojection  (SetProj=True)
Apply spatial selection
  ├─ 'Rectangle'    → util.extract_data_rect
  ├─ 'Polygon'      → util.extract_data_poly
  ├─ 'Union/Intersection' → util.modify_polygon + extract_data_poly
  └─ 'Lines'        → numpy.where on line column
Optional full-data output  (MergeOut=True)
Loop over flight lines
  └─ direction correction  (CorrectDirection=True)
  └─ NaN check  (CheckNaN=True)
  └─ aesys.write_aempy → _FL<n>_data.npz
```

---

## Key Parameters

| Variable | Description |
|---|---|
| `DataSelect` | `'Rectangle'`, `'Polygon'`, `'Union'`, `'Intersection'`, `'Lines'` |
| `LineList` | Explicit flight-line numbers for `'Lines'` mode |
| `FlightlineAngle` | Expected flight direction (degrees from North) |
| `Spread` | Tolerance around `FlightlineAngle` before reversal |
| `CorrectDirection` | Flip flight line if direction deviates |
| `CheckNaN` | Skip lines with any NaN values |
| `LinesMin` | Minimum number of sites per output flight line |
| `MergeOut` | Also write a full merged data file |

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `datetime`, `warnings`, `getpass` unused | Removed |
