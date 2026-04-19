# PROJECT_viz_kmz_readme.md

## Script: `PROJECT_viz_kmz.py`

**Purpose:** AEMpyX KMZ/KML flight-line export for Google Earth.
Reads flight-line NPZ data, optionally links associated PNG plots, and
writes a compressed KMZ file with per-site markers and flight-line
folder structure.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `csv`, `inspect`, `getpass` |
| Third-party | `numpy`, `simplekml`, `datetime` |
| AEMpyX modules | `version`, `aesys`, `util` |

---

## Workflow

```
Set directories, SearchStrng
util.get_filelist → data_files
Open simplekml.Kml
Loop over flight-line files
  └─ aesys.read_aempy → Data
  └─ util.project_utm_to_latlon → lat, lon
  └─ add site markers every MarkEvery steps
  └─ link PNG plot (AddImages=True)
  └─ mark start/end/centre points
Optional: add special locations (AddSpecial=True)
kml.savekmz → KMZFile.kmz
```

---

## Key Parameters

| Variable | Description |
|---|---|
| `MarkEvery` | Mark every N-th site with a data icon |
| `AddImages` | Link flight-line PNG plots as balloon images |
| `AddSpecial` | Add custom points from `Special.dat` CSV |
| `MarkStartPoints` / `MarkEndPoints` / `MarkCenterPoints` | Toggle line endpoint markers |
| `ImageWidth` | Width of embedded PNG in pixels |

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| `imstring` | HTML string used mismatched single-quote delimiters causing `SyntaxError` | Rewrote with escaped quotes |
| `imstring` | Used in `MarkStartPoints`/`EndPoints`/`CenterPoints` blocks before assignment when `AddImages=False` | Initialised to `''` before loop |
