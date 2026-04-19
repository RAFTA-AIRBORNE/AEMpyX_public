# PROJECT_viz_data_area_readme.md

## Script: `PROJECT_viz_data_area.py`

**Purpose:** AEMpyX spatial area data visualisation and interpolation.
Reads flight-line data, interpolates selected components onto a regular
grid (using KDTree nearest-neighbour, griddata, or RBF), and produces
colour-mapped area plots.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `inspect`, `getpass` |
| Third-party | `numpy`, `matplotlib`, `scipy` (interpolate, spatial) |
| AEMpyX modules | `version`, `util`, `aesys` |

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `warnings`, `shapely`, `viz` unused | Removed |
