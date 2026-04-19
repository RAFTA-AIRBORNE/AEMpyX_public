# PROJECT_viz_models_area_readme.md

## Script: `PROJECT_viz_models_area.py`

**Purpose:** AEMpyX spatial area model visualisation.
Loads merged inversion results and produces colour-mapped depth-slice
or layer-property maps over a geographic area.

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
| Third-party | `numpy`, `matplotlib`, `scipy` (interpolate, spatial), `matplotlib.axis`, `mpl_toolkits.axes_grid1` |
| AEMpyX modules | `version`, `util`, `aesys`, `inverse` |

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `skgstat`, `shapely`, `viz` unused | Removed |
