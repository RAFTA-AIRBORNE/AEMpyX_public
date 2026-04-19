# viz_readme.md

## Module: `viz.py`

**Purpose:** Visualisation utilities for AEMpyX results.
Produces publication-quality plots of data, model ensembles,
depth profiles, flight-line sections, and geographic exports
(GeoTIFF, KMZ).

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

Last change: vr Apr 2026

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `sys`, `math` |
| Third-party | `numpy`, `numpy.ma`, `matplotlib` (pyplot, colors, cm, axes_grid1, ticker), `simplekml`, `aesys`, `inverse`, `util` |

---

## Naming Conventions

- Input parameters use CapitalizedCamelCase (`Data`, `Model`, `Fontsizes`).
- Figure and axis objects follow matplotlib conventions.
- Boolean output flags: `OutInfo`, `out`, `save`.

---

## Function Reference

### Ensemble plots

**`plot_model_ensemble(Data, Model, Ctrl, ...)`**
Plot a model ensemble (e.g. from MCMC or uncertainty propagation)
as a shaded percentile band with median and mean overlaid.
Supports log or linear resistivity axis.

**`plot_data_ensemble(Data, Ctrl, ...)`**
Plot data ensemble (observed, calculated, and spread) per component
across a flight line or site range.

---

### Single-realisation plots

**`plot_data(Data, Ctrl, ...)`**
Plot observed vs. calculated data for a single inversion result,
with normalised RMS annotation.

**`plot_model(Model, Ctrl, ...)`**
Plot a 1-D layered resistivity model with optional prior and
sensitivity shading.

**`plot_depth_prof(Data, Model, Ctrl, ...)`**
Plot resistivity depth profiles along a flight line as a
pseudo-section (colour-filled column plot).

**`plot_matrix(Matrix, Ctrl, ...)`**
Plot a matrix (e.g. covariance, Jacobian, resolution) as a
colour-mapped image.

---

### System-specific flight-line plots

**`plot_data_aem05(Data, Ctrl, ...)`**
Plot per-component data for the AEM05 FD system over a flight line.

**`plot_data_genesis(Data, Ctrl, ...)`**
Plot per-component data for the Genesis TD system over a flight line.

**`plot_flightline_aem05(Data, Model, Ctrl, ...)`**
Combined flight-line section plot for AEM05: data panels + resistivity
pseudo-section.

**`plot_flightline_genesis(Data, Model, Ctrl, ...)`**
Combined flight-line section plot for Genesis: data panels + resistivity
pseudo-section.

---

### Geographic export

**`save_geotiff(filename, data, extent, crs, ...)`**
Save a raster array as a GeoTIFF with defined spatial extent and CRS.

**`make_kml(kmzfile, ...)`**
Create a KMZ file (Google Earth overlay) from a raster or flight-line
data set using `simplekml`.

**`gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels)`**
Create a matplotlib figure sized to match a geographic bounding box
for Google Earth export.
Returns `(fig, ax)`.

---

## Usage Notes

- All plot functions accept a `save` or `SaveFig` parameter; when set,
  figures are written to file rather than displayed interactively.
- Colour maps, font sizes, and axis limits are controlled through a
  `Ctrl` dict whose keys vary by function — see individual docstrings.
- `\Omega` in axis labels must use raw strings (`r'...'`); corrected
  in two functions (Apr 2026).

---

## Known Issues / Notes

- Removed unused imports: `os`, `inspect`, `process_time`, `datetime`,
  `warnings`, `cycler`; removed duplicate `import numpy`. (Apr 2026)
- `\Omega` escape sequence in `set_xlabel` calls fixed to raw strings.
  (Apr 2026)

---

## Change Log

| Date | Author | Note |
|---|---|---|
| 2020–2021 | vrath / duygu | Initial version (AEMpyX). |
| 2022–2023 | vrath / duygu | Genesis support; ensemble plots; geo-export. |
| Apr 2026 | vrath | Removed unused imports; escape fixes; provenance; readme. |
