# PROJECT_interpolate_data_area_readme.md

## Script: `PROJECT_interpolate_data_area.py`

**Purpose:** AEMpyX spatial data interpolation and area mapping.
Reads AEM flight-line data, interpolates selected components onto a
regular mesh or profile using KD-tree nearest-neighbour lookup, and
reports data statistics. Intended as a precursor to area plotting.

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
| Third-party | `numpy`, `matplotlib` (pyplot, ticker), `mpl_toolkits.axes_grid1`, `scipy` (interpolate, spatial) |
| AEMpyX modules | `version`, `util`, `aesys` |

---

## Workflow

```
Set system  (AEM_system block)
Set data source  (InDatDir, DataFile, DataSet)
  ├─ DataSet = 'mesh'  → define grid dimensions (numIndexes)
  └─ DataSet = 'prof'  → read profile coordinates from CSV (NumIndexes)
Loop over data files
  └─ aesys.read_aempy → Data
  └─ extract E, N, Z coordinates
  └─ if 'mesh':
       build KDTree on data points
       query mesh points → nearest-neighbour index
       DI = Dats[index].reshape(XI.shape)
  └─ loop over components (CompList)
       extract D = Data[:, comp]
       report min/max of raw and interpolated data
```

---

## Key Parameters

| Variable | Description |
|---|---|
| `AEM_system` | `'aem05'` or `'genesis'` |
| `DataSet` | `'mesh'` (2-D grid) or `'prof'` (profile) |
| `numIndexes` | `[nX, nY]` grid dimensions for `'mesh'` mode |
| `NumIndexes` | `[nPts]` number of points for `'prof'` mode |
| `XYFact` | Coordinate scaling factor (default `1e-3` → km) |
| `CompList` | List of component names to process (e.g. `['P1','Q1',…,'ALT']`) |
| `step` | Decimation step for data reading |
| `InDatDir` | Input data directory |
| `DataFile` | Filename or glob pattern |

---

## DataSet Modes

### `'mesh'` mode
Builds a regular `(nX × nY)` grid over the data extent and uses
`scipy.spatial.KDTree` to find the nearest data point for each grid
cell. The interpolated array `DI` has shape `(nX, nY)`.

### `'prof'` mode
Reads a profile coordinate CSV (one `x,y` pair per line). Currently
interpolates onto a 1-D linspace between data extent limits using
`NumIndexes[0]` points. Profile-specific interpolation (along-profile
distance) can be added.

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| `DI` variable | Used in `nanmin`/`nanmax` but never assigned | Added `DI = Dats[index].reshape(XI.shape)` after KDTree query |
| `NumIndexes` vs `numIndexes` | Profile branch defined `NumIndexes` but linspace used `numIndexes` | Unified to `NumIndexes` |
| Imports | `datetime`, `warnings`, `getpass`, `matplotlib.axis`, `mpl_toolkits.axes_grid1`, `shapely`, `viz`, `inverse` unused | Removed |

---

## Usage Notes

- The `PLM` and `ALT` components (last two in `CompList`) are excluded
  from the data loop via `range(len(CompList)-2)`.
- `DataFile` may be a glob pattern; the script currently wraps it in a
  single-element list — full `util.get_data_list` integration can be
  substituted for wildcard expansion.
- Compatible with jupytext (`.py` ↔ `.ipynb`).
