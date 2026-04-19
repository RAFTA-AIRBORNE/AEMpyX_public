# PROJECT_inv_dataset_lcp_readme.md

## Script: `PROJECT_inv_dataset_lcp.py`

**Purpose:** AEMpyX Laterally Correlated Procedure (LCP).
Applies spatial lateral correlation to a set of independently
inverted 1-D models, after Christensen & Tølbøll (2009).
Operates layer-wise using overlapping tiles and Gaussian
covariance weighting.

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
| Third-party | `numpy`, `scipy` (interpolate, spatial, linalg), `datetime` |
| AEMpyX modules | `version`, `util`, `aesys`, `inverse` |

---

## Workflow

```
Set system (AEM_system)
Load merged model file (inverse.merge_model_sets)
Loop over model files
  └─ numpy.load → e, n, d, m, c arrays
  └─ optional log10 transform (ParaTrans=1)
  Step 1: tile-based lateral correlation
    └─ overlapping tiles (TileSize × TileOverlap)
    └─ scipy.spatial.distance.pdist → inter-site distances
    └─ Gaussian covariance: cov_s = exp(dist/Scale)
    └─ layer-wise or full-matrix MAP update
    └─ mod_cor accumulated
  Step 3: (placeholder) forward recheck or re-inversion
  numpy.savez_compressed → corrfile with mod_cor added
```

---

## Key Parameters

| Variable | Description |
|---|---|
| `TileSize` | Tile side length (m) |
| `TileOverlap` | Fractional overlap between adjacent tiles (0–1) |
| `TileMinSites` | Minimum sites in a tile to apply LCP |
| `CovarThresh` | Correlation length (m); `Scale = 0.5 * CovarThresh` |
| `LayerWise` | If `True`, correlation applied per layer independently |
| `ParaTrans` | `1` = log10-transform models before correlation |
| `ReCalc` | `'fwd'` or `'inverse'` — post-correlation step (placeholder) |
| `MergeModels` | If `True`, merge all input files before processing |

---

## Output Files

```
<corrfile>   (modified in-place via numpy.savez_compressed)
    All original keys from input NPZ
    mod_cor  : laterally correlated model array (nsites × nlyr)
    header   : updated header string with date
```

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `shapely`, `viz`, `inverse` unused initially | Removed `shapely`, `viz`; `inverse` re-added for correct `merge_model_sets` call |
| `util.merge_model_sets` | Function lives in `inverse`, not `util` | Changed to `inverse.merge_model_sets` |
| `dictout=True` arg | Not supported by `inverse.merge_model_sets` | Removed |
| Step 3 timing | `elapsed` printed before re-assignment; duplicate `elapsed = process_time()` | Replaced with `start_step3` / `elapsed_step3` / `elapsed_total` |

---

## Reference

Christensen, N. B. & Tølbøll, R. J. (2009).
A lateral model parameter correlation procedure for one-dimensional
inverse modelling.
*Geophysical Prospecting*, 57, 919–929.
doi: 10.1111/j.1365-2478.2008.00756.x
