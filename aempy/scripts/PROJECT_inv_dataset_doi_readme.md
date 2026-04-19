# PROJECT_inv_dataset_doi_readme.md

## Script: `PROJECT_inv_dataset_doi.py`

**Purpose:** Depth-of-investigation (DoI) analysis for AEMpyX inversion results.
Loads two or more sets of inversion results (produced with different
prior resistivities) and computes a per-site DoI estimate using one
of two methods. Results are stored as a single NPZ file.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `inspect`, `time`, `warnings` |
| Third-party | `numpy` |
| AEMpyX modules | `version`, `util`, `inverse` |

---

## DoI Methods

### `'Oldenburg1999'`

Based on the ratio of model difference to reference difference
(Oldenburg & Li, 1999):

```
v_doi = |Δlog₁₀ m| / |Δlog₁₀ m_ref|
```

where `Δ` is taken between the two inversion results. High values
indicate poor constraint (model is just reflecting the prior).

### `'Variance'`

Computes the ensemble variance across model runs:

```
m_ano = log₁₀(m_i) - mean(log₁₀(m))
v_doi = sqrt(var(m_ano))
```

High variance indicates sensitivity to the choice of prior, i.e.,
shallow DoI at that depth.

---

## Workflow

```
Set Method ('Oldenburg1999' or 'Variance')
Set file list  (util.get_data_list or explicit)
Read first result file → initialise output NPZ (Fileout)
Loop over result files
  └─ numpy.load → site_modl, site_merr, mod_ref, mod_act
  └─ inverse.get_nlyr, inverse.extract_mod
  └─ stack modls, merrs, mrefs
Compute DoI
  ├─ 'Oldenburg1999': log10 diff ratio
  └─ 'Variance':      log10 anomaly variance
Write DoI arrays to Fileout via util.add_object_npz
```

---

## Key Parameters

| Variable | Description |
|---|---|
| `Method` | `'Oldenburg1999'` or `'Variance'` |
| `FileList` | `'search'` (glob) or `'set'` (explicit list) |
| `SearchStrng` | Glob pattern for result files (e.g. `'*results.npz'`) |
| `InModDir` | Directory containing inversion result files |
| `Fileout` | Output NPZ path |

---

## Output File (`Fileout`)

Initial keys (written from first result file):

```
fline, m_act, m_ref
site_x, site_y, site_z, site_dact, site_dobs, site_derr
doimethod
```

Keys added per method:

| Key | Oldenburg1999 | Variance |
|---|:---:|:---:|
| `doi_files` | ✓ | ✓ |
| `site_doi_mavg` | ✓ | ✓ |
| `site_doi_models` | ✓ | ✓ |
| `site_doi_mrefs` | ✓ | ✓ |
| `site_doi_merrs` | ✓ | ✓ |
| `site_doi_doi` | ✓ | ✓ |
| `site_doi_mdif` | ✓ | — |
| `site_doi_rdif` | ✓ | — |
| `site_doi_ano` | — | ✓ |
| `site_doi_var` | — | ✓ |
| `doi_meth` | ✓ | ✓ |

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `scipy` imported but unused | Removed |
| `Method` | Assigned `'Oldenburg1999'` then immediately overridden by `'Variance'` — dead assignment | Collapsed to single assignment with comment |
| Duplicate jupytext header | Two `# ---` jupytext blocks (incl. stray `# ! /usr/bin/python`) | Removed duplicate |

---

## References

Oldenburg, D.W. & Li, Y. (1999).
Estimating depth of investigation in DC resistivity and IP surveys.
*Geophysics*, 64(2), 403–416.

---

## Usage Notes

- Requires **at least two** result files for `'Oldenburg1999'` (difference
  of two priors). `'Variance'` generalises to any number ≥ 2.
- All result files must have the same model dimensions (`modshp`);
  the script exits if they differ.
- Results are appended to `Fileout` incrementally via `util.add_object_npz`,
  so a partial run will leave a valid but incomplete NPZ.
- Compatible with jupytext (`.py` ↔ `.ipynb`).
