# PROJECT_pre_process_flightline_readme.md

## Script: `PROJECT_pre_process_flightline.py`

**Purpose:** AEMpyX flight-line data pre-processing pipeline.
Applies a sequence of filters, flagging, NaN handling, and SVD-based
noise filtering to raw AEM data, writing intermediate and final outputs.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `copy`, `inspect` |
| Third-party | `numpy` |
| AEMpyX modules | `version`, `util`, `prep`, `aesys` |

---

## Processing Pipeline (per flight line)

```
aesys.read_aempy
├─ prep.filter_column    altitude + DEM lowpass (Butterworth n=4)
├─ prep.filter_column    PLM lowpass
├─ prep.insert_flag      PLM threshold (plm)
├─ prep.insert_flag      data < 0 threshold (less)
├─ prep.insert_flag      altitude > max threshold (great)
├─ [optional: write NaN file  OutNaN=True]
├─ prep.handle_gaps      delete / impute NaN rows
└─ loop k=1..SingValMax:
     prep.calc_svd_decomp   (k-PC reconstruction)
     aesys.write_aempy      _k<n>_data.npz
     [optional: write residual  OutRes=True]
```

---

## Key Parameters

| Variable | Description |
|---|---|
| `SingValMax` | Maximum number of SVD components to test |
| `OutNaN` | Write data-with-NaN file before gap handling |
| `OutRes` | Write SVD residual file for each k |
| `plmthresh` | PLM flagging threshold |

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| `prep.calc_svd_decomp` | `k=` → `K=`, `out_full=` → `OutFull=` (capitalised params) | Fixed |
| `prep.insert_flag` | Positional args changed to keyword args matching updated API | Fixed |
| `prep.filter_column` | `method=` → `Method=`, `columns` → `Columns=` | Fixed |
| `prep.handle_gaps` | `columns` positional → `Columns=` | Fixed |
| Imports | `getpass` unused | Removed |
