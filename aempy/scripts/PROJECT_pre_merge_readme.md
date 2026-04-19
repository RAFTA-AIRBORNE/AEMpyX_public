# PROJECT_pre_merge_readme.md

## Script: `PROJECT_pre_merge.py`

**Purpose:** AEMpyX data and model set merging.
Concatenates multiple flight-line NPZ files (data or inversion results)
into a single merged NPZ, with optional quality thresholding for model sets.

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
| AEMpyX modules | `version`, `util`, `inverse` |

---

## Workflow

```
Set DataType ('data' or 'models')
Set file list (util.get_data_list)
├─ DataType = 'data':
│   inverse.merge_data_sets → merged NPZ
└─ DataType = 'models':
    inverse.merge_model_sets(qthresh=Thresh) → merged NPZ
```

`Thresh = ['smp'|'rms', value]` filters sites by spread-map value or
normalised RMS before merging.

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Path fix | `InDatDir = OutDatDir + '/'` (wrong variable) | `OutDatDir = OutDatDir + '/'` |
| `merge_data_sets` | `file_list=` → `infile_list=` (correct kwarg name) | Fixed |
| Imports | `process_time`, `datetime`, `warnings`, `aesys` unused | Removed |
