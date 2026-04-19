# PROJECT_inv_dataset_regpar_readme.md

## Script: `PROJECT_inv_dataset_regpar.py`

**Purpose:** AEMpyX regularisation parameter sweep.
Runs a set of fixed-τ inversions in parallel (one job per τ value)
to explore the effect of regularisation strength. Results can be
compared to choose an optimal parameter or to construct an L-curve.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `inspect`, `copy`, `multiprocessing` |
| Third-party | `numpy`, `scipy`, `joblib` (parallel) |
| AEMpyX modules | `version`, `aesys`, `util`, `inverse` |

---

## Workflow

```
Set system, errors, active channels
Set file list (util.get_data_list or explicit)
Set TestFix = [τ₁, τ₂, …] regularisation values to test
Build base ctrl_dict (TikhOpt)
Loop over data files
  ├─ Parallel: for each τ, deep-copy ctrl_dict, set Tau1=τ
  │            joblib.Parallel → inverse.run_inv_flightline per τ
  └─ Serial:   loop over τ, run inverse.run_inv_flightline
Output: one _results.npz per (file, τ) combination
```

---

## Key Parameters

| Variable | Description |
|---|---|
| `TestFix` | Array of τ values to sweep (e.g. `[1, 3, 10, 30, 100, 300, 1000]`) |
| `RegFun` | `'fix'` — fixed regularisation (τ is set from `TestFix`) |
| `Parallel` | `True` to parallelise across τ values with joblib |
| `Njobs` | Number of parallel workers |
| `Nlyr` | Number of model layers |
| `outstrng` | Encodes all key inversion parameters; `_t<τ>` appended per run |

---

## Output Files

```
<OutResDir>/<file><outstrng>_t<τ>_results.npz   (one per τ)
<OutResDir>/<file><outstrng>_t<τ>_ctrl.npz
```

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `datetime`, `process_time`, `time`, `getpass` unused | Removed |
| Serial branch | `ctrl_tmp['inversion'][3][0]` — `ctrl_tmp` only defined in Parallel branch | Changed to `str(reg[0])` using loop variable |
