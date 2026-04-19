# PROJECT_inv_dataset_readme.md

## Script: `PROJECT_inv_dataset.py`

**Purpose:** AEMpyX dataset inversion for the AEM05 FD system.
Inverts a set of flight-line NPZ files using Tikhonov regularisation
(`TikhOpt`), with optional parallel execution via `joblib`.
Dispatches to `inverse.run_inv_flightline`.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `inspect`, `multiprocessing` |
| Third-party | `numpy`, `scipy`, `joblib` (optional, parallel) |
| AEMpyX modules | `version`, `aesys`, `util`, `inverse` |

---

## Workflow

```
Set system, data errors, active channels  (AEM_system block)
Set file list  (util.get_data_list or explicit list)
Set regularisation  (RegFun, Tau0, Tau1 ranges)
Set model  (inverse.init_1dmod, layer geometry, priors, bounds)
Build regularisation matrices  (inverse.diffops, scipy.sparse.block_diag)
Build ctrl_dict
Run inversion
  ├─ Parallel: joblib.Parallel → inverse.run_inv_flightline per file
  └─ Serial:   loop → inverse.run_inv_flightline per file
```

---

## Key Parameters

| Variable | Description |
|---|---|
| `AEM_system` | `'aem05'` or `'genesis'` |
| `RunType` | `'TikhOpt'` (Tikhonov with automatic τ selection) |
| `RegFun` | `'gcv'`, `'lcc'`, `'mle'`, `'fix'` |
| `Tau0` / `Tau1` | Regularisation parameter search grids |
| `RegShift` | Log-shift applied to optimal τ |
| `SetPrior` | `'update'` (rolling), `'set'` (uniform), `'read'` (from file) |
| `Nlyr` | Number of layers |
| `DzStart` / `DzEnd` | First and last layer thickness (m); log-spaced |
| `Guess_r` | Initial resistivity prior (Ohm·m) |
| `Guess_s` | Log-scale prior standard deviation |
| `Parallel` | `True` to use joblib parallelism |
| `Njobs` | Number of parallel workers (`-1` = all CPUs) |
| `Uncert` | If `True`, compute and store posterior covariance |

---

## `ctrl_dict` Structure

```python
ctrl_dict = {
    'system'    : [AEM_system, FwdCall],
    'header'    : [titstrng, ''],
    'inversion' : [RunType, RegFun, Tau0, Tau1, Maxiter, ThreshFit,
                   LinPars, SetPrior, Delta, RegShift],
    'covar'     : [L0, Cm0, L1, Cm1],
    'uncert'    : [Uncert],
    'data'      : [DataTrans, data_active, DatErr_add, DatErr_mult, Direction],
    'model'     : [ParaTrans, mod_act, mod_apr, mod_var, mod_bnd],
}
```

---

## Output Files

```
<OutResDir>/<file><outstrng>_results.npz   (from run_inv_flightline)
<OutResDir>/<file><outstrng>_ctrl.npz      (from run_inv_flightline)
```

The `outstrng` encodes key inversion parameters for traceability.

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `datetime`, `process_time`, `time`, `getpass` unused | Removed |

---

## Usage Notes

- `SearchStrng = '*data.npz'` searches `InDatDir` for matching files.
- Set `FileList = 'set'` and populate `dat_files` explicitly for single-file runs.
- `ThreshFit = [target_rms, dm_thresh, dd_thresh, 'rms'|'smp']` controls
  the convergence criterion inside `run_inv_flightline`.
- Compatible with jupytext (`.py` ↔ `.ipynb`).
