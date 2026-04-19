# PROJECT_inv_dataset_genesis_readme.md

## Script: `PROJECT_inv_dataset_genesis.py`

**Purpose:** AEMpyX dataset inversion for the Genesis TEM system.
Runs a self-contained per-site inversion loop (not via
`run_inv_flightline`) supporting `TikhOpt`, `TikhOcc`, and `MAP`
run types. Results are written per flight line.

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
| Third-party | `numpy`, `scipy` |
| AEMpyX modules | `version`, `aesys`, `util`, `inverse` |

---

## Workflow

```
Set system (AEM_system = 'genesis')
Set file list  (util.get_data_list or explicit)
Set regularisation  (RegFun, Tau0, Tau1)
Set model  (inverse.init_1dmod, layer geometry, priors, bounds)
Build ctrl_dict  (based on RunType)
Loop over flight-line files
  └─ aesys.read_aempy → DataObs
  └─ loop over sites
       ├─ set Model dict, Data dict
       ├─ inverse.run_tikh_opt / run_tikh_occ / run_map
       └─ accumulate site results
  └─ numpy.savez_compressed → _results.npz
```

---

## Supported Run Types

| `RunType` | Method |
|---|---|
| `'TikhOpt'` | Tikhonov with automatic regularisation parameter selection |
| `'TikhOcc'` | Occam smoothest-model inversion |
| `'MAP'` | Maximum a-posteriori with explicit spatial covariance |

---

## Key Parameters

| Variable | Description |
|---|---|
| `AEM_system` | `'genesis'` (or `'aem05'` with adjusted settings) |
| `RunType` | See table above |
| `RegFun` | `'gcv'`, `'lcc'`, `'mle'`, `'fix'` |
| `SetPrior` | `'set'` (uniform), `'update'` (rolling), `'read'` (file) |
| `Nlyr` | Number of layers |
| `data_active[0:11] = 0` | Disables horizontal component (Genesis-specific) |
| `Uncert` | If `True`, store posterior Jacobian and covariance |
| `Direction` | `'normal'` or `'reverse'` site-loop order |

---

## `ctrl_dict` Structure (TikhOpt)

```python
ctrl_dict = {
    'system'    : [AEM_system, FwdCall],
    'header'    : [titstrng, ''],
    'inversion' : [RunType, RegFun, Tau0, Tau1, Maxiter, ThreshRMS,
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
<OutResDir>/<name><outstrng>_ctrl_results.npz   (ctrl dict)
<OutResDir>/<name><outstrng>_results.npz
    fl_data, fl_name, header
    mod_ref, mod_act, dat_act
    site_modl, site_sens, site_merr
    site_dobs, site_dcal, site_derr
    site_nrms, site_num, site_log
    site_x, site_y, site_gps, site_alt, site_dem
    [site_jacd, site_pcov]   if Uncert=True
```

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `copy` unused | Removed |
| `alg.run_*` calls | `alg` module does not exist in AEMpyX | Changed to `inverse.run_*` |
| `mod_bnd == None` | Deprecated identity comparison | Changed to `is not None` |
| `mod_bnd` condition | Logic inverted (`not ... or`) | Fixed to `is not None and size != 0` |
| `block_diag([Cmi for Cmi in range(7)])` | Loop variable shadows outer `Cmi` | Changed to `for _ in range(7)` |
| `block_diag([CmS for CmS in range(7)])` | Same issue for `CmS` | Changed to `for _ in range(7)` |
| Duplicate `start = time.time()` | Reset timer inside loop before file open | Removed early duplicate |
| `ctrl_dict` | Missing `'data'`, `'model'`, `'header'` keys required by `run_*` | Added all required keys |
| `Uncert` result save | Both branches wrote to same `_results.npz` overwriting first save | Second save (with Uncert) supersedes; first is now omitted |

---

## Usage Notes

- This script manages its own site loop — it does **not** call
  `run_inv_flightline`. Per-site `Model` and `Data` dicts are built
  explicitly and passed directly to `inverse.run_tikh_opt` etc.
- The MAP block requires `inverse.set_zcenters` and `inverse.covar`;
  it is a rho-only workaround for now.
- Compatible with jupytext (`.py` ↔ `.ipynb`).
