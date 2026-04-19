# inverse_readme.md

## Module: `inverse.py`

**Purpose:** 1-D AEM inversion engine for AEMpyX.
Provides regularised least-squares inversion (Tikhonov variants),
uncertainty estimation, model and data transforms, prior management,
and flight-line / area batch processing.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

Created: Tue Nov 26 2024
Last change: vr Apr 2026

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `time`, `datetime`, `warnings`, `copy` |
| Third-party | `numpy`, `scipy` (`linalg`, `sparse`, `special`, `fftpack`), `numba` (`njit`, `prange`), `joblib` |
| Project-internal | `aesys`, `core1d`, `util` |

---

## Naming Conventions

- Major input parameter names use CapitalizedCamelCase (`Ctrl`, `Model`, `Data`).
- Short dict keys and loop variables remain lowercase (`mod_apr`, `dat_obs`).
- Boolean flags follow the pattern `out` / `Out*` (`OutInfo`, `out`).

---

## `ctrl` Dictionary Structure

```python
ctrl = {
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

### `SetPrior` options

| Value | Behaviour |
|---|---|
| `'set'` | Use `mod_apr` uniformly for all sites. |
| `'read'` | Load site-specific priors from `prior_file`. |
| `'upd'` | Rolling prior: update from previous site's result. |

---

## Top-level Workflow

**`run_inv_flightline(data_dir, data_file, ctrl, prior_file, result_dir, result_strng, results_out, out)`**
Batch inversion of all sites in an AEM flight line or area.
Reads data via `aesys.read_aempy`; writes compressed NPZ result files.
Dispatches per-site to `run_tikh_opt`, `run_tikh_lcurve`, `run_tikh_doi`,
`run_tikh_occ`, or `run_map` depending on `ctrl['inversion'][0]`.

---

## Selected Function Reference

### Inversion runners (single site)

**`run_tikh_opt(Ctrl, Model, Data, OutInfo)`**
Tikhonov inversion with L-curve or GCV regularisation parameter selection.

**`run_tikh_lcurve(Ctrl, Model, Data, OutInfo)`**
Tikhonov inversion with explicit L-curve search.

**`run_tikh_doi(Ctrl, Model, Data, OutInfo)`**
Tikhonov inversion with depth-of-investigation analysis.

**`run_tikh_occ(Ctrl, Model, Data, OutInfo)`**
Occam-style smoothest-model inversion.

**`run_map(Ctrl, Model, Data, OutInfo)`**
Maximum a-posteriori (MAP) inversion.

### Error and transform utilities

**`set_errors(dat_obs, daterr_add, daterr_mult)`**
Compute data error vector from additive and multiplicative components.
Returns `(dat_err, comment)`.

**`transform_data(d_vec, e_vec, d_trn, d_state)`**
Apply or reverse data transformation (log, linear, etc.).
Returns `(d_vec, e_vec, d_state)`.

**`transform_parameter(m_vec, m_trn, m_state, mode)`**
Apply or reverse model parameter transformation.
Returns `(m_vec, m_state)`.

### Prior management

**`load_prior(prior_file, m_ref, m_apr, m_act)`**
Load site priors from NPZ file and map to active parameter vector.
Returns `site_prior`.

**`set_prior(pval, flightline)`**
Broadcast a prior value array to match the shape of `flightline` data.
Returns `prior`.

**`insert_mod(M, m, m_act)`**
Insert active-parameter subset `m` back into full model vector `M`.
Returns `M`.

### Model utilities

**`get_nlyr(mod_ref)`**
Extract number of layers from reference model vector.
Returns `nlyr`.

**`init_layers(nlyr, start, end, logspace, out)`**
Generate layer thicknesses, node coordinates, and cell-centre depths.
Returns `(dz, z_node, z_cent)`.

**`calc_mad(data, center)`**
Median absolute deviation of `data` about `center`.

### Dataset utilities

**`merge_data_sets(infile_list, outfile_name, aem_system, dictout, out)`**
Concatenate multiple AEM data NPZ files into one.
Returns `merged_data` (array or dict).

**`merge_model_sets(infile_list, outfile_name, qthresh, out)`**
Concatenate inversion result NPZ files with optional quality filtering.
`qthresh = ['smp'|'rms', threshold_value]`.
Returns `merged_models` (dict).

---

## Output Files

| File | Contents |
|---|---|
| `*_results.npz` | `site_modl`, `site_sens`, `site_nrms`, `site_smap`, `site_conv`, `site_x/y/gps/alt/dem`, `site_pcov`, `mod_ref`, `mod_act` |
| `*_ctrl.npz` | Serialised `ctrl` dictionary for full reproducibility |

---

## Known Issues / Notes

- `init_layers`: `logspace=False` branch corrected from
  `numpy.linalg.space` (does not exist) to `numpy.linspace`. (Apr 2026)
- `set_prior`: bare `error()` call replaced with `sys.exit()`;
  `numpy.size()` with no argument replaced with `numpy.size(pval)`. (Apr 2026)
- `merge_model_sets`: `'smp'` and `'rms'` quality filters are evaluated
  in separate `if`-blocks — a site failing both criteria may be appended
  twice. Use mutually exclusive criteria or filter post-merge if needed.
- Numba-compiled (`@njit`) functions incur a first-call compilation cost.
- Unescaped apostrophe in `sys.exit('...don't...')` fixed. (Apr 2026)
- Invalid escape sequence `\&` in a docstring reference fixed. (Apr 2026)

---

## Change Log

| Date | Author | Note |
|---|---|---|
| Nov 2020 | vrath | Initial version (AEMpyX). |
| Jun 2024 | vrath | Bloomsday refactor; `run_inv_flightline` restructured. |
| Nov 2024 | vrath | Module header dated; parallel `joblib` path added. |
| Apr 2026 | vrath | Bug fixes; unused imports removed; provenance; readme. |
