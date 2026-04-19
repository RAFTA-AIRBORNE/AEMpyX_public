# PROJECT_inv_dataset_rto_readme.md

## Script: `PROJECT_inv_dataset_rto.py`

**Purpose:** AEMpyX randomize-then-optimize (RTO) uncertainty estimation.
Runs Tikhonov inversion with RTO sampling on a selected subset of sites
to generate a posterior ensemble and derive uncertainty statistics
(mean, variance, median, MAD, percentiles).

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

## Site Selection (`Sample`)

| `Sample[0]` | Behaviour |
|---|---|
| `'rand'` | Random `Num_samples` sites |
| `'step'` | Every `Step`-th site from `Start` to `Stop` |
| `'list pos'` | Explicit list of site indices (`Samplist`) |
| `'list dist'` | Sites nearest to given along-profile distances (`Distlist`) |

---

## Workflow

```
Set system, errors, data_active
Set file list
Set Sample strategy → site_list
Build ctrl_dict (TikhOpt)
  └─ Ctrl['rto']    = [NSamples, Percentiles]
  └─ Ctrl['output'] = ['ens']  (to store full ensemble)
Loop over files
  └─ aesys.read_aempy → DataObs
  └─ compute site_r (along-profile distance from origin)
  └─ build site_list from Sample
  └─ loop over site_list  (num_site counter)
       ├─ inverse.run_rto → results
       └─ accumulate: site_modl, site_jacd,
                      rto_avg/var/med/mad/prc [/ens]
  └─ numpy.savez_compressed → <Fileout>.npz
  └─ util.add_object_npz → append site_rto_ens if requested
```

---

## RTO Parameters

| Variable | Description |
|---|---|
| `NSamples` | Number of RTO samples per site |
| `Percentiles` | Quantiles to compute (linear scale) |
| `Ctrl['rto']` | `[NSamples, Percentiles]` — consumed by `inverse.run_rto` |
| `Ctrl['output']` | `['ens']` to store full ensemble in output NPZ |

---

## Output Keys

```
fl_data, fl_name, fl_orig, header
mod_ref, mod_act, dat_act
site_modl, site_sens, site_merr
site_dobs, site_dcal, site_derr
site_nrms, site_num, site_log
site_jacd, site_pcov
site_rto_avg, site_rto_var, site_rto_med, site_rto_mad, site_rto_prc
[site_rto_ens]   if 'ens' in Ctrl['output']
```

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `copy`, `jit`, `prange` unused | Removed |
| `import alg` | `alg` module does not exist | Removed; `alg.run_rto` → `inverse.run_rto` |
| `mod_bnd == None` | Deprecated comparison; inverted logic | `is not None and size != 0` |
| `Ctrl['output'] = ['ens ']` | Trailing space in key string | Removed space → `'ens'` |
| `ctrl_dict` | Missing `'data'`, `'model'`, `'header'`; used old `'transform'` | Standardised to full key set |
| `# %logstart -o` | Jupyter magic at module level | Removed |
