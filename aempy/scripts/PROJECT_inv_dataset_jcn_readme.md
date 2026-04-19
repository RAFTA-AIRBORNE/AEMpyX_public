# PROJECT_inv_dataset_jcn_readme.md

## Script: `PROJECT_inv_dataset_jcn.py`

**Purpose:** AEMpyX jackknife (JCN) uncertainty estimation.
Runs Tikhonov inversion with jackknife resampling on a selected
subset of sites to estimate model uncertainty from data perturbation.

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
Build ctrl_dict (TikhOpt + jackknife output)
  └─ Ctrl['output'] = ['ens'] to store ensemble
Loop over files
  └─ aesys.read_aempy → DataObs
  └─ compute site_r (along-profile distance from origin)
  └─ build site_list from Sample
  └─ loop over site_list
       ├─ inverse.run_jcn → results
       └─ accumulate: site_modl, site_jacd, jcn_avg/var/med/mad [/ens]
  └─ numpy.savez_compressed → <Fileout>.npz
  └─ util.add_object_npz → append site_jcn_ens if requested
```

---

## Output Keys

```
fl_data, fl_name, fl_orig, header
mod_ref, mod_act, dat_act
site_modl, site_sens, site_merr
site_dobs, site_dcal, site_derr
site_nrms, site_num, site_log
site_jacd, site_pcov
site_jcn_avg, site_jcn_var, site_jcn_med, site_jcn_mad
[site_jcn_ens]   if 'ens' in Ctrl['output']
```

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `copy`, `random` unused | Removed |
| `import alg` | `alg` module does not exist | Removed; `alg.run_jcn` → `inverse.run_jcn` |
| `mod_bnd == None` | Deprecated comparison; inverted logic | `is not None and size != 0` |
| `block_diag([Cmi for Cmi in range(7)])` | Loop variable shadows outer `Cmi` | `for _ in range(7)` |
| `block_diag([CmS for CmS in range(7)])` | Same issue | `for _ in range(7)` |
| `site_num = numpy.array([])` inside loop | Reset every iteration — always took first-site branch | Moved to first-site check via `if ii == site_list[0]` |
| `ctrl_dict` | Missing `'data'`, `'model'`, `'header'` keys; used old `'transform'` | Standardised to full key set |
| `# %logstart -o` | Jupyter magic at module level | Removed |
