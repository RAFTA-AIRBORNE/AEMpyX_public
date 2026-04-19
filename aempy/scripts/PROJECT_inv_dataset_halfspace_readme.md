# PROJECT_inv_dataset_halfspace_readme.md

## Script: `PROJECT_inv_dataset_halfspace.py`

**Purpose:** AEMpyX halfspace (single-layer) inversion.
Inverts each site independently as a uniform halfspace, producing
fast, coarse estimates used as priors for multi-layer inversions.
Runs its own site loop (not via `run_inv_flightline`).

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `inspect`, `warnings` |
| Third-party | `numpy`, `scipy`, `time` (`process_time`) |
| AEMpyX modules | `version`, `aesys`, `util`, `inverse` |

---

## Workflow

```
Set system, errors, active channels
Set file list  (util.get_data_list or explicit)
Build Nlyr=1 model (halfspace), ctrl_dict
Loop over flight-line files
  └─ aesys.read_aempy → DataObs
  └─ loop over sites
       ├─ build Model, Data dicts
       ├─ inverse.run_tikh_opt
       └─ accumulate site results
  └─ numpy.savez_compressed → _halfspace.npz
```

---

## Key Parameters

| Variable | Description |
|---|---|
| `Nlyr` | Fixed at 1 (halfspace) |
| `RunType` | `'TikhOpt'` with `RegFun='fix'` (minimal regularisation) |
| `Guess_r` | Prior halfspace resistivity (Ohm·m) |
| `OutStrng` | `'_halfspace'` — appended to each output filename |

---

## Output Files

```
<OutResDir>/<name>_halfspace_ctrl.npz
<OutResDir>/<name>_halfspace.npz
    fl_data, fl_name, header, aem_system
    mod_ref, mod_act, dat_act
    site_modl, site_sens, site_merr
    site_dobs, site_dcal, site_derr
    site_nrms, site_num, site_log
    site_x, site_y, site_gps, site_alt, site_dem
```

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `copy`, `getpass`, `datetime`, `process_time` unused | Removed; `process_time` re-added correctly |
| `mod_bnd == None` | Deprecated comparison; inverted logic | `is not None and size != 0` |
| Duplicate `start = process_time()` | Timer reset before file open, then again inside loop | Removed first occurrence |
| `header=numpy.array(...),` | Trailing comma made `header` a tuple | Comma removed |
| `# %logstart -o` | Jupyter magic at module level | Removed |
