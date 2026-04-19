# PROJECT_viz_synth_readme.md

## Script: `PROJECT_viz_synth.py`

**Purpose:** AEMpyX synthetic inversion ensemble visualisation.
Loads results from `PROJECT_inv_synth.py` and plots model and data
ensembles (percentile bands + true model overlay) as paired figures.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `inspect`, `getpass`, `warnings` |
| Third-party | `numpy`, `matplotlib` (pyplot, colors, cm, collections, patches, backends.backend_pdf) |
| AEMpyX modules | `version`, `aesys`, `util`, `viz`, `inverse` |

---

## Workflow

```
Set directories, plot parameters
util.get_filelist → data_files
Loop over result files
  └─ numpy.load → m_act, m_ref, m_ens, d_ens, r_ens
  └─ inverse.calc_stat_ens → ens_modl, ens_dcal, ens_nrms
  └─ optional: load m_true, d_true
  └─ matplotlib subplots (1 × 2)
       ax[0]: viz.plot_model_ensemble + viz.plot_model (true)
       ax[1]: viz.plot_data_ensemble  + viz.plot_data  (true)
  └─ savefig (PlotFormat list)
Optional PDF catalogue
```

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| `import viz` | Duplicated twice | Removed duplicate |
