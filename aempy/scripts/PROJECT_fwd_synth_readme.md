# PROJECT_fwd_synth_readme.md

## Script: `PROJECT_fwd_synth.py`

**Purpose:** AEMpyX 1-D forward modelling.
Computes synthetic AEM responses for a parameterised model, with
optional perturbation of responses to generate data ensembles for
subsequent inversion studies.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `sys`, `os`, `inspect` |
| Third-party | `numpy` |
| AEMpyX modules | `version`, `util`, `inverse`, `aesys` |

---

## Workflow

```
Set system (AEM_system)
  └─ aesys.get_system_params  → FwdCall, NN, ...
Set model (nlyr, Model_base, VarPar, VarInd)
Loop over parameter values (VarPar)
  └─ inverse.transform_parameter  → m_current
  └─ inverse.calc_fwdmodel        → d_ref  (noiseless response)
  └─ inverse.set_errors × Nsamples → data_obs  (perturbed ensemble)
Save results as compressed NPZ
  └─ SplitData=True  → one file per model variant
  └─ SplitData=False → single file
```

---

## Key Parameters

| Variable | Description |
|---|---|
| `AEM_system` | `'aem05'` or `'genesis'` |
| `nlyr` | Number of layers |
| `Model_base` | Base model vector (resistivity, thicknesses, …) |
| `VarPar` | List of values to iterate for the variable parameter |
| `VarInd` | Index in the model vector of the iterated parameter |
| `Alt` | List of altitudes; use `VarInd == numpy.size(m_i)` to iterate altitude |
| `Nsamples` | Number of perturbed realisations per model variant |
| `PerturbDat` | If `True`, add noise to each sample |
| `SplitData` | If `True`, write one NPZ per model variant |
| `DataTrans` | 0 = raw, 1 = ln, 2 = asinh |
| `DatErr_add` | Additive error floor (ppm) |
| `DatErr_mult` | Multiplicative error fraction |

---

## Output Files

```
<OutDir>/<FWDBaseName>[_model<i>_<N>samples].npz
    model   : model vector for this variant
    data    : (1 + Nsamples) × (3 + ndata) array
              col 0: model index
              col 1: sample index (-1 = noiseless reference)
              col 2: altitude
              col 3+: data channels
    para    : parameter record [mod_num, VarInd, VarVal, DataTrans,
                                 DatErr_add, DatErr_mult]
```

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| `calc_fwdmodel` call | `alt=Alt` (list) passed instead of `alt=alt` (scalar) | Changed to `alt=alt` |
| `set_errors` call | Positional args; now uses keyword args | `daterr_add=`, `daterr_mult=` |
| Imports | `time`, `warnings`, `process_time`, `datetime`, `getpass`, `core1d` unused | Removed |

---

## Usage Notes

- Altitude loop: set `VarInd = numpy.size(m_i)` and `VarPar = Alt` to sweep altitude
  rather than a model parameter.
- `DataTrans = 1` requires strictly positive data; use `DataTrans = 2` (asinh)
  for TDEM data that may contain sign changes.
- Compatible with jupytext (`.py` ↔ `.ipynb` conversion via `formats: py:light,ipynb`).
