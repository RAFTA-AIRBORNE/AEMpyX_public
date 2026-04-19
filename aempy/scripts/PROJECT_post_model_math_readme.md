# PROJECT_post_model_math_readme.md

## Script: `PROJECT_post_model_math.py`

**Purpose:** Post-processing arithmetic on AEMpyX inversion results.
Currently implements `'twoway'` averaging: loads normal-direction and
reverse-direction inversions of the same flight line, averages the
log-transformed models, and computes the absolute difference.

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
| AEMpyX modules | `version`, `aesys`, `util` |

---

## Method: `'twoway'`

For each result file `*_normal*results.npz`, the script expects a
companion file `*_reverse*results.npz` in the same directory.

```
m_avg = 0.5 * (log(m_normal) + log(m_reverse))   → back-transformed
m_dif = |log(m_normal) - log(m_reverse)|
```

The averaged model replaces `site_modl`; `site_diff` is added as a
new key. All other fields are copied from the normal-direction file.

---

## Output Files

```
<OutResDir>/<n>_average.npz
    All keys from normal-direction result file
    site_modl  : geometric-mean model
    site_diff  : absolute log-difference between directions
    header     : updated with date
```

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `from datetime import datetime` duplicated twice | Deduplicated |
| Imports | `import getpass` duplicated twice | Deduplicated |
| Imports | `process_time`, `time`, `warnings`, `scipy`, `inverse` unused | Removed |
| `m_dif` | `numpy.abs(m_normal + m_revers)` — sum instead of difference | Fixed to `m_normal - m_revers` |
| `numpy.exp(m)` | Applied to bare `m` (undefined at that point) | Changed to `numpy.exp(m_avg)` |
