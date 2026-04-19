# PROJECT_viz_site_uncert_readme.md

## Script: `PROJECT_viz_site_uncert.py`

**Purpose:** AEMpyX site-level uncertainty visualisation.
Plots posterior model and data uncertainty (from JCN, RTO, or
covariance outputs) for individual inversion sites.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `inspect`, `getpass`, `random` |
| Third-party | `numpy`, `matplotlib`, `scipy.interpolate` |
| AEMpyX modules | `version`, `util`, `aesys`, `inverse`, `viz` |

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `functools`, `cycler`, `scipy.linalg` unused | Removed |
