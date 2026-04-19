# PROJECT_viz_models_flightline_readme.md

## Script: `PROJECT_viz_models_flightline.py`

**Purpose:** AEMpyX flight-line model visualisation.
Plots 1-D resistivity model pseudo-sections and data fits along
individual flight lines, with optional PDF catalogue output.

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
| Third-party | `numpy`, `matplotlib` (pyplot, colors, cm, backends.backend_pdf), `datetime` |
| AEMpyX modules | `version`, `util`, `aesys`, `inverse` |

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `viz` imported but only appears in a commented-out line | Removed |
