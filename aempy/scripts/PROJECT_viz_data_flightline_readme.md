# PROJECT_viz_data_flightline_readme.md

## Script: `PROJECT_viz_data_flightline.py`

**Purpose:** AEMpyX flight-line data visualisation.
Plots observed (and optionally calculated) AEM data per flight line,
with optional PDF catalogue output.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys`, `inspect`, `getpass` |
| Third-party | `numpy`, `matplotlib` (pyplot, backends.backend_pdf), `cycler` |
| AEMpyX modules | `version`, `util`, `aesys`, `viz` |

---

## Workflow

```
Set system, plot parameters
Set file list (util.get_data_list)
Optional: open PDF catalogue (PdfPages)
Loop over flight-line files
  └─ aesys.read_aempy → Data
  └─ viz.plot_flightline_aem05  or  viz.plot_flightline_genesis
  └─ matplotlib.pyplot.savefig  (PlotFmt list)
  └─ catalog.savefig  if PDFCatalog
Close PDF catalogue
```

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Imports | `process_time`, `warnings`, `prep` unused | Removed |
