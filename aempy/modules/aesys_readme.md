# aesys_readme.md

## Module: `aesys.py`

**Purpose:** AEM system definitions and data I/O for AEMpyX.
Handles system parameter look-up, raw survey data ingestion,
and reading/writing of the internal AEMpyX data format in
ASCII, compressed NPZ, and NetCDF4.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

Created: Sun Nov 1 17:08:06 2020
Last change: vr Apr 2026

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `os`, `sys` |
| Third-party | `numpy`, `netCDF4` |
| Project-internal | `util` |

---

## Naming Conventions

- Input parameters use CapitalizedCamelCase (`File`, `System`, `OutInfo`).
- Internal variables use lowercase (`nn`, `ff`, `data_in`).
- Boolean output-control parameters follow the pattern `OutInfo`.

---

## Supported AEM Systems

| System key | Description | n_data |
|---|---|---|
| `aem05` | Tellus AEM05 FD system | 8 |
| `sglo` | AEM05 overlap variant a1/nm | 8 |
| `sglt` | AEM05 test-line variant DKV 2062 | 8 |
| `genesis` | Genesis TD system | 22 |
| `cggo` | Genesis overlap variant a1/nm | 22 |
| `cggt` | Genesis test-line variant | 22 |

---

## Function Reference

### System parameters

**`get_system_params(System, OutInfo)`**
Return forward-call string, column counts `nn`, ASCII format string,
column header string, and `miscpars` (frequencies/time-gates, units,
component dictionary) for the named AEM system.
Returns `(fwdcall, nn, fmt, col, miscpars)`.

`nn` layout: `[n_total, n_meta, n_data, n_optn]`
where meta = line + x + y + gps + alt + dem (6 columns).

---

### Survey data ingestion

**`read_survey_data(DatFile, Survey, OutInfo, Invalid, EPSG_in, EPSG_out)`**
Read raw XYZ survey deliveries and map columns to the internal AEMpyX
format. Supported survey codes: `fm`, `a1`–`a7`, `wf`, `tb`, `cav`/`cv`.
Applies coordinate reprojection via `util.project_utm_to_utm` where needed.
Returns `Data` (numpy array).

---

### Header utilities

**`get_header(file, headstr, OutInfo)`**
Extract comment lines beginning with `headstr` from a text file.
Returns `Header` (list of strings).

**`grow_header(Header, Addstr, OutInfo)`**
Append `Addstr` to an existing header list and normalise to a flat string.
Returns `Header` (str).

**`print_header(Header)`**
Print header to stdout.

---

### Unified read / write

**`read_aempy(File, Format, System, OutInfo)`**
Dispatch to the appropriate reader based on file extension
(`.npz`, `.asc`, `.nc4`).
Returns `(Data, Header, System)`.

**`write_aempy(File, Data, Format, System, Header, OutInfo)`**
Dispatch to the appropriate writer based on file extension or `Format`.

---

### Format-specific I/O

| Function | Format | Direction |
|---|---|---|
| `write_aempy_asc` | ASCII | write |
| `read_aempy_asc` | ASCII | read |
| `write_aempy_npz` | Compressed NPZ | write |
| `read_aempy_npz` | Compressed NPZ | read |
| `write_aempy_ncd` | NetCDF4 | write |
| `read_aempy_ncd` | NetCDF4 | read |

**`merge_data_files(File_list, Merged, MergedHeader, AEM_system, OutInfo)`**
Concatenate a list of NPZ data files into a single output file.
Returns merged `Data` array.

---

## Internal Data Format

```
col  0       line number
col  1       Easting  (m)
col  2       Northing (m)
col  3       GPS time / fiducial
col  4       radar altitude (m)
col  5       DEM / terrain height (m)
col  6 ..    data channels (n_data columns)
col  6+n_data .. optional channels (PLM, flags, ...)
```

---

## Known Issues / Notes

- `read_aempy_ncd`: variable access uses dict subscript `ncin.variables['key']`
  (not a function call). Fixed Apr 2026.
- `read_survey_data` `cav` branch: `numpy.float()` removed in Python 3.12;
  replaced with `float()`. Fixed Apr 2026.
- Smart/curly quotes (U+2018/2019) in `grow_header` caused Python 3.12
  parse failures; replaced with ASCII apostrophes. Fixed Apr 2026.
- `compdict` entry for `'ALT'` was malformed in `sglo` and `sglt` system
  definitions; corrected to `[col, value, label]` ordering. Fixed Apr 2026.

---

## Change Log

| Date | Author | Note |
|---|---|---|
| Nov 2020 | vrath | Initial version (AEMpyX). |
| Jun 2021 | duygu | Edits and extensions. |
| Apr 2023 | duygu | Further edits. |
| Apr 2026 | vrath | Standardised docstrings; capitalised parameters; provenance; readme. Bug fixes (see above). |
