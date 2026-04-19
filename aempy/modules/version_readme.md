# version_readme.md

## Module: `version.py`

**Purpose:** Version string management for AEMpyX.
Single source of truth for the package version number and release
timestamp, consumed by `util.print_title` and run-control scripts.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

Created: Feb 2021
Last change: vr Apr 2026

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `sys`, `os`, `inspect`, `datetime` |

Note: `sys`, `os`, and `inspect` are used solely for the
`sys.path.append` call that ensures the module is importable
from its own directory. `datetime` is used in `versionstrg`.

---

## Function Reference

**`versionstrg()`**
Return the current version string and a formatted release timestamp.

```python
version, release_date = versionstrg()
# version      -> '0.99.99'
# release_date -> 'MM/DD/YYYY, HH:MM:SS'  (current time at call)
```

---

## Usage Notes

- The `sys.path.append` at module level modifies the interpreter path
  at import time. This is an established pattern in AEMpyX run scripts;
  be aware of potential side-effects when importing in other contexts.
- `release_date` reflects the time of the function call, not a fixed
  release date. For reproducible logging, call `versionstrg()` once at
  the start of a run and store the result.

---

## Change Log

| Date | Author | Note |
|---|---|---|
| Feb 2021 | vrath | Initial version (AEMpyX). |
| Apr 2026 | vrath | Provenance block added; readme generated. |
