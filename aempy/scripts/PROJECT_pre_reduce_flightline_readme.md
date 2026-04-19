# PROJECT_pre_reduce_flightline_readme.md

## Script: `PROJECT_pre_reduce_flightline.py`

**Purpose:** AEMpyX flight-line data reduction.
Supports three actions on processed flight-line NPZ files:
decimation (block-averaging), splitting into interleaved subsets,
and profile-distance cutting.

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
| AEMpyX modules | `version`, `util`, `prep`, `aesys` |

---

## Actions

| `Action` | Description |
|---|---|
| `'decimate'` | Block-average with `prep.reduce_data(Method=['mean'|'median', Window])` |
| `'split'` | Interleave: extract every `Step`-th sample starting at 0, 1, …, Step-1 |
| `'cut'` | Trim to a profile-distance interval `[start_m, end_m]` (single file only) |

---

## Bugs Fixed (Apr 2026)

| Location | Issue | Fix |
|---|---|---|
| Path fix | `InDatDir = OutDatDir + '/'` (wrong variable) | `OutDatDir = OutDatDir + '/'` |
| `prep.get_profile` | Function does not exist in `prep` | Replaced with inline numpy cumulative distance: `cumsum(sqrt(dx²+dy²))` |
| `prep.reduce_data` | Positional arg → `DataVec=D` explicit keyword | Fixed |
| Imports | `warnings`, `getpass`, `inverse` unused | Removed |
