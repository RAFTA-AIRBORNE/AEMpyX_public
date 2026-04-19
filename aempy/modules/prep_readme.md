# prep_readme.md

## Module: `prep.py`

**Purpose:** Data preparation utilities for AEM data processing.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

Created: Sat Nov 21 16:05:35 2020
Last change: vr Apr 2026

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `sys` |
| Third-party | `numpy`, `scipy` (`signal`, `linalg`, `interpolate`) |
| Project-internal | `aesys` |

---

## Naming Conventions

- Input parameters use CapitalizedCamelCase (e.g. `Columns`, `Method`, `OutInfo`).
- Internal / loop variables use lowercase (e.g. `cols`, `meth`, `blocksize`).
- Boolean output-control parameters follow the pattern `Out*` (`OutInfo`, `OutFull`).

---

## Function Reference

**`filter_column(M, Columns, Method, OutInfo)`**
Apply a digital filter to selected columns of a data matrix.
`Method = ['butter', order, cutoff_freq]`.
Returns `(M, headstring)`.

**`insert_flag(M, Criterion, ThreshVal, Columns, Flag, InCol, System)`**
Replace bad or out-of-range values in selected columns with `Flag` (default `numpy.nan`).
`Criterion`: `'neg'` | `'less'` | `'great'` | `'plm'` | `'nan'`.
Returns `(M, NaNindex)`.

**`handle_gaps(M, Columns, Impute, System)`**
Impute or remove NaN-flagged entries.
`Impute[0]`: `'nan'` | `'delete'` | `'ave'` | `'med'` | `'noi'` | `'spl'` | `'spln'`.
Returns `M`.

**`calc_svd_decomp(D, Columns, ErrD, K, ThreshMSE, OutFull, OutInfo)`**
PCA/SVD-based filter after Minsley et al. (2012).
Returns `P` — or `(P, U, S, V, MSE, FRO)` when `OutFull=True`.

**`reduce_data(DataVec, System, Method, ErrEst, OutInfo)`**
Block-average a multi-column array column by column.
`Method = ['mean'|'med'|'dec', blocksize]`.
Returns `(dat_out [, err_out], comment)`.

**`process_column(DataVec, Method, ErrEst, OutInfo)`**
Block-average a single 1-D vector.
Returns `(dat_out, dat_err, comment)`.

---

## Usage Notes

- All functions exit via `sys.exit()` on fatal parameter errors.
- Blocksize is silently forced to the next odd integer when even.
- `'spl'` gap-filling requires columns 1 and 2 of `M` to hold X and Y
  coordinates respectively.
- `calc_svd_decomp` `ThreshMSE` parameter is reserved and not yet active.

---

## Change Log

| Date | Author | Note |
|---|---|---|
| Nov 2020 | vrath | Initial version (AEMpyX). |
| Jul 2021 | vrath | `insert_flag`: added `'plm'` and `'nan'` criteria. |
| Sep 2021 | vrath | `handle_gaps`: added `'spl'`/`'spln'` method. |
| Feb 2023 | vrath | `reduce_data`: `ErrEst` option; odd-blocksize enforcement. |
| Apr 2026 | vrath | Standardised docstrings; capitalised parameters; provenance; readme. Bug fixes. |
