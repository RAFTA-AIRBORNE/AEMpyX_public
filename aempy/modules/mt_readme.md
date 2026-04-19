# mt_readme.md

## Module: `mt.py`

**Purpose:** Magnetotelluric (MT) processing utilities for AEMpyX.
Covers data and model I/O in ASCII, NetCDF4, and Fortran binary formats;
1-D MT forward modelling; 3-D model construction and manipulation;
and spatial filtering identical to the routines in `post.py`.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

Last change: vr Apr 2026

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `sys` |
| Third-party | `numpy`, `scipy.io` (`FortranFile`), `scipy.ndimage` (`convolve`, `uniform_filter`, `gaussian_filter`, `median_filter`), `netCDF4` |

---

## Function Reference

### Data I/O

**`decode_h2(strng)`**
Decode a two-character hex string to an integer (used for header parsing).

**`read_jac(JacFile, OutInfo)`**
Read a Jacobian matrix from a Fortran binary file.
Returns `(Jac, nsite, ncomp, nparam)`.

**`read_data_jac(DatFile, OutInfo)`**
Read MT data and associated Jacobian from file.
Returns `(Dat, Jac, Site, Comp, Head)`.

**`read_data(DatFile, OutInfo)`**
Read MT data from ASCII file.
Returns `(Dat, Site, Comp, Head)`.

**`write_data(DatFile, Dat, Site, Comp, Head, OutInfo)`**
Write MT data to ASCII file.

**`write_jac_ncd(NCFile, Jac, Dat, Site, Comp, OutInfo)`**
Write Jacobian and data to NetCDF4 file.

**`write_data_ncd(...)`**
Write MT data to NetCDF4 file.

**`write_model_ncd(...)`**
Write 3-D resistivity model to NetCDF4 file.

**`write_model(ModFile, dx, dy, dz, mval, reference, OutInfo)`**
Write 3-D model to ASCII file.

**`read_model(ModFile, trans, volumes, OutInfo)`**
Read 3-D model from ASCII file.
Returns `(x, y, z, rho [, vol])`.

---

### 1-D MT forward modelling

**`mt1dfwd(freq, sig, d, inmod, outdat)`** *(two implementations)*
1-D MT forward model using the propagator-matrix method.
`inmod`: `'res'` or `'r'` (resistivity) / `'sig'` (conductivity).
`outdat` / `out`: `'both'` | `'imp'` | `'rho'` | `'pha'`.
Returns apparent resistivity and/or phase, or impedance Z.

**`calc_rhoa_phas(freq, Z)`**
Convert impedance tensor `Z` to apparent resistivity and phase.
Returns `(rhoa, phase)`.

---

### 3-D model construction

**`insert_body(rho, dx, dy, dz, reference, body_rho, body_shape, body_params, Out)`**
Insert a resistivity anomaly of defined shape into a background model.

**`cells3d(dx, dy, dz, center, reference)`**
Compute 3-D cell centre or corner coordinates from cell-size vectors.
Returns `(x, y, z)`.

**`in_ellipsoid(x, y, z, params)`**
Return boolean mask: True where points lie inside the defined ellipsoid.

**`in_box(x, y, z, params)`**
Return boolean mask: True where points lie inside the defined box.

**`clip_model(x, y, z, rho, ...)`**
Clip a 3-D model to spatial bounds.

**`linear_interpolation(p1, p2, x0)`**
Linearly interpolate between two points at coordinate `x0`.

---

### Rotation matrices

**`rotx(theta)` / `roty(theta)` / `rotz(theta)`**
Return 3×3 rotation matrices about the x, y, z axes respectively.
`theta` in degrees.

---

### Spatial filtering (mirrors `post.py`)

**`medfilt3D(M, kernel_size, boundary_mode, maxiter, Out)`**
Iterated median filter on an n-D array.

**`anidiff3D(M, ckappa, dgamma, foption, maxiter, Out)`**
Convenience wrapper for anisotropic diffusion.

**`anisodiff3D(stack, niter, kappa, gamma, step, option)`**
3-D Perona–Malik anisotropic diffusion.

**`shock3d(M, dt, maxiter, filt, boundary_mode, signfunc)`**
Shock filter (morphological sharpening) in 3-D.

**`gauss3D(Kshape, Ksigma)`**
Normalised 3-D Gaussian kernel.

---

### Model preparation

**`prepare_model(rho, rhoair)`**
Clip or pad a resistivity model to ensure the air layer has
resistivity `rhoair`. Returns `rho`.

---

## Known Issues / Notes

- `mt1dfwd` is defined twice with slightly different signatures
  (`outdat` vs `out`, `inmod='res'` vs `inmod='r'`). The second
  definition overrides the first at import time.
- Removed unused imports: `os`, `inspect`, `norm` (numpy.linalg),
  `laplace`, `gaussian_gradient_magnitude`, local `sleep`. (Apr 2026)

---

## Change Log

| Date | Author | Note |
|---|---|---|
| 2020–2021 | vrath / duygu | Initial version (AEMpyX). |
| 2023 | vrath | Extended with shock filter and model utilities. |
| Apr 2026 | vrath | Removed unused imports; provenance; readme. |
