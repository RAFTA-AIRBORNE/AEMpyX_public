# post_readme.md

## Module: `post.py`

**Purpose:** Post-processing utilities for AEM inversion results.
Provides spatial filtering, anisotropic diffusion, cross-gradient
coupling, and quality-control selection of inversion sites.

---

## Provenance

AEMpyX project.

**Authors:** Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

Created: Apr 4, 2021
Last change: vr Apr 2026

---

## Dependencies

| Type | Package |
|---|---|
| Standard library | `sys` |
| Third-party | `numpy`, `scipy.linalg`, `scipy.sparse`, `scipy.signal`, `scipy.ndimage` (`uniform_filter`, `gaussian_filter`, `median_filter`), `differint` |

---

## Naming Conventions

- Input parameters use CapitalizedCamelCase where standardised.
- Boolean output-control parameters follow the pattern `Out` / `OutInfo`.

---

## Function Reference

**`fractrans(m, x, a)`**
Compute the fractional derivative of order `a` of vector `m`
over domain `x`, using the `differint` package.
Returns `mm` (numpy array).

**`crossgrad(m1, m2, mesh, Out)`**
Compute the cross-gradient coupling functional between two model
vectors `m1` and `m2` on a shared `mesh`.
Returns `(cgm, cgnm)` — the cross-gradient vector and its
normalised form. 1-D input is rejected; 2-D yields a scalar
cross-product, 3-D yields a 3-component vector.

Reference: Rosenkjaer et al. (2015) *Geothermics* 57, 258–274;
Schnaidt (2015) PhD thesis, University of Adelaide.

**`medfilt3D(M, kernel_size, boundary_mode, maxiter, Out)`**
Apply an iterated median filter to an n-D array `M`.
`kernel_size` sets the filter window per axis.
Returns filtered array `G`.

**`anidiff3D(M, ckappa, dgamma, foption, maxiter, Out)`**
Convenience wrapper for `anisodiff3D` with named diffusion parameters.
Returns filtered array `G`.

**`anisodiff3D(stack, niter, kappa, gamma, step, option, ploton)`**
3-D Perona–Malik anisotropic diffusion.
`option=1`: gradient-magnitude gating; `option=2`: Lorentzian gating.
Returns `stackout`.

Reference: Perona & Malik (1990) *IEEE TPAMI* 12(7), 629–639.
Original MATLAB: P. Kovesi (UWA); Python port: A. Muldal (Oxford).

**`gauss3D(Kshape, Ksigma)`**
Generate a normalised 3-D Gaussian kernel of shape `Kshape`
and standard deviation `Ksigma`.
Returns kernel array `K`.

**`get_good_sites(q_val, q_thresh, out)`**
Return indices of sites where quality metric `q_val < q_thresh`.
Returns `good` (index tuple from `numpy.where`).

---

## Known Issues / Notes

- `crossgrad`: `m1.dim` corrected to `m1.ndim` (Apr 2026).
- `crossgrad`: `numpy.zeros_like(g1, cgdim)` (invalid signature)
  corrected to `numpy.zeros((numpy.size(g1), cgdim))` (Apr 2026).
- `fractrans`: `m == None` corrected to `m is None` (Apr 2026).

---

## Change Log

| Date | Author | Note |
|---|---|---|
| Apr 2021 | vrath | Initial version (AEMpyX). |
| Jul 2023 | vrath | `crossgrad` added. |
| Apr 2026 | vrath | Removed unused imports; bug fixes; provenance; readme. |
