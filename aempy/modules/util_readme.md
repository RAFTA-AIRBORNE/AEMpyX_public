# util_readme.md

## Module: `util.py`

**Purpose:** General-purpose utility functions for AEM data workflows.

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
| Standard library | `os`, `sys`, `ast`, `inspect`, `warnings`, `fnmatch`, `random` |
| Third-party | `numpy`, `pyproj` (`CRS`, `Transformer`), `shapely` |

---

## Naming Conventions

- Input parameters use CapitalizedCamelCase where standardised;
  legacy lowercase signatures preserved where not yet updated.
- Boolean output-control parameters follow the pattern `out` / `Out*`.

---

## Function Reference

### Path and environment

**`to_ospath(inpath, opsys)`**
Convert path separators to the native OS convention.
Returns `outpath`.

**`check_env(envar, action)`**
Verify that a conda / shell environment variable is set.
Exits on failure when `action='error'`.

### File utilities

**`sample_list(in_list, method, out)`**
Sub-sample a list. `method[0]`: `'pass'` | `'rand'` | `'step'`.
Returns `out_list`.

**`get_data_list(how, sort, fullpath, out)`**
Build a list of data files by `'search'` (wildcard) or `'read'` (list file).
Returns list of filenames.

**`get_filebase(file)`**
Return `(name, ext)` for a given file path.

**`get_filelist(searchstr, searchpath, sortedlist, fullpath)`**
Generate a filename list from a directory using `fnmatch` wildcards.
Returns list of filenames.

**`get_files(SearchString, SearchDirectory)`**
Thin wrapper around `fnmatch.filter`. Returns `FileList`.

**`splitall(path)`**
Split a path into all its component parts. Returns `allparts`.

**`change_filename(old_filename, how)`**
Modify a filename by appending, prepending, or replacing a substring.
`how[0]`: `'app'` | `'prep'` | `'repl'`. Returns `new_filename`.

**`strcount(keyword, fname)`**
Count occurrences of `keyword` in a text file. Returns count.

**`strdelete(keyword, fname_in, fname_out, out)`**
Delete all lines containing `keyword` from a file.

**`strreplace(key_in, key_out, fname_in, fname_out)`**
Replace `key_in` with `key_out` in every line of a file.

### Array utilities

**`check_finite(arr_in, out)`**
Remove rows containing any non-finite value from a 2-D array.
Returns `arr_out`.

**`unique(list, out)`**
Return unique elements of a list preserving order.

**`nearly_equal(a, b, sig_fig)`**
Test near-equality to `sig_fig` significant figures. Returns bool.

**`stack_ragged(array_list, axis)`**
Stack a list of ragged arrays into a single array with an index.
Returns `(stacked, idx)`.

**`save_stacked_arrays(fname, array_list, axis, compressed)`**
Save stacked ragged arrays to NPZ.

**`load_stacked_arrays(fname, axis)`**
Load and split stacked ragged arrays from NPZ.

**`add_object_npz(filein, xkeys, xobjects, fileout)`**
Add new key/object pairs to an existing NPZ file.

**`del_object_npz(filein, xkeys, fileout)`**
Delete keys from an NPZ file.

### Geometry

**`get_nearest_point(point, line, tol)`**
Find the index of the nearest point on a profile. Returns `nearest`.

**`get_direction_angle(p1, p2)`**
Direction angle (degrees, from North) and distance from `p1` to `p2`.
Returns `(ang, length)`.

**`project_to_line(x, y, line)`**
Project a point onto a line defined by two points.
Returns `(xn, yn)`.

**`segment_distance(p, p1, p2, axis, return_t, segment)`**
N-dimensional point-to-segment distance.
Returns `dist` — or `(dist, t)` when `return_t=True`.

**`find_nearest(site0, sitevec)`**
Brute-force nearest-neighbour in 2-D. Returns `minpos`.

**`gen_searchgrid(Points, XLimits, dX, YLimits, dY, Out)`**
Generate a 2-D search grid mapping data points to cells.
Returns object array `p[ix, iy]` of point index lists.

**`gen_grid_latlon(LatLimits, nLat, LonLimits, nLon, out)`**
Generate equidistant 1-D grids in latitude/longitude.
Returns `(Lat, Lon)`.

**`gen_grid_utm(XLimits, nX, YLimits, nY, out)`**
Generate equidistant 1-D grids in UTM metres.
Returns `(X, Y)`.

### Coordinate projection (pyproj)

| Function | Converts |
|---|---|
| `project_wgs_to_geoid` | WGS84 ellipsoid height → geoid height |
| `project_utm_to_geoid` | UTM + ellipsoid height → geoid height |
| `project_gk_to_latlon` | Gauss–Krüger → WGS84 lat/lon |
| `get_utm_zone` | lat/lon → UTM EPSG code |
| `project_latlon_to_utm` | WGS84 lat/lon → UTM |
| `project_utm_to_latlon` | UTM → WGS84 lat/lon |
| `project_latlon_to_itm` | WGS84 lat/lon → Irish Transverse Mercator |
| `project_itm_to_latlon` | ITM → WGS84 lat/lon |
| `project_itm_to_utm` | ITM → UTM |
| `project_utm_to_itm` | UTM → ITM |
| `project_utm_to_utm` | UTM zone → UTM zone |

### Polygon / spatial selection

**`modify_polygon(Polygons, Operator, Params, Out)`**
Apply `'rotation'`, `'intersection'`, or `'union'` to shapely polygons.

**`extract_data_poly(Data, PolyPoints, method, Out)`**
Extract rows of `Data` whose (x, y) coordinates fall inside a polygon.
`method`: `'env'` (envelope) | `'con'` (convex hull) | `'shp'` (given shape).

**`point_inside_polygon(x, y, poly, method)`**
Test whether point (x, y) lies inside `poly`.
`method`: `'shapely'` (default) or pure-Python ray-casting.

**`extract_data_rect(Data, Corners, Out)`**
Extract rows of `Data` within an axis-aligned rectangle
`Corners = [Xmin, Ymin, Xmax, Ymax]`.

### Miscellaneous

**`fractrans(m, x, a)`**
Fractional derivative of order `a` of `m` over domain `x` (`differint`).

**`anisodiff(img, niter, kappa, gamma, step, option)`**
2-D Perona–Malik anisotropic diffusion.

**`anisodiff3(stack, niter, kappa, gamma, step, option)`**
3-D Perona–Malik anisotropic diffusion.

**`ms2ohmm(ms)`**
Convert conductivity in mS/m to resistivity in Ohm·m.

**`print_title(version, fname, form, out)`**
Print and return a formatted title string with version and file
modification date.

**`list_functions(filename)` / `list_functions_in(this_module)`**
List function names defined in a file or module.

**`make_pdf_catalog(WorkDir, PdfList, FileName)`**
Concatenate PDF site plots into a single catalogue file (requires `fitz`).

**`export_to_vtk(points, scale, modmesh, methmesh, exportfile)`**
Write a 3-D model point cloud to VTK format (requires `evtk`).

---

## Known Issues / Notes

- `get_filebase`: return statement is inside an `else` branch that is
  unreachable; use `os.path.splitext(os.path.basename(file))` directly.
- Deprecated `pyproj.transform()` in `project_gk_to_latlon` replaced
  with `Transformer.from_crs`. (Apr 2026)
- Duplicate definitions of `project_wgs_to_geoid` and
  `project_utm_to_geoid` removed. (Apr 2026)
- `list_functions_in`: `import this_module` replaced with correct
  `inspect.getmembers` pattern. (Apr 2026)
- `(p1 == None)` comparisons replaced with `is None` throughout. (Apr 2026)

---

## Change Log

| Date | Author | Note |
|---|---|---|
| Nov 2020 | vrath | Initial version (AEMpyX). |
| Jan 2022 | vrath | `anisodiff`: simplified for Python 3. |
| Jan 2023 | vrath | `get_data_list`, `sample_list`: extended modes. |
| Apr 2023 | vrath | `get_filelist`: updated. |
| Apr 2026 | vrath | Bug fixes; removed duplicates; provenance; readme. |
