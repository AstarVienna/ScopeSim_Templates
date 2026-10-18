# Units and conventions

Most of the bugs fixed in version 0.10 came from different templates making
different silent assumptions. This page collects the conventions that all
templates now follow.

## Positions

- Positions are **offsets in arcsec from the field centre**. This applies to the
  `x`/`y` columns of table fields, and to the `CRVAL1`/`CRVAL2` of image fields
  (in degrees, as FITS requires). ScopeSim treats image coordinates as offsets
  from the optical axis. An image whose `CRVAL` is an absolute RA/Dec, such as
  53° for the HUDF, therefore lands far outside the field of view.
- `x` increases to the right and `y` upwards on the detector. Templates built
  from real sky positions (`galactic_centre`, `s301`,
  `hubble_ultra_deep_field`, ...) map **East to +x** and North to +y, so the
  scene appears mirrored with respect to the usual "East left" sky
  orientation.
- The `ra`/`dec` arguments that some templates accept describe the field
  centre. They are stored in `src.meta` and are *not* used to place the
  source. For image templates (`source_from_array`, `galaxy`, ...) leave
  them at 0, because they end up in `CRVAL`.
- Pixel scales are in arcsec per pixel. Plain floats are taken as arcsec and
  `astropy` Quantities are converted.
- Image arrays are indexed `data[y, x]` (NumPy order). Functions that take an
  image size as a tuple say whether they mean `(nx, ny)` or a NumPy shape.

## Angles

- Position angles of the `galaxy`, `galaxy3d` and `elliptical` profiles are in
  **degrees, counter-clockwise from the +x axis**. In `galaxy3d` this angle
  applies to the light, velocity and dispersion maps alike.
- `hubble_ultra_deep_field` converts the catalogue's astronomical position
  angle (from North through East) to this convention itself.

## Brightness

A bare number passed as a brightness is interpreted as follows. Pass an
`astropy` Quantity (`u.mag`, `u.ABmag`, `u.STmag`, `u.Jy`, ...) to be
explicit.

| template | unitless amplitude means |
|---|---|
| `star`, `stars`, `star_field`, `star_grid` | Vega magnitude (`u.mag`) |
| `point_source`, `uniform_source`, `source_from_array`, `source_from_file` | AB magnitude |
| `galaxy`, `galaxy3d`, `elliptical` | AB magnitude |
| `hubble_ultra_deep_field` | AB magnitudes from the catalogue (F160W) |

Image templates normalise the image to sum to 1, so the amplitude is the
brightness of the whole image, unless the docstring says it is a surface
brightness (e.g. `uniform_source`, which is per arcsec²).

## Randomness

- Every stochastic template takes `seed`. It may be an `int` (reproducible), a
  `numpy.random.Generator` (to share one random stream between several calls)
  or `None`.
- With `seed=None` (the default) every call gives a **new realisation**. To
  make every template reproducible without passing `seed` each time, set a
  package-wide default:

  ```python
  from scopesim_templates.rc import __config__
  __config__["!random.seed"] = 42
  ```

- Templates never reseed or consume NumPy's *global* random state
  (`np.random.seed`, `np.random.random`, ...), which belongs to your code.
  The test suite enforces this.

## Time

- Templates with orbital motion (`galactic_centre`,
  `galactic_centre_background`, `s301`) take an epoch as an
  `astropy.time.Time`, a `datetime`, an ISO or Julian-year string (e.g.
  `"2031-10-20"`, `"J2031.8"`), or a float MJD.
- `globular_cluster` takes `time` in years since its (seed-determined)
  initial orbital phase. `omega_cen`/`evolve_omega_cen` evolve a realisation
  to a given epoch; see their docstrings.
- To animate one realisation, keep the `seed` fixed and change only the
  time.

## Metadata

Every template records how it was called in `src.meta["function_call"]`,
and its name in `src.meta["object"]`. A longer human-readable label, where
there is one, is in `src.meta["description"]`.
