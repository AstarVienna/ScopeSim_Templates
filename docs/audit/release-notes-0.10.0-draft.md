# ScopeSim_Templates 0.10.0 — draft release notes

*Draft for the "notes" field of the **Prepare release** workflow (rule: `minor`, 0.9.1.dev1 → 0.10.0).
Branch `kl/globular_cluster_imbh`. Not yet merged to `main` or released.*

## ⚠ Behaviour changes

- **No default random seed.** Stochastic templates (`star_field`, `cluster`, `globular_cluster`, `nebulosity`,
  VIKING fields, ...) now give a new realisation on every call. Previously most of them silently used
  `seed=9001`. For reproducible output pass `seed=<int>` (or a `numpy.random.Generator`), or set
  `scopesim_templates.rc.__config__["!random.seed"]`.
- **The global NumPy RNG is no longer touched.** Templates used to call `np.random.seed(...)`, which reset your
  own random stream. They now use private generators.
- **Wrong output was fixed, so some results change:**
  - `source_from_array`, `source_from_file`, `galaxy` and `uniform_source`: CRPIX1/2 were swapped for non-square images.
  - MICADO `line_list` and METIS `laser_spectrum`: the image was transposed relative to its WCS.
  - MICADO `flatlamp` was off-centre by ½ pixel.
  - `stars(library=...)` with non-pickles libraries gave stars of one type the first star's brightness.
  - VIKING fields: every galaxy had the first galaxy's flux.
  - `galaxy`/`galaxy3d`: `theta` is now in degrees for the light profile (it was radians, while the
    velocity and dispersion maps used degrees). `galaxy3d` no longer ignores `z` and `amplitude` for
    `Spextrum` input.
  - `hubble_ultra_deep_field`: galaxies now actually land on the detector (positions were absolute RA/Dec),
    and position angles are converted correctly.
- `spiral_two_component(fluxes=...)` raises `NotImplementedError` instead of silently ignoring the argument;
  `offset=` is now applied.
- `galaxy`/`galaxy3d` accept `str`, `Spextrum` and `synphot.SourceSpectrum` SEDs alike.

## New templates

- `stellar.galactic_centre_background`: a hypothetical main-sequence population on Keplerian orbits around
  Sgr A*, up to the GRAVITY 2022 mass limit (`total_mass="gravity2022"`, ~3000 Msun).
  `galactic_centre(..., background_mass=...)` adds it to the S-star field.
- `stellar.s301`, `stellar.s301_orbit`: the S301 star (Mang et al. 2026), either orbit solution, at any epoch.
- `stellar.omega_cen`, `stellar.evolve_omega_cen`: an ω Centauri model with a central IMBH and Keplerian time
  evolution of one realisation.
- `stellar.globular_cluster`: a globular cluster around a central IMBH, with time evolution.
- `extragalactic.hubble_ultra_deep_field`: reworked. It uses a vendored 3152-galaxy catalogue (van der Wel+2012
  × Skelton+2014), flux-conserving Sérsic profiles, optional star-forming clumps, and field selection helpers.
  A 60″ field builds in about 1 s.
- `misc.nebulosity`: procedural, semi-realistic nebula images. This replaces the NGC 6888 FITS image.

## Fixes

- `micado.pinhole_mask()` (always raised `TypeError`) and `from scopesim_templates.micado import *` (raised
  `TypeError`) work again.
- `load_stars_source(ra=...)`, `load_galaxies_source(ns=...)` and `star_field(seed=...)` no longer crash.
- `utils.samplers` works on NumPy ≥ 2.4 (`np.trapz` was removed).
- `stellar.imf`: the multiplicity path works again (`np.object`), and caller-supplied mass limits are no longer
  mutated.
- Templates no longer modify their inputs (`stars` spec-type lists, FITS headers in `source_from_imagehdu*`,
  wavelength arrays in `pinhole_mask`). `source_from_file` closes its file.
- Recording the function call in `meta` can no longer crash source creation.
- The Kepler solver used by `galactic_centre`/`globular_cluster` is robust up to e → 1.

## Deprecations (nothing removed)

- `micado.empty_sky` → use `micado.darkness` (kept for MicadoWISE).
- `source_from_image` now emits a `DeprecationWarning` → use `source_from_array`.
- VIKING `random_seed=` → `seed=`.

## Packaging and docs

- The wheel no longer ships tests or EPS figures (1.74 MB, previously 3.45 MB); the sdist still includes the tests.
  `matplotlib` is no longer a runtime dependency. `scopesim_templates.__version__` is available.
- New **Units and conventions** documentation page; README lists all templates.
- Tests: global-RNG guard on every test; regression tests for each fixed bug; an end-to-end smoke test that
  observes every public template with ScopeSim; slow tests are marked `slow`.
