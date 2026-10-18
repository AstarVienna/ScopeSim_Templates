# ScopeSim Templates

[![Poetry](https://img.shields.io/endpoint?url=https://python-poetry.org/badge/v0.json)](https://python-poetry.org/)
[![Build Status](https://github.com/AstarVienna/ScopeSim_Templates/actions/workflows/tests.yml/badge.svg)](https://github.com/AstarVienna/ScopeSim_Templates/actions/workflows/tests.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/scopesim-templates/badge/?version=latest)](https://scopesim-templates.readthedocs.io)
[![Codecov](https://img.shields.io/codecov/c/github/AstarVienna/ScopeSim_Templates/branch/main?logo=codecov)](https://app.codecov.io/gh/AstarVienna/ScopeSim_Templates/tree/main)

[![PyPI - Version](https://img.shields.io/pypi/v/ScopeSim-Templates)](https://pypi.org/project/ScopeSim-Templates/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/ScopeSim-Templates)
![GitHub Release Date](https://img.shields.io/github/release-date/AstarVienna/ScopeSim_Templates)
![dev version](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2FAstarVienna%2FScopeSim_Templates%2Fmain%2Fpyproject.toml&query=%24.project.version&label=dev%20version&color=teal)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This package contains a number of templates to generate `Source` objects to be used in simulations with the [ScopeSim Simulator](https://scopesim.readthedocs.io/en/latest/)

The Documentation can be found here: https://scopesim-templates.readthedocs.io/en/latest

## Installation

The best way to install the software is to use `pip`

``` bash
pip install scopesim_templates
```

To install the development version you can clone the repository

``` bash
git clone https://github.com/AstarVienna/ScopeSim_Templates/
cd ScopeSim_Templates
pip install -e .
```

## Requirements

- numpy
- astropy
- [synphot](https://synphot.readthedocs.io/en/latest/index.html)
- [ScopeSim](https://scopesim.readthedocs.io)
- [speXtra](https://spextra.readthedocs.io)
- [pyckles](https://pyckles.readthedocs.io)

## The `Source` object

The above functions are created to ease the creation of standard sources but the power of the `Source` object doesn't end there and can be used to create endless source possibilities.

In a nutshell a `Source` object contains a spatial and spectral description of the sources. The spectral description is in the form of `synphot` spectra and the spatial description can be an `astropy` table referencing the spectra or a fits image. `Source` can also accept data cubes. The `speXtra` package contains an extensive library of spectral templates that can be used with the sources. Please check the relevant documentation.

## `Source` templates included

All templates return a `scopesim.Source`. Positions are offsets in arcsec from
the field centre. Stochastic templates take a `seed` argument (an `int` or a
`numpy.random.Generator`); without it every call gives a new realisation (see
[Units and conventions](https://scopesim-templates.readthedocs.io/en/latest/conventions.html)).

- `stellar`
  - `star`, `stars`: one or more stars with given spectral types and magnitudes
  - `star_field`: a field of stars at random positions
  - `star_grid`: a field of stars on a regular grid
  - `cluster`: a young (age = 0) cluster of a given mass
  - `globular_cluster`: a globular cluster around a central IMBH, with time evolution
  - `omega_cen`, `evolve_omega_cen`: an ω Centauri model (King profile, Jeans
    velocities, IMBH, RGB stars) and its Keplerian time evolution
  - `galactic_centre`: the Gillessen+2017 S-stars around Sgr A* at any epoch,
    optionally with a background population (`background_mass=`)
  - `galactic_centre_background`: a hypothetical unseen main-sequence population
    around Sgr A* (up to the GRAVITY 2022 limit of ~3000 Msun)
  - `s301`, `s301_orbit`: the S301 star (Mang et al. 2026) at any epoch, for either
    orbit solution

- `extragalactic`
  - `galaxy`: a Sérsic model with a user-selected SED from the `speXtra` database
  - `galaxy3d`: like `galaxy`, plus a velocity field and velocity dispersion field
  - `spiral_two_component`: a two-component spiral with a young disk SED and an old bulge SED
  - `elliptical`: a Sérsic model using the Brown et al. (2014) SEDs
  - `hubble_ultra_deep_field`: the HUDF from published catalogues as Sérsic
    profiles, optionally with star-forming clumps

- `misc`
  - `point_source`: like `star`, but using any SED from the `speXtra` database
  - `uniform_source`: a uniform source with any SED from `speXtra`
  - `source_from_array`: a source from a 2D `numpy` array
  - `source_from_file`: a source from an image in a FITS file
  - `source_from_imagehdu`: a source from an `ImageHDU` with an arbitrary flux and scale
  - `source_from_imagehdu_with_flux`: a source from an `ImageHDU` whose flux per pixel is known
  - `source_from_cube`: a source from a 3D datacube
  - `nebulosity`: a procedural, semi-realistic nebula image (use with `source_from_array`)

- `calibration`
  - `lamp`: a calibration lamp, i.e. a homogeneous source with emission lines
  - `flat_field`: a flat field
  - `empty_sky`: an empty sky with no sources

- `micado`, `metis`: instrument-specific calibration sources and test fields
  (flat lamps, pinhole masks, line lists, lasers, VIKING fields)

Please see the documentation for how to use each source, and contact us (raise
an issue or submit a pull request) if more specialized sources are needed.
