# ScopeSim Templates

[![Build Status](https://github.com/AstarVienna/ScopeSim_Templates/actions/workflows/tests.yml/badge.svg)](https://github.com/AstarVienna/ScopeSim_Templates/actions/workflows/tests.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/scopesim-templates/badge/?version=latest)](https://scopesim-templates.readthedocs.io/en/latest)
[![Codecov](https://img.shields.io/codecov/c/github/AstarVienna/ScopeSim_Templates/branch/main?logo=codecov)](https://app.codecov.io/gh/AstarVienna/ScopeSim_Templates/tree/main)
[![PyPI - Version](https://img.shields.io/pypi/v/ScopeSim-Templates)](https://pypi.org/project/ScopeSim-Templates/)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This packages contain a number of templates to generate `Source` objects to be used in simulations with the [ScopeSim Simulator](https://scopesim.readthedocs.io/en/latest/)

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
- [ScopeSim](https://scopesim.readthedocs.io/en/latest/)
- [speXtra](https://spextra.readthedocs.io/en/latest/)
- [pyckles](https://pyckles.readthedocs.io/en/latest/)

## The `Source` object

The above functions are created to easy the creation of standard sources but the power of the `Source` object doesn't end there and can be used to create endless sources possibilities.

In a nutshell a `Source` object contains a spacial and spectral description of the sources. The spectral description are in the form of `synphot` spectra and the spacial description can be an `astropy` table referencing the spectra or a fits image. `Source` can also accept datacubes. The `speXtra` package contains an extensive library of spectral templates that can be used with the sources. Please check the relevant documentation.

## `Source` templates included

Currently, the package covers the most typical sources used in astronomy:

- `stellar`
  - `star`: Places a single star on the field
  - `stars`: Places a list of stars on the field
  - `star_field`: Creates field of stars with random positions
  - `star_grid`: Creates a field of stars with regular positions
  - `cluster`: Creates an age=0 cluster with a user selected mass

- `extragalactic`  
  - `galaxy`: A simple sersic model with a user selected SED from the `speXtra` database
  - `galaxy3D`: A more complex model that includes a velocity field and velocity dispersion field
  - `spiral_two_component`: Simple two component model with an outer spiral young SED and an old SED bulge
  - `elliptical`: Another sersic model using the Brown SEDs

- `misc`  
  - `point_source`: similar to `star` but using any SED from the `speXtra` database
  - `uniform_source`: Creates a uniform source with any SED from `speXtra`
  - `source_from_file`: Load the source from a fits file. Depending on the characteristics othe functions may be more suitable
  - `source_from_imagehdu`: creates a source from an `ImageHDU` with an arbitrary flux and scale
  - `source_from_imagehdu_with_flux`: creates a source from an `ImageHDU` where the flux/pixel is known
  - `source_from_array`: General function to create a source from a 2D `numpy` array
  - `source_from_cube`: Wrapper to create a source from a 3D datacube

- `calibration`  
  - `lamp`: Simulates a calibration lamp, i.e. a homogenous source with emissions lines
  - `flat_field`: Simulates a flat field
  - `empty_sky`: To simulate a sky without no other sources

Please see the documentation how to use each particular source and contact us (raise an issue or submit a pull request) if more specialized sources are needed.
