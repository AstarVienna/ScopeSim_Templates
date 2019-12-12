.. ScopeSim_templates documentation master file, created by
   sphinx-quickstart on Mon Nov 11 12:47:33 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the ScopeSim Templates documentation!
================================================

This package provides templates and helper functions for creating
:class:`scopesim.source.source.Source` objects.

This package has been released on PyPi::

   pip install scopesim_templates

Available templates
-------------------

Currently the package only includes the following functions:

* in ``basic.stars``:
   * :func:`stars`
   * :func:`stars_field`
   * :func:`star_grid`
* in ``basic.galaxy``
   * :func:`spiral_two_component`
* in ``basic.misc``
* :func:`empty_sky`


We will be porting the old templates from ``SimCADO`` in the very near future


Input required for a Source object
----------------------------------

* A spatial description
* A spectral description

Spatial description
-------------------

Series of arrays
++++++++++++++++
Point sources can be described by (x, y) positions [arcsec] from the
centre of the field of view.

* ``x, y``
* ``weight``
* ``ref``


Astropy Table
+++++++++++++


FITS ImageHDU
+++++++++++++
Extended souces should provide an intensity map in the form of an
``fits.HDUImage`` object. The header keywords should contain information
regarding position or the centre of the image (``CRPIXn``) relative to the
centre of the field of view (``CRVALn``) [degrees], and pixel scale
[``CDELTn`` | ``CDn_m``] [degress / pixel]



Spectral description
--------------------

Series of arrays
++++++++++++++++
* ``wave``
* ``flux``


``synphot.SourceSpectrum``
++++++++++++++++++++++++++


