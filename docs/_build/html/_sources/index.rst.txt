.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
      th { display:none; }
    </style>


.. image:: _static/logos/logo_long_scopesim_templates_t.png
    :width: 600 px
    :alt: Welcome to the ScopeSim_Templates Documentation!
    :align: center


|logo| Another tool from the `A* Vienna software team <https://astarvienna.github.io/>`_

.. |logo| image:: https://raw.githubusercontent.com/AstarVienna/astarvienna.github.io/main/logos/star_small_t.png
   :height: 30px
   :align: middle


A library of templates and helper functions for creating
:class:`scopesim.source.source.Source` objects that can be used to run `ScopeSim` simulations.

In short :class:`scopesim.source.source.Source` objects contain a description of the spatial and
spectral information of the source. For more information see :ref:`here <Source Object>`.



Installation
------------

This package has been released on PyPi::

   pip install scopesim_templates


From basic to advanced helper functions
---------------------------------------
ScopeSim Templates is a python package, and is therefore by nature infinitely extendable.

As it is impossible for us to know all the details about your specific science case, we provide a
small selection of basic objects (star cluster, elliptical galaxy, etc).
Feel free to start with these to get started with ScopeSim.

However if your needs outgrow the basic objects, we encourage you to extended the objects to fit your
specific science case. In this case **we strongly encourage you to get in contact with us adding your code
in the form of a subpackage**. You can do this either by opening an issue on Github, or by emailing one of the developers.


Available subpackages
---------------------

Documentation for all the helper functions contained in each package can be found in the API documentation for each package.

* ``stellar``:
   * :func:`star`: Places a single star on the field
   * :func:`stars`: Places a list of stars on the field
   * :func:`cluster`:  Creates an age=0 cluster with a user selected mass
   * :func:`stars_field`: Creates field of stars with random positions
   * :func:`star_grid`: Creates a field of stars with regular positions

* ``extragalactic``
   * :func:`galaxy`: A simple sersic model with a user selected SED from the ``speXtra`` database
   * :func:`galaxy3D`:  A more complex model that includes a velocity field and velocity dispersion field
   * :func:`spiral_two_component`:  Simple two component model with an outer spiral young SED and an old SED bulge
   * :func:`elliptical`: Another sersic model using the Brown SEDs

* ``misc``
   * :func:`point_source`: similar to :func:`star` but using any SED from the ``speXtra`` database
   * :func:`uniform_source`: Creates a uniform source with any SED from ``speXtra``
   * :func:`source_from_image_hdu`: creates a source from an ``ImageHDU`` with an arbitrary flux and scale
   * :func:`source_from_imagehdu_with_flux`: creates a source from an ``ImageHDU`` where the flux/pixel is known
   * :func:`source_from_file`: Load the source from a fits file. Depending on the characteristics other
               functions may be more suitable
   * :func:`source_from_array`: General function to create a source from a 2D ``numpy`` array
   * :func:`source_from_cube`:  Wrapper to create a source from a 3D datacube

* ``calibration``:
   * :func:`empty_sky`: To simulate a sky without no other sources
   * :func:`lamp`: Simulates a calibration lamp, i.e. a homogenous source with emissions lines
   * :func:`flat_field`: Simulates a flat field



Contact
-------
If you find an issue with ScopeSim Templates, please let us know via the
`Github issues page <https://github.com/AstarVienna/ScopeSim_Templates/issues>`_


.. toctree::
   :maxdepth: 2

   notebooks/starting.ipynb
   notebooks/stellar.ipynb
   notebooks/extragalactic.ipynb
   source_object
   modules



More...
=======

.. toctree::
   :maxdepth: 1
   :titlesonly:
