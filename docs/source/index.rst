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


From basic to advanced helper functions
---------------------------------------
ScopeSim Templates is a python package, and is therefore by nature infinitely extendable.

As it is impossible for us to know all the details about your specific science case, we provide a small selection of basic objects (star cluster, elliptical galaxy, etc).
Feel free to start with these to get started with ScopeSim.

However if your needs outgrow the basic objects, we encourage you to extended the objects to fit your specific science case.
In this case **we strongly encourage you to get in contact with us adding your code in the form of a subpackage**.
You can do this either by opening an issue on Github, or by emailing one of the developers.


Available subpackages
---------------------

Documentation for all the helper funcitons contained in each package can be found in the API documentation for each package.

* Basic
* Advanced
* MICADO calibration


The Basic subpackage
--------------------

* in ``basic.stars``:
   * :func:`stars`
   * :func:`stars_field`
   * :func:`star_grid`
* in ``basic.galaxy``
   * :func:`spiral_two_component`
* in ``basic.misc``
* :func:`empty_sky`


