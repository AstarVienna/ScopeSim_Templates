Basic astronomical objects
==========================
The ``basic`` subpackage contains a series of basic astronomical objects that can be used for prelimiinry simulations with all ScopeSim instrument packages.

Currently the subpackage contains the following functions:

* `.basic`
   * ``.stars``:
      * :func:`star`
      * :func:`stars`
      * :func:`cluster`
      * :func:`stars_field`
      * :func:`star_grid`
   * ``.galaxy``
      * :func:`spiral_two_component`
      * :func:`elliptical`
   * ``.misc``
      * :func:`empty_sky`
      

Examples
--------

A young star cluster with a core radius r_c=1pc, M=1000 Msun and located in the LMC (d=50kpc):

.. jupyter-execute::
    :raises:
    
    from scopesim_templates.basic.stars import cluster
    src = cluster(mass=1E3, distance=50000, core_radius=1):
    
Lets have a look inside the object:

.. jupyter-execute::
    :raises:
    
    src.fields[0]

Here we can see the spatial information is in the form of an astropy ``Table``.    
    
The columns ``x`` and ``y`` show the position of each star in ``arcsec`` relative to the centre of the field of view.  
    
The column ``ref`` connects each star in this table to a spectrum in the following list:
    
.. jupyter-execute::
    :raises:
    
    src.spectra
    
When ``ScopeSim`` ingests this ``Source`` object, it will look primarily at these three columns.

Now for a graphical representation of the cluster as it will be seen by ``ScopeSim``:

.. jupyter-execute::
    :raises:
    
    from matplotlib import pyplot
    %matplotlib inline
    
    plt.scatter(src.fields[0]["x"], src.fields[0]["y"])
    
    
    
    
    


    

