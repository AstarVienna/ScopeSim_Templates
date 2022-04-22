.. _stellar:

Stellar Module Examples
=======================

A young star cluster with a core radius r_c=1pc, M=1000 Msun and located in the LMC (d=50kpc):

.. code-block:: python
    
    from scopesim_templates.stellar import cluster
    src = cluster(mass=1E3, distance=50000, core_radius=1):
    
Lets have a look inside the object:

.. code-block:: python
    
    src.fields[0]

Here we can see the spatial information is in the form of an astropy ``Table``.    
    
The columns ``x`` and ``y`` show the position of each star in ``arcsec`` relative to the centre of the field of view.  
    
The column ``ref`` connects each star in this table to a spectrum in the following list:
    
.. code-block:: python

    src.spectra
    
When ``ScopeSim`` ingests this ``Source`` object, it will look primarily at these three columns.

Now for a graphical representation of the cluster as it will be seen by ``ScopeSim``:

.. code-block:: python
    
    from matplotlib import pyplot
    %matplotlib inline
    
    plt.scatter(src.fields[0]["x"], src.fields[0]["y"])
    

