# TODO 


## Reorganization and funtions missing (marked in bf)

Comments: 
    - ``agn_template()`` ``quasar_template()`` (and a ``SNe_template``) can be used with a generalized point_source source
    - nebula sources may only be wrappers around ``source_from_image()``



- stellar
    - star()
    - stars()
    - star_grid()
    - star_field()

- stellar_clusters
    - basic_cluster()
    - **art_pop_wrapper()**
    - **astrometric_cluster()**
    - **galactic_center()**

- galaxies
    - elliptical()
    - galaxy()
    - galaxy3d()
    - **vela_wrapper()**
    - **agn_template()**
    - **quasar_template()**
    
- galaxy_clusters
    - **realistic_background_galaxy_distribution()**
    - **mock_HUDF()**
    - **galaxy_cluster()**  (NFW + LF + mass-size distro)

- nebula
    - **orion(distance, magnitude)**
    - **PNe nebula**

- solar_system_objects
    - 

- misc
    - from_image
    - **from_cube**   (poor_man_source_from_cube does exist)
    - from_point_source_table
    - **from_moments**

- backgrounds
    - empty_sky()
    - **background_stars(ra, dec, catalogue="simbad")**
    - **background_galaxies(ra, dec, catalogue="simbad")**  (catalogs may not have the resolution, maybe a galaxy correlation function will be enough or connect to simulations)

- calibration
    - dark_frame()
    - **pinhole_mask()**
    

