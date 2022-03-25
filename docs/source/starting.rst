.. _start:

***************
Getting Started
***************


``scopesim.source.Source`` objects are composed of a spatial description and a spectral one. Spatial description
can be ``astropy.table.Table`` objects for point sources or ``astropy.fits.ImageHDU`` for extended sources.
Spectral description is provided as ``synphot.SourceSpectrum`` and compatible objects

For example, the creation of ``scopesim.source.Source`` objects might require quite a bit of interaction from the
user

.. jupyter-execute::


    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import synphot
    from scopesim import Source

    # creation of a
    x, y = np.meshgrid(np.arange(100), np.arange(100))
    img = np.exp(-1 * ( ( (x - 50) / 5)**2 + ( (y - 50) / 5)**2 ) )

    hdr = fits.Header(dict(NAXIS=2,
                           NAXIS1=img.shape[0]+1,
                           NAXIS2=img.shape[1]+1,
                           CRPIX1=0,
                           CRPIX2=0,
                           CRVAL1=0,
                           CRVAL2=0,
                           CDELT1=0.2/3600,
                           CDELT2=0.2/3600,
                           CUNIT1="DEG",
                           CUNIT2="DEG",
                           CTYPE1='RA---TAN',
                           CTYPE2='DEC--TAN'))
    hdu = fits.ImageHDU(data=img, header=hdr)

    wave = np.arange(1000, 35000, 10 )
    bb = synphot.models.BlackBody1D(temperature=5000)
    sp = synphot.SourceSpectrum(synphot.Empirical1D, points=wave, lookup_table=bb(wave))

    src = Source(image_hdu=hdu, spectra=sp)

    plt.imshow(src.fields[0].data)
    src.spectra[0].plot()