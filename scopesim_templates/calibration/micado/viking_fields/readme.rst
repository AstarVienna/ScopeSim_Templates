Task from Ric
=============

Requirements
------------
Sources
- 3 catalogues of stars:
    - illumination correction
    - standard stars, and
    - science exposures
    - all 1x1 arcmin fields
- 3 catalogues of galaxies:
    - all 1x1 arcmin fields
    - interchangable with star fields

- stars catalogues
    - RA, Dec, P_*, Hmag
    - Coordinates, Probability of a star from VIKING, Brightness
- bg galaxies catalogues:
    - H_mag Reff(mas) type (1 Disk, 0 Elliptical)
    - Sersic profiles, n=2 (disk), n=4 (elliptical)
    - Random position
    - Random ellipticity and orientation for Disks

Observations
- MCAO, H-band
- Dithered exposures for each star catalogue (dx, dy) = (7, 3) arcsec
- DIT set to avoid saturation (probs ~15sec)
- NDIT=1, illum, std_star
- NDIT>1, sci


Emails
------

Can you please create some simulations for the pipeline to work with
(and maybe send me the python script that does it, so I can then copy
and adapt it for future use).
These are for the 3 sets of on-sky fields for illumination correction,
standard stars, and science exposures.
It seemed best to base them on real data, so i have used regions of the
VIKING survey, and the attached file hopefully illustrates what to
expect for each. I have produced plots of the part of the VIKING field I
used (left side of figures) and a zoom=in on an arcmin field selected on
the right. In both cases, bright (H<15mag) stars are marked in blue,
other stars in green, and galaxies in red. The associated cat_xxx.dat
files show the centre of each field and the sources in a square arcmin
centered on that.

What I'd like for each is a set of dithered exposures as described in
Yixian's google doc, which I think she has shared with you.
The dithers should be a few arcsec so we can create a sky background
from the on-source frames; and the DIT should be set to avoid
saturation. For illum and stdstar probably we can live with NDIT=1; but
for the science field we may need NDIT>1 to go deep enough.

We can simply adopt the MCAO standard AO performance. And for the
galaxies, please use arbitrary parameters, and put R_eff at 0.3-1.5
arcsec or so (actually I have no idea how big the galaxies are, but this
seems reasonable for our purpose).

I think that's it. The idea is to have then quasi-realistic without
putting too much effort into it. (And the pipeline team will have fun
identifying what's in the fields from a rather larger catalogue, in
order to do the flux calibration)



Trilegal is easy, since most of the parameters I leave as default.
If you want to have a look, it is at http://stev.oapd.inaf.it/cgi-bin/trilegal_1.6
I only touch the first 2 boxes, giving the coordinates of the pointing, the size of the field to simulate, the filter set to use (I don't worry about translating between different sorts of JHK, I just assume they are near enough for what we want), and the limiting magnitude.



So attached are some interim catalogues and plots where I've added a realistic distribution of stars from Trilegal - the small purple crosses in the right panels of the plots.
This is really an education because there are actually rather fewer stars than you'd expect. That's really all there is (and in the south galactic pole 'science' field, there is only 1 extra star down to 32 mag). Of course, there will be many more galaxies, but this will take a bit more time. next week.



So I've looked up galaxy population statistics based on Cresci+06, Maihara+01, Yan+98, an d Totani & Yoshii 2000.
I hesitate to give you code to calculate a galaxy population because it is rather too crude and hand-wavey a representation for me to feel comfortable sharing it, but I can give you the distributions of magnitude, size, and type I have generated with it. Which are ok for us during pipeline development & testing, but not for a wider audience.
Each file is a random selection for a 1armcin field, and can go with any of the catalogue files i sent (the interim ones, with the extra few stars) since i assume the galaxy distribution is homogeneous.
The files have 3 columns for H-band magnitude (also Vega system), the effective radius in mas, and the type (1=disk, 0=spheroidal).
You can then create random inclinations (axis ratios) for the disky ones yourself.
And I also let you distribute them randomly over a 1arcmin field

So you can see there are a lot of faint galaxies, but only when you go to really faint magnitudes.



There's always an after-thought, which I could have explained at the time. This one is about the realism of the faint galaxy populations I made.
In terms of number counts, density, magnitudes, that is all realistic. The limitation is in the details of the galaxy shapes. We use simple sersic profiles (with n=2 or n=4 typically). But real galaxies at high z are quite different from this, and will have a lot of prominent structure overlaying the general profile shape. And in many cases they might be completely different. This is too complex to simulate in a simple way.
So if we were to want to go further in the direction of realism - which I don't think is needed for testing the pipeline - we would have to make use of high resolution cosmological simulations like IllustrisTNG (TNG50) or NewHorizon. Even TNG50 has a best resolution of ~75pc and 10^5Msun, which means the smallest strctures they can see are down to a few hundred parsecs. NewHorizon is a bit better at ~35pc resolution. But I still expect that creating an equivalent of an H or K-band image would be a significant amount of work for somebody. Would be interesting to get a better idea of what sort of detail & structure we might actually see, and nice for showing in presentations, but I don't really see much other real use.
