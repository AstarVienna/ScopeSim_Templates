# Galactic Centre and omega Cen science cases

Three time-evolvable templates ported from the MAVIS 2026 workshop
feasibility studies. All of them take a `seed` (an `int`, a
`numpy.random.Generator`, or `None`; `None` gives a new realisation on every
call) and never touch the global NumPy random state. Positions are arcsec
offsets from the target centre (`x` towards increasing RA, `y` north); the
target's RA/Dec is in `src.meta`.

Epochs follow the convention of `galactic_centre()`: an `astropy.time.Time`,
a `datetime`, an ISO string (`"2031-10-22"`), a Julian-epoch string
(`"J2031.81"`), or a bare number, which is read as an **MJD**.

## Unseen stars around Sgr A*

`galactic_centre_background()` samples a Kroupa IMF up to an old-population
turn-off and puts every star on a random Keplerian orbit about Sgr A*.
There is only an upper limit on this mass: GRAVITY Collaboration 2022
(A&A 657, L12) allow up to ~3000 Msun inside S2's apocentre.

```python
from scopesim_templates.stellar import galactic_centre, galactic_centre_background

bg = galactic_centre_background("2031-10-22", total_mass="gravity2022",
                                mag_limit=24, seed=42)

# S-stars plus the background, as two fields of one Source
src = galactic_centre("2031-10-22", background_mass=3000, seed=42,
                      background_kwargs={"mag_limit": 24})
```

Re-using the same integer `seed` at different epochs follows the same
background stars along their orbits.

## S301

`s301()` places S301 (Mang et al. 2026) on either of the two published,
mirror-degenerate orbit solutions and shifts its spectrum to the apparent
radial velocity (Doppler, transverse Doppler and gravitational redshift).
Sgr A*'s NIR counterpart is added as a power law by default.

```python
from astropy.time import Time
from scopesim_templates.stellar import s301, s301_orbit

orbit = s301_orbit(1)
t_peri = orbit.next_pericentre("2026-01-01")          # ~2031.81
src = s301(Time(t_peri + 1 / 365.25, format="jyear"), solution=1)
orbit.apparent_radial_velocity(t_peri + 1 / 365.25)   # km/s
```

## omega Centauri

`omega_cen()` is a literature-calibrated toy model of omega Cen (analytic
King profile, Jeans velocities including the IMBH, Kroupa IMF truncated at
the turn-off plus an RGB normalised to the observed M/L). It complements the
simpler, generic `globular_cluster()`.

```python
from scopesim_templates.stellar import omega_cen, evolve_omega_cen

src = omega_cen(fov=10, imbh_mass=1e5, mag_limit=27, seed=1)  # epoch J2030
later = evolve_omega_cen(src, "J2045")                       # same stars

# or directly
later = omega_cen(fov=10, imbh_mass=1e5, mag_limit=27, seed=1, time="J2045")
```

At a fixed `seed`, changing `imbh_mass` changes only the stellar velocities,
so different black-hole masses can be compared star for star.
