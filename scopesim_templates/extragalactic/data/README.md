# Vendored data

## `hudf_catalogue.fits`

Galaxy catalogue for the Hubble Ultra Deep Field, used by
`scopesim_templates.extragalactic.hubble_ultra_deep_field`.

* 3152 galaxies within a 2.5 arcmin cone around the HUDF centre
  (RA, Dec = 53.1625, -27.79142 deg; Beckwith et al. 2006, AJ 132, 1729).
* Built by `scripts/make_hudf_catalogue.py` (VizieR cone searches plus a
  0.5 arcsec positional cross-match; median separation 0.03 arcsec).

| column | meaning | source |
|---|---|---|
| `id` | van der Wel+2012 GOODS-S identifier | vdW12 |
| `ra`, `dec` | J2000, deg | vdW12 |
| `z` | spectroscopic redshift where available, else EAZY peak photo-z | S14 |
| `z_spec` | True if `z` is spectroscopic | S14 |
| `m_F160W` | GALFIT total F160W magnitude, AB | vdW12 |
| `r_eff_arcsec` | half-light semi-major axis, arcsec | vdW12 |
| `sersic_n` | Sersic index | vdW12 |
| `axis_ratio` | b/a | vdW12 |
| `position_angle_deg` | major-axis PA, deg, 0 = North, 90 = East | vdW12 |
| `log_mstar` | log10 stellar mass / Msun (FAST) | S14 |
| `m_F606W` | total F606W magnitude, AB (NaN if undetected) | S14 |
| `match_sep_arcsec` | cross-match separation | - |

Only F160W fits with GALFIT quality flag `Q = 0` are kept.

### Attribution

If you use this catalogue, please cite the original works:

* van der Wel, A., Bell, E. F., Häussler, B., et al. 2012, ApJS, 203, 24,
  "Structural Parameters of Galaxies in CANDELS"
  (VizieR `J/ApJS/203/24`).
* Skelton, R. E., Whitaker, K. E., Momcheva, I. G., et al. 2014, ApJS, 214,
  24, "3D-HST WFC3-selected Photometric Catalogs in the Five CANDELS/3D-HST
  Fields" (VizieR `J/ApJS/214/24`).

This research made use of the VizieR catalogue access tool, CDS, Strasbourg,
France.
