.. _scopesim_templates-api:

**************************
ScopeSim_Templates Package
**************************


``scopesim.stellar`` module
===========================

Templates to simulate stellar sources

.. autosummary::
   :nosignatures:
   :toctree: generated/

   scopesim_templates.stellar.star
   scopesim_templates.stellar.stars
   scopesim_templates.stellar.star_field
   scopesim_templates.stellar.star_grid
   scopesim_templates.stellar.cluster

``scopesim_templates.extragalactic`` module
===========================================

Templates to simulate extragalactic sources

.. autosummary::
   :nosignatures:
   :toctree: generated/

   scopesim_templates.extragalactic.galaxy
   scopesim_templates.extragalactic.galaxy3d
   scopesim_templates.extragalactic.spiral_two_component
   scopesim_templates.extragalactic.elliptical

``scopesim_templates.misc`` module
==================================

Templates that could be used to simulate more general sources

.. autosummary::
   :nosignatures:
   :toctree: generated/

   scopesim_templates.misc.point_source
   scopesim_templates.misc.uniform_source
   scopesim_templates.misc.source_from_file
   scopesim_templates.misc.source_from_array
   scopesim_templates.misc.source_from_imagehdu
   scopesim_templates.misc.source_from_imagehdu_with_flux
   scopesim_templates.misc.source_from_cube


``scopesim_templates.calibration`` module
==========================================

Simple templates that could be used to simulate calibration frames. Make sure to turn off the corresponding
effects during the simulation

.. autosummary::
   :nosignatures:
   :toctree: generated/

   scopesim_templates.calibration.empty_sky
   scopesim_templates.calibration.flat_field
   scopesim_templates.calibration.lamp

