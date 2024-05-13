from importlib import metadata
version = metadata.version(__package__)
date = '2024-02-06 18:00:00 GMT'
yaml_descriptions = """
- version : 0.5.2
  date : 2024-05-13
  comment : Patch to add METIS calibration sources.
  changes :
  - Add scopesim_templates/metis/laser.py and pinhole_mask.py

- version : 0.5.1
  date : 2024-02-06
  comment : Patch to remove Python 3.8 and SystemDict.
  changes :
  - Drop support for Python 3.8 #77
  - Use astar_utils.NestedMapping instead of scopesim.system_dict.SystemDict #76
  - Fix scaling of spectrum #74
  - Add config file for auto-generated release notes #68
  - Comment out print statements, perhaps replace with logging later #72

- version : 0.5.0
  date : 2023-10-19
  comment : Fix off-by-one error.
  changes :
  - Fix off-by-one error.
  - Use pyproject.toml.
  - Various cleanups and fixes.

- version : 0.4.4
  date : 2023-03-08
  comment : Fix metadata bug in 2D galaxy.

- version : 0.4.3
  date : 2023-02-15
  comment : Add pinhole mask and minor bugfixes

- version : 0.4.2
  date : 2022-07-12
  comment : Updated DLC server to point to scopesim.univie.ac.at

- version : 0.4.1
  date : 2022-04-22
  comment : Updated package structure
"""