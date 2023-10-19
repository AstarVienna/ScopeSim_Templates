from importlib import metadata
version = metadata.version(__package__)
date = '2023-10-19 08:00:00 GMT'
yaml_descriptions = """
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