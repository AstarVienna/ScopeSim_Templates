---
substitutions:
  logo: |-
    ```{image} https://raw.githubusercontent.com/AstarVienna/astarvienna.github.io/main/logos/star_small_t.png
    :align: middle
    :height: 30px
    ```
---

```{raw} html
<style media="screen" type="text/css">
  h1 { display:none; }
  th { display:none; }
</style>
```

```{image} _static/logos/logo_long_scopesim_templates_t.png
:align: center
:alt: Welcome to the ScopeSim_Templates Documentation!
:width: 600 px
```

# ScopeSim Templates

Another tool from the [A\* Vienna software team](https://astarvienna.github.io/)

A library of templates and helper functions for creating
{class}`scopesim.source.source.Source` objects that can be used to run `ScopeSim` simulations.

In short {class}`scopesim.source.source.Source` objects contain a description of the spatial and
spectral information of the source. For more information see [here](./source_object.md).

## Installation

This package has been released on PyPi:

```bash
pip install scopesim_templates
```

## From basic to advanced helper functions

ScopeSim Templates is a python package, and is therefore by nature infinitely extendable.

As it is impossible for us to know all the details about your specific science case, we provide a small selection of basic objects (star cluster, elliptical galaxy, etc).
A comprehensive list of those can be found below in the API reference (grouped into subpackages).
Feel free to start with these to get started with ScopeSim.

However if your needs outgrow the basic objects, we encourage you to extended the objects to fit your specific science case. In this case **we strongly encourage you to get in contact with us adding your code in the form of a subpackage**.
You can do this either by opening an issue on Github, or by emailing one of the developers.

## Contact

If you find an issue with ScopeSim Templates, please let us know via the
[Github issues page](https://github.com/AstarVienna/ScopeSim_Templates/issues)

## Contents

```{toctree}
:maxdepth: 2

notebooks/starting.md
notebooks/stellar.md
notebooks/extragalactic.md
source_object
```

## API reference

```{eval-rst}
.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:
   :caption: Package Contents

   scopesim_templates.stellar
   scopesim_templates.extragalactic
   scopesim_templates.calibration.calibration
   scopesim_templates.misc.misc
```
