#!/usr/bin/env python3
"""
ScopeSim: A python package to simulate telescope observations
"""

from datetime import datetime
from distutils.core import setup
from setuptools import find_packages

# Version number
with open('scopesim_templates/version.py') as f:
    __version__ = f.readline().split("'")[1]


def setup_package():
    setup(name='ScopeSim_Templates',
          version=__version__,
          description="On-sky source templates for ScopeSim",
          author="Kieran Leschinski",
          author_email="kieran.leschinski@unive.ac.at",
          url="https://github.com/astronomyk/ScopeSim_Templates",
          package_dir={'scopesim_templates': 'scopesim_templates'},
          packages=find_packages(),
          include_package_data=True,
          install_requires=["numpy>=1.13",
                            "scipy>0.17",
                            "astropy>1.1.2",
                            "matplotlib>1.5.0",
                            "requests>2.0",
                            "pyyaml>3",
                            "synphot>0.1",
                            "pyckles",
                            "scopesim",
                            ],
          classifiers=["Programming Language :: Python :: 3",
                       "License :: OSI Approved :: MIT License",
                       "Operating System :: OS Independent",
                       "Intended Audience :: Science/Research",
                       "Topic :: Scientific/Engineering :: Astronomy", ]
          )


if __name__ == '__main__':
    setup_package()
