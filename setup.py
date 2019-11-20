#!/usr/bin/env python3
"""
ScopeSim: A python package to simulate telescope observations
"""

from datetime import datetime
from distutils.core import setup
from setuptools import find_packages

# Version number
MAJOR = 0
MINOR = 1
ATTR = 'dev0'

VERSION = '%d.%d%s' % (MAJOR, MINOR, ATTR)


def setup_package():
    setup(name='ScopeSim',
          version=VERSION,
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
                            ],
          classifiers=["Programming Language :: Python :: 3",
                       "License :: OSI Approved :: MIT License",
                       "Operating System :: OS Independent",
                       "Intended Audience :: Science/Research",
                       "Topic :: Scientific/Engineering :: Astronomy", ]
          )


if __name__ == '__main__':
    setup_package()
