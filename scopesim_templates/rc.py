from os import path as pth
import sys
import yaml

if pth.exists("F:/Work/ScopeSim/"):
    sys.path.append("F:/Work/ScopeSim/")

import scopesim

Source = scopesim.source.source.Source
SystemDict = scopesim.system_dict.SystemDict
ter_curve_utils = scopesim.effects.ter_curves_utils
im_plane_utils = scopesim.optics.image_plane_utils
scopesim_utils = scopesim.utils

PKG_DIR = pth.abspath(pth.dirname(__file__))
with open(pth.join(PKG_DIR, "defaults.yaml")) as fname:
    __config__ = SystemDict(yaml.full_load(fname))
