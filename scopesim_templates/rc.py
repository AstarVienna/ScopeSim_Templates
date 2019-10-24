from os import path as pth
import yaml

try:
    import scopesim
except:
    import sys
    sys.path.append("C:/Work/ScopeSim/")
    import scopesim

Source = scopesim.source.source.Source
SystemDict = scopesim.system_dict.SystemDict

PKG_DIR = pth.abspath(pth.dirname(__file__))
with open(pth.join(PKG_DIR, "defaults.yaml")) as fname:
    __config__ = SystemDict(yaml.load(fname))
