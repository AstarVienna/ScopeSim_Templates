%matplotlib inline
import os
import datetime
 
import shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

import scopesim as sim
import scopesim_templates as sim_tp
from scopesim_templates.misc.misc import source_from_array

def donut_ellipse(a1, b1, ecc, inc, width, height):
    
    array = np.zeros((height, width), dtype=np.float32)
    
    x0 = width // 2
    y0 = height // 2
    
    # Outer ellipse
    inc_rad = np.radians(inc)
    cos_inc = np.cos(inc_rad)
    sin_inc = np.sin(inc_rad)
    
    for y in range(height):
        for x in range(width):
            norm_x = (x - x0) / a1
            norm_y = (y - y0) / b1
            transformed_x = norm_x * cos_inc - norm_y * sin_inc
            transformed_y = norm_x * sin_inc + norm_y * cos_inc
            if transformed_x**2 / (1 -ecc**2) + transformed_y**2 <= 1:
                array[y, x] = 1
                
    # Inner ellipse
    ratio = b1 / a1
    a2 = a1 * ratio
    b2 = b1 * ratio
    
    for y in range(height):
        for x in range(width):
            norm_x = (x - x0) / a2
            norm_y = (y - y0) / b2
            transformed_x = norm_x * cos_inc - norm_y * sin_inc
            transformed_y = norm_x * sin_inc + norm_y * cos_inc
            if transformed_x**2 / (1 -ecc**2) + transformed_y**2 <= 1:
                array[y, x] = 0
    
    src = source_from_array(arr=array, sed="bosz/lr/a0v",
                            pixel_scale=0.0057*u.arcsec, amplitude=10*u.ABmag,
                            filter_curve="L")
    return src