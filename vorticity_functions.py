"""
Rosalie Cormier, 2023

Vorticity-related functions
"""

import numpy as np
#import xarray
import dask.array as daskarray
from matplotlib import pyplot as plt
#%matplotlib inline
#from matplotlib.colors import SymLogNorm
#from xmitgcm import open_mdsdataset
import xgcm
#from mpl_toolkits.basemap import Basemap, cm, shiftgrid
import xarray as xr

def comp_vorticity(grid_llc, mean_u, mean_v, dxC, dyC, cell_area, k_val):
    
    zeta = (-grid_llc.diff(mean_u * dxC, 'Y') + grid_llc.diff(mean_v * dyC, 'X')) / cell_area
    zeta = zeta.isel(k=k_val).squeeze()
    
    return zeta