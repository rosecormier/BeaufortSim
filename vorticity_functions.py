"""
Rosalie Cormier, 2023

Vorticity-related functions
"""

import numpy as np
import dask.array as daskarray
from matplotlib import pyplot as plt
import xgcm
import xarray as xr

def comp_vorticity(grid_llc, mean_u, mean_v, dxC, dyC, cell_area, k_val):
    
    zeta = (-grid_llc.diff(mean_u * dxC, 'Y') + grid_llc.diff(mean_v * dyC, 'X')) / cell_area
    #zeta = zeta.isel(k=k_val).squeeze()
    zeta = zeta.squeeze()
    return zeta

def comp_OkuboWeiss(grid_llc, mean_u, mean_v, dxC, dyC, cell_area, k_val):
    
    omega = comp_vorticity(grid_llc, mean_u, mean_v, dxC, dyC, cell_area, k_val)
    
    normal_strain = (grid_llc.diff(mean_u * dyC, 'X') - grid_llc.diff(mean_v * dxC, 'Y')) / cell_area
    normal_strain = normal_strain.isel(k=k_val).squeeze()
    
    shear_strain = (grid_llc.diff(mean_v * dyC, 'X') + grid_llc.diff(mean_u * dxC, 'Y')) / cell_area
    shear_strain = shear_strain.isel(k=k_val).squeeze()
    
    W = normal_strain**2 + shear_strain**2 - omega**2
    
    return W