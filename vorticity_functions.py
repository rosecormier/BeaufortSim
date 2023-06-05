"""
Rosalie Cormier, 2023

Vorticity-related functions
"""

import numpy as np
import dask.array as daskarray
import xgcm
import xarray as xr
import ecco_v4_py as ecco

from ecco_general import ecco_resample

def comp_vorticity(grid_llc, mean_u, mean_v, dxC, dyC, cell_area):
    return (-grid_llc.diff(mean_u * dxC, 'Y') + grid_llc.diff(mean_v * dyC, 'X')).squeeze() / cell_area

def comp_normal_strain(grid_llc, mean_u, mean_v, dxG, dyG, cell_area):
    return (grid_llc.diff(mean_u * dyG, 'X') + - grid_llc.diff(mean_v * dxG, 'Y')) / cell_area

def comp_shear_strain(grid_llc, mean_u, mean_v, dxC, dyC, cell_area):
    return (grid_llc.diff(mean_v * dyC, 'X') + grid_llc.diff(mean_u * dxC, 'Y')) / cell_area

def comp_total_strain(grid_llc, ds_grid, mean_u, mean_v, dxC, dyC, dxG, dyG, rA, rAz, k_val, latmin, latmax, lonmin, lonmax, resolution):
    
    normal = comp_normal_strain(grid_llc, mean_u, mean_v, dxG, dyG, rA)
    normal_strain_field = ecco_resample(ds_grid, normal.isel(k=k_val).squeeze(), latmin, latmax, lonmin, lonmax, resolution)[4] 
    
    shear = comp_shear_strain(grid_llc, mean_u, mean_v, dxC, dyC, rAz)
    shear_strain_field = ecco_resample(ds_grid, shear.isel(k=k_val).squeeze(), latmin, latmax, lonmin, lonmax, resolution)[4]
    
    return normal_strain_field + shear_strain_field

def comp_OkuboWeiss(grid_llc, ds_grid, mean_u, mean_v, dxC, dyC, dxG, dyG, rA, rAz, k_val, latmin, latmax, lonmin, lonmax, resolution):
    
    omega = comp_vorticity(grid_llc, mean_u, mean_v, dxC, dyC, rAz)
    omega_field = ecco_resample(ds_grid, omega.isel(k=k_val).squeeze(), latmin, latmax, lonmin, lonmax, resolution)[4]
    
    normal_strain = comp_normal_strain(grid_llc, mean_u, mean_v, dxG, dyG, rA)
    normal_strain_field = ecco_resample(ds_grid, normal_strain.isel(k=k_val).squeeze(), latmin, latmax, lonmin, lonmax, resolution)[4]

    shear_strain = comp_shear_strain(grid_llc, mean_u, mean_v, dxC, dyC, rAz)
    shear_strain_field = ecco_resample(ds_grid, shear_strain.isel(k=k_val).squeeze(), latmin, latmax, lonmin, lonmax, resolution)[4]
    
    W = normal_strain_field**2 + shear_strain_field**2 - omega_field**2
    
    return W