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

def comp_vorticity(grid_llc, u_mean, v_mean, dx, dy, cell_area):
    return (-grid_llc.diff(u_mean * dx, 'Y') + grid_llc.diff(v_mean * dy, 'X')).squeeze() / cell_area

def comp_normal_strain(grid_llc, u_mean, v_mean, dx, dy, cell_area):
    return (grid_llc.diff(u_mean * dy, 'X') - grid_llc.diff(v_mean * dx, 'Y')) / cell_area

def comp_shear_strain(grid_llc, u_mean, v_mean, dx, dy, cell_area):
    return (grid_llc.diff(v_mean * dy, 'X') + grid_llc.diff(u_mean * dx, 'Y')) / cell_area

def comp_total_strain(ds_grid, normal_strain, shear_strain, latmin, latmax, lonmin, lonmax, resolution):

    normal_strain_field = ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    shear_strain_field = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    
    return normal_strain_field + shear_strain_field

def comp_OkuboWeiss(omega, normal_strain, shear_strain):
    return normal_strain**2 + shear_strain**2 - omega**2