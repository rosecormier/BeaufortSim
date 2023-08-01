"""
Vorticity-related functions.

Rosalie Cormier, 2023
"""

#from os.path import join

import ecco_v4_py as ecco

from functions_ecco_general import ecco_resample#, load_dataset
from functions_geostrophy import comp_f

##############################

def comp_vorticity(grid_llc, u_mean, v_mean, dx, dy, cell_area):
    return (-grid_llc.diff(u_mean * dx, 'Y') + grid_llc.diff(v_mean * dy, 'X')).squeeze() / cell_area

##############################

def comp_normal_strain(grid_llc, u_mean, v_mean, dx, dy, cell_area):
    return (grid_llc.diff(u_mean * dy, 'X') - grid_llc.diff(v_mean * dx, 'Y')) / cell_area

##############################

def comp_shear_strain(grid_llc, u_mean, v_mean, dx, dy, cell_area):
    return (grid_llc.diff(v_mean * dy, 'X') + grid_llc.diff(u_mean * dx, 'Y')) / cell_area

##############################

def comp_total_strain(ds_grid, normal_strain, shear_strain, latmin, latmax, lonmin, lonmax, resolution):

    normal_strain_field = ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    shear_strain_field = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    
    return normal_strain_field + shear_strain_field

##############################

def comp_OkuboWeiss(omega, normal_strain, shear_strain):
    return normal_strain**2 + shear_strain**2 - omega**2

##############################

def comp_local_Ro(omega, y):
    
    f = comp_f(y)
    
    return abs(omega) / f

##############################

def get_OW_field(ds_grid, ds_vel, k, lats_lons, resolution, zeta_field):
    
    """
    Computes OW in a form that can be visualized.
    """
    
    latmin, latmax, lonmin, lonmax = lats_lons #Set spatial bounds

    xgcm_grid = ecco.get_llc_grid(ds_grid)
    
    normal_strain = comp_normal_strain(xgcm_grid, ds_vel['UVEL'], ds_vel['VVEL'], ds_grid.dxG, ds_grid.dyG, ds_grid.rA).isel(k=k).squeeze()
    shear_strain = comp_shear_strain(xgcm_grid, ds_vel['UVEL'], ds_vel['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz).isel(k=k).squeeze()
    
    normal_strain= ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    shear_strain = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    
    OW = comp_OkuboWeiss(zeta_field, normal_strain, shear_strain) #Compute OW 
    
    return OW