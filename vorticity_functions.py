"""
Rosalie Cormier, 2023

Vorticity-related functions
"""

import numpy as np
import dask.array as daskarray
#from matplotlib import pyplot as plt
import xgcm
import xarray as xr
import ecco_v4_py as ecco

def comp_vorticity(grid_llc, mean_u, mean_v, dxC, dyC, cell_area, k_val):
    
    zeta = (-grid_llc.diff(mean_u * dxC, 'Y') + grid_llc.diff(mean_v * dyC, 'X')) / cell_area
    zeta = zeta.squeeze()
    
    return zeta

def comp_normal_strain(grid_llc, mean_u, mean_v, dxG, dyG, hFacW, hFacS, drF, cell_area, k_val):
    
    u_transport = mean_u * dyG * hFacW * drF
    v_transport = mean_v * dxG * hFacS * drF
    
    du_dx = grid_llc.diff(u_transport, 'X') / cell_area
    dv_dy = grid_llc.diff(v_transport, 'Y') / cell_area
    
    normal_strain = du_dx - dv_dy
    
    return normal_strain

def comp_shear_strain(grid_llc, mean_u, mean_v, dxC, dyC, cell_area, k_val):
    
    shear_strain = (grid_llc.diff(mean_v * dyC, 'X') + grid_llc.diff(mean_u * dxC, 'Y')) / cell_area
    
    return shear_strain

def comp_total_strain(grid_llc, ds_grid, mean_u, mean_v, dxC, dyC, dxG, dyG, hFacW, hFacS, drF, rA, rAz, k_val, latmin, latmax, lonmin, lonmax, resolution):
    
    normal = comp_normal_strain(grid_llc, mean_u, mean_v, dxG, dyG, hFacW, hFacS, drF, rA, k_val)
    normal_strain_field = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, normal.isel(k=k_val).squeeze(), latmin, latmax, resolution, lonmin, lonmax, resolution, fill_value=np.NaN, mapping_method='nearest_neighbor', radius_of_influence=120000)[4]
    
    shear = comp_shear_strain(grid_llc, mean_u, mean_v, dxC, dyC, rAz, k_val)
    shear_strain_field = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, shear.isel(k=k_val).squeeze(), latmin, latmax, resolution, lonmin, lonmax, resolution, fill_value=np.NaN, mapping_method='nearest_neighbor', radius_of_influence=120000)[4]
    
    return normal_strain_field + shear_strain_field

def comp_OkuboWeiss(grid_llc, ds_grid, mean_u, mean_v, dxC, dyC, dxG, dyG, hFacW, hFacS, drF, rA, rAz, k_val, latmin, latmax, lonmin, lonmax, resolution):
    
    omega = comp_vorticity(grid_llc, mean_u, mean_v, dxC, dyC, rAz, k_val)
    omega_field = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, omega.isel(k=k_val).squeeze(), latmin, latmax, resolution, lonmin, lonmax, resolution, fill_value=np.NaN, mapping_method='nearest_neighbor', radius_of_influence=120000)[4]
    
    normal_strain = comp_normal_strain(grid_llc, mean_u, mean_v, dxG, dyG, hFacW, hFacS, drF, rA, k_val)
    normal_strain_field = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, normal_strain.isel(k=k_val).squeeze(), latmin, latmax, resolution, lonmin, lonmax, resolution, fill_value=np.NaN, mapping_method='nearest_neighbor', radius_of_influence=120000)[4]

    shear_strain = comp_shear_strain(grid_llc, mean_u, mean_v, dxC, dyC, rAz, k_val)
    shear_strain_field = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, shear_strain.isel(k=k_val).squeeze(), latmin, latmax, resolution, lonmin, lonmax, resolution, fill_value=np.NaN, mapping_method='nearest_neighbor', radius_of_influence=120000)[4]
    
    W = normal_strain_field**2 + shear_strain_field**2 - omega_field**2
    
    return W