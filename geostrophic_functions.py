"""
Rosalie Cormier, 2023, based on code by Andrew Delman
"""

import xarray as xr
import numpy as np
import ecco_v4_py as ecco
import xgcm

Omega = (2 * np.pi) / 86164 #Earth angular velocity

def to_radians(angle):
    
    """
    Converts degrees to radians.
    """
    
    return angle * np.pi / 180

def comp_geos_vel(ecco_ds_grid, pressure, dens):
    
    """
    Computes derivatives of pressure.
    
    ecco_ds_grid = grid dataset
    pressure = pressure data to differentiate
    dens = density data
    ds_vel = DataSet containing velocities
    """
    
    xgcm_grid = ecco.get_llc_grid(ecco_ds_grid)
    
    #Compute derivatives of pressure in x and y
    
    d_press_dx = (xgcm_grid.diff(pressure, axis="X", boundary='extend')) / ecco_ds_grid.dxC
    d_press_dy = (xgcm_grid.diff(pressure, axis="Y", boundary='extend')) / ecco_ds_grid.dyC
    
    #Convert DataArray content from dask to np arrays
    
    d_press_dx.data = d_press_dx.values
    d_press_dy.data = d_press_dy.values
    
    #Interpolate to centres of grid cells
    press_grads_interp = xgcm_grid.interp_2d_vector({'X': d_press_dx, 'Y': d_press_dy}, boundary='extend')
    
    dp_dx, dp_dy = press_grads_interp['X'], press_grads_interp['Y']
    dp_dx.name = 'dp_dx'
    dp_dy.name = 'dp_dy'
    
    #Compute RHS of geostrophic-balance equations
    
    GB_RHS_1, GB_RHS_2 = dp_dx / dens, - dp_dy / dens
    GB_RHS_1 = GB_RHS_1.where(ecco_ds_grid.maskC) #Mask land areas
    GB_RHS_2 = GB_RHS_2.where(ecco_ds_grid.maskC) #Mask land areas
    
    #Compute Coriolis param. from latitudes of grid cell centres
    
    lat = to_radians(ecco_ds_grid.YC)
    f = 2 * Omega * np.sin(lat)
    
    #Compute geostrophic velocity components
    
    u_g, v_g = GB_RHS_2 / f, GB_RHS_1 / f
    u_g = u_g.where(ecco_ds_grid.maskC) #Mask land areas
    v_g = v_g.where(ecco_ds_grid.maskC) #Mask land areas
    
    return u_g, v_g


def mask_delta_u(mask_threshold, u_complex):
    
    u, v = np.real(u_complex), np.imag(u_complex)
    speed = np.sqrt(u**2 + v**2)
    
    mask_small_speed = (speed < mask_threshold)

    return mask_small_speed

def comp_delta_u_norm(ecco_ds_grid, k_val, u_complex, u_g_complex, mask=None):
    
    """
    Computes Delta-u diagnostic for geostrophic balance.
    
    ecco_ds_grid = grid DataSet
    k_val = depth index of interest
    """

    u, v = np.real(u_complex), np.imag(u_complex)
    u_g, v_g = np.real(u_g_complex), np.imag(u_g_complex)
    
    u_diff = u - u_g
    v_diff = v - v_g
    vel_diff_complex = u_diff + (1j * v_diff)
    vel_complex = u + (1j * v)
    vel_diff_abs = np.abs(vel_diff_complex)
    vel_abs = np.abs(vel_complex)
    vel_diff_norm = vel_diff_abs / vel_abs

    if mask is not None:
        vel_diff_norm = np.where(mask, np.nan, vel_diff_norm)
    
    return vel_diff_norm