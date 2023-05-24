"""
Rosalie Cormier, 2023, based on code by Andrew Delman
"""

import numpy as np
import ecco_v4_py as ecco
import xgcm

Omega = (2 * np.pi) / 86164 #Earth angular velocity

def to_radians(angle):
    
    """
    Converts to radians.
    
    angle = angle in degrees
    """
    
    return angle * np.pi / 180

def comp_geos_vel(ecco_ds_grid, pressure, dens, ds_vel):
    
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
    d_press_dy = (xgcm_grid.diff(pressure, axis="Y", boundary='extent')) / ecco_ds_grid.dyC
    
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
    
    #Interpolate velocities to centres of grid cells
    
    ds_vel.UVEL.data, ds_vel.VVEL.data = ds_vel.UVEL.values, ds_vel.VVEL.values
    vel_interp = xgcm_grid.interp_2d_vector({'X': ds_vel.UVEL, 'Y': ds_vel.VVEL}, boundary='extend')
    u, v = vel_interp['X'], vel_interp['Y']
    
    #Compute Coriolis param. from latitudes of grid cell centres
    
    lat = to_radians(ecco_ds_grid.YC)
    f = 2 * Omega * np.sin(lat)
    
    #Compute geostrophic velocity components
    
    u_g, v_g = GB_RHS_2 / f, GB_RHS_1 / f
    u_g = u_g.where(ecco_ds_grid.maskC) #Mask land areas
    v_g = v_g.where(ecco_ds_grid.maskC) #Mask land areas
    
    return u_g, v_g

def comp_delta_u_norm(u, v, u_g, v_g):
    
    """
    Computes Delta-u diagnostic for geostrophic balance.
    
    u, v = components of velocity
    u_g, v_g = components of geostrophic velocity
    """
    
    num_x, num_y = u - u_g, v - v_g
    norm_num = np.sqrt(num_x**2 + num_y**2)
    denom_x, denom_y = u, v
    norm_denom = np.sqrt(denom_x**2 + denom_y**2)
    
    return norm_num / norm_denom