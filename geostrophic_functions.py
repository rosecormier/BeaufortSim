"""
Rosalie Cormier, 2023
"""

import ecco_v4_py as ecco
import xgcm

def pressure_derivs(ecco_ds_grid, pressure):
    
    """
    Computes derivatives of pressure.
    
    ecco_ds_grid = grid dataset
    pressure = pressure data to differentiate
    """
    
    xgcm_grid = ecco.get_llc_grid(ecco_ds_grid)
    
    #Compute derivatives of pressure in x and y
    
    d_press_dx = (xgcm_grid.diff(pressure, axis="X", boundary='extend')) / ecco_ds_grid.dxC
    d_press_dy = (xgcm_grid.diff(pressure, axis="Y", boundary='extent')) / ecco_ds_grid.dyC
    
    #Convert DataArray content from dask to np arrays
    
    d_press_dx.data = d_press_dx.values
    d_press_dy.data = d_press_dy.values
    
    #Interpolate to centres of grid cells
    press_grads_interp = xgcm_grid.interp_2d_vector({"X": d_press_dx, "Y": d_press_dy}, boundary='extend')
    
    dp_dx, dp_dy = press_grads_interp['X'], press_grads_interp['Y']
    dp_dx.name = 'dp_dx'
    dp_dy.name = 'dp_dy'
    
    return dp_dx, dp_dy