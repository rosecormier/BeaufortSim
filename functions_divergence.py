"""
Functions to compute velocity divergence.

Rosalie Cormier, 2023
"""

from functions_ecco_general import ecco_resample

##############################

def comp_2d_divergence(grid_llc, u, v, dx, dy, cell_area, ds_grid, lats_lons, resolution):
    
    latmin, latmax, lonmin, lonmax = lats_lons

    dx.data = dx.values
    dx = grid_llc.interp(dx, 'X')
    
    dy.data = dy.values
    dy = grid_llc.interp(dy, 'Y')
    
    du_dx = grid_llc.diff(u * dy, 'X', boundary='extend').squeeze()
    du_dx = grid_llc.interp(du_dx, 'X')
    
    dv_dy = grid_llc.diff(v * dx, 'Y', boundary='extend').squeeze()
    dv_dy = grid_llc.interp(dv_dy, 'Y')
   
    return (du_dx + dv_dy) / cell_area