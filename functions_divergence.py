"""
Functions to compute velocity divergence.

Rosalie Cormier, 2023
"""

import numpy as np
import ecco_v4_py as ecco

##############################

def comp_2d_divergence(grid_llc, u, v, dx, dy, cell_area):
    return (grid_llc.diff(u * dy, 'X') + grid_llc.diff(v * dx, 'Y')).squeeze() / cell_area

##############################

def comp_2d_Ek_divergence(grid_llc, u_Ek, v_Ek, dx, dy, cell_area, ds_grid, lats_lons, resolution):

    latmin, latmax, lonmin, lonmax = lats_lons
    
    dx.data = dx.values
    dx = grid_llc.interp(dx, 'X')
    
    dy.data = dy.values
    dy = grid_llc.interp(dy, 'Y')
    
    du_dx = grid_llc.diff(u_Ek * dy, 'X', boundary='extend').squeeze()
    du_dx = grid_llc.interp(du_dx, 'X')
    
    dv_dy = grid_llc.diff(v_Ek * dx, 'Y', boundary='extend').squeeze()
    dv_dy = grid_llc.interp(dv_dy, 'Y')
    
    div_u_Ek = (du_dx + dv_dy) / cell_area
   
    #Convert to useful field
    lon_centers, lat_centers, lon_edges, lat_edges, field = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, div_u_Ek, latmin, latmax, resolution, lonmin, lonmax, resolution, fill_value=np.NaN, mapping_method='nearest_neighbor', radius_of_influence=120000)

    return lon_centers, lat_centers, lon_edges, lat_edges, field 