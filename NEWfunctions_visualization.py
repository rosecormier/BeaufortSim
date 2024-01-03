"""
Auxiliary functions for visualizing ECCO data.

Rosalie Cormier, 2024
"""

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py as ecco

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from os.path import join

from functions_comp_data_meta import load_comp_data_file
from functions_ecco_general import load_ECCO_data_file, load_grid, scalar_to_grid, vector_to_grid
from functions_field_variables import get_field_variable, get_vector_comps, field_is_primary

##############################

plt.rcParams['font.size'] = 16
plt.rcParams['text.usetex'] = True

##############################

def cbar_label(scalar_attr):
    
    """
    Returns label for plot colorbar.
    """

    cbar_label_dict = {'pressure': r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$', \
                       'density': r'Density anomaly $(kg/{m}^3)$', \
                       'vertical_vel': 'Velocity (m/s)', \
                       'vorticity': 'Vorticity (1/s)', \
                       'normal_strain': r'Normal strain $(1/s^2)$', \
                       'shear_strain': r'Shear strain $(1/s^2)$', \
                       '2D_div_vel': 'Horizontal velocity divergence (1/s)'}
    label = cbar_label_dict[scalar_attr]
    return label

##############################

#def pcolormesh_quiver_title(ecco_ds_grid, k_plot, datestr, scalar_attr, xvec_attr, resid=False):
    
#    """
#    Returns title for contourf-quiver plot.
#    """
"""    
    ds_grid = ecco_ds_grid.copy()
   
    depth = - ds_grid.Z[k_plot].values
    depthstr = str(depth) + ' m depth'
        
    scalar_dict = {'PHIHYDcR': 'Pressure anomaly', \
                  'OW': 'Okubo-Weiss parameter', \
                  'OW_geos': r'Okubo-Weiss parameter (computed from $\vec{u}_g$)', \
                  'Ro_l': 'Local Rossby number', \
                  'ZETA': 'Vorticity', \
                  'WVEL': 'Vertical component of water velocity', \
                  'DIVU': 'Divergence of horizontal water velocity', \
                  'DIVUEk': 'Divergence of Ekman current'}
    scalar_str = scalar_dict[scalar_attr]
    
    vector_dict = {'UVEL': 'horizontal water velocity', \
                  'UG': 'geostrophic water velocity', \
                  'UEk': 'Ekman current'}
    vector_str = vector_dict[xvec_attr]
    
    if resid:
        title = scalar_str + ' and ' + vector_str + \
            ' residuals (relative to annual mean) \n in Arctic Circle at {}, {} \n'.format(depthstr, \
                                                                                           datestr)

    else:
        title = scalar_str + ' and ' + vector_str + ' in BGR \n at {}, {} \n'.format(depthstr, \
                                                                                               datestr)
    
    return title
"""
##############################

def pcolormesh_k_title(ds_grid, k_plot, variable, datestr):
 
    depth = -ds_grid.Z[k_plot].values
    depthstr = str(depth) + ' m depth'
        
    variable_dict = {'WVEL': 'Vertical velocity', \
                    'Delta_u': r'$|\Delta \vec{u}|_n$', \
                    'ZETA': 'Vorticity', \
                    'zetanorm': r'Vorticity, normalized by $f_{mean}$', \
                    'OW': 'Okubo-Weiss parameter', \
                    's': 'Strain', \
                    'zeta_geos': r'Vorticity (computed from $\vec{u}_g$), normalized by $f_{mean}$', \
                    'OW_geos': r'Okubo-Weiss parameter (computed from $\vec{u}_g$)', \
                    'Ro_l': 'Local Rossby number', \
                    'geos_metric': r'Metric for geostrophy $\frac{||\vec{u} - \vec{u}_g||}{|\vec{u}|| + ||\vec{u}_g||}$', \
                    'PHIHYDcR': 'Hydrostatic pressure anomaly', \
                    'DIVU': 'Divergence of horizontal water velocity', \
                    'DIVUEk': 'Divergence of Ekman current'}
    variable_name = variable_dict[variable]
    
    title = variable_name + ' in BGR at {}, {} \n'.format(depthstr, datestr)
    
    return title

##############################

def get_quiver(ax, ds_grid, vector_ds, vector_comps, depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res, quiv_scale=0.3):
    
    """
    Resamples to lat-lon grid and gets quiver object given an ax.
    """
    
    #curr_ds_grid = ds_grid.load()

    lons, lats, lon_edges, lat_edges, vec_E_comp, vec_N_comp = vector_to_grid(ds_grid, vector_ds, vector_comps, depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res)
    
    #lons, lats, lon_edges, lat_edges, vec_x_comp = vector_to_grid(curr_ds_grid, vec_E_comp, latmin, latmax, lonmin, lonmax, resolution)
    #vec_y_comp = ecco_resample(curr_ds_grid, vec_N_comp, latmin, latmax, lonmin, lonmax, resolution)[4]
            
    skip = (slice(0, -1, 1), slice(0, -1, 1))
        
    quiv = ax.quiver(lons[skip], lats[skip], vec_E_comp[skip], vec_N_comp[skip], color='k', \
                     transform=ccrs.PlateCarree(), scale=quiv_scale, scale_units='width', regrid_shape=30)
    
    return quiv

##############################

def plot_geography(ax, labels=True):
    
    """
    Adds land, coastlines, and grid to an ax.
    """
    
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    lines = ax.gridlines(draw_labels=labels)
    
    if labels:
        
        lines.xlocator = mticker.FixedLocator(np.arange(-180,-90,10))
        lines.xformatter = LONGITUDE_FORMATTER
        lines.xpadding = 10
        lines.xlabel_style = {'size': 12, 'color': 'grey', 'rotation': 0}
        lines.ylocator = mticker.FixedLocator(np.arange(70,86,2))
        lines.yformatter = LATITUDE_FORMATTER
        lines.ylabel_style = {'size': 12, 'color': 'grey', 'rotation': 0}
        lines.x_inline = False
        lines.y_inline = False
        lines.xlabels_top = False
        lines.ylabels_right = False
        
        plt.draw()
        
    return ax

##############################

def get_scalar_bounds(scalar_field, scalar_bounds=None):
    
    """
    Set vmin, vmax for a scalar field.
    """
    
    if scalar_bounds == None:
        vmin, vmax = np.nanmin(scalar_field), np.nanmax(scalar_field) #Default to actual bounds on data
    else:
        vmin, vmax = scalar_bounds[0], scalar_bounds[1] #If scalar bounds are explicitly given
    
    return vmin, vmax

##############################

def get_pcolormesh(ax, lon_centers, lat_centers, scalar, cmap, vmin, vmax, logscale=False):
    
    """
    Create pcolormesh object given an axis.
    """
    
    if logscale:
        color = ax.pcolormesh(lon_centers, lat_centers, scalar, transform=ccrs.PlateCarree(), cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    elif not logscale:
        color = ax.pcolormesh(lon_centers, lat_centers, scalar, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
    
    return ax, color

##############################

def ArcCir_pcolormesh(scalar_field_name, date_string, datdir_primary, datdir_secondary, time_ave_type, plot_plane_type, spatial_bounds, resolutions, outfile, extend='both', vector_field_name=None):
    
    """
    Create pcolormesh plot of a scalar variable in a subdomain of the Arctic.
    Optional - overlay quiver plot of a vector variable.
    """
    
    ds_grid = load_grid(datdir_primary) #Load the ECCO grid
    
    if plot_plane_type == 'depth_const': #will add options 'latitude_const' and 'longitude_const'
        depth, latmin, latmax, lonmin, lonmax = spatial_bounds[0], spatial_bounds[1], spatial_bounds[2], spatial_bounds[3], spatial_bounds[4]
        lat_res, lon_res = resolutions[0], resolutions[1]
    
    #Load the scalar DataSet
    
    if field_is_primary(scalar_field_name): #Call function that loads ECCO data
        scalar_ds = load_ECCO_data_file(scalar_field_name, date_string, datdir_primary, time_ave_type)

    elif not field_is_primary(scalar_field_name): #Call function that loads computed data
        scalar_ds = load_comp_data_file(scalar_field_name, date_string, datdir_secondary, time_ave_type)

    lons, lats, lon_edges, lat_edges, scalar_field = scalar_to_grid(ds_grid, scalar_ds, get_field_variable(scalar_field_name), depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res)
    
    vmin, vmax = get_scalar_bounds(scalar_field) #Set scalar bounds

    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-135))    
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
        
    #Create pcolormesh object #tba - make a function for cmap
    ax, color = get_pcolormesh(ax, lons, lats, scalar_field, 'viridis', vmin, vmax) 
    
    if vector_field_name is not None:

        #Load the vector DataSet
        
        if field_is_primary(vector_field_name): #Call function that loads ECCO data
            vector_ds = load_ECCO_data_file(vector_field_name, date_string, datdir_primary, time_ave_type)

        elif not field_is_primary(vector_field_name): #Call function that loads computed data
            vector_ds = load_comp_data_file(vector_field_name, date_string, datdir_secondary, time_ave_type)
        
        #vec_E_comp = vector_to_grid(ds_grid, vector_ds, get_vector_comps(vector_field_name), depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res)[4]
        #vec_N_comp = vector_to_grid(ds_grid, vector_ds, get_vector_comps(vector_field_name), depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res)[5]
        
        quiv = get_quiver(ax, ds_grid, vector_ds, get_vector_comps(vector_field_name), depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res)
        #quiv = get_quiver(ax, ds_grid, vec_E_comp, vec_N_comp, latmin, latmax, lonmin, lonmax, resolution, quiv_scale)
        
    ax = plot_geography(ax)
    ax.set_title(pcolormesh_k_title(ds_grid, int(depth), get_field_variable(scalar_field_name), date_string)) #need to fix title function
    fig.colorbar((color), ax=ax, label=cbar_label(scalar_field_name), extend=extend, location='bottom')
    
    plt.savefig(outfile)
    plt.close()