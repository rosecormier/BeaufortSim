"""
Auxiliary functions for visualizing ECCO data.

Rosalie Cormier, 2023
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

#from functions_divergence import comp_2d_Ek_divergence
#from functions_ecco_general import get_month_name, get_scalar_in_xy, ds_to_field, comp_temp_mean, ecco_resample, load_dataset, get_vector_partner
#from functions_field_variables import get_field_vars, get_variable_str
#from functions_vorticity import get_OW_field, comp_local_Ro

from functions_comp_data_meta import load_comp_data_file
from functions_ecco_general import load_ECCO_data_file, load_grid, ds_to_field
from functions_field_variables import get_field_variable, field_is_primary

##############################

plt.rcParams['font.size'] = 16
plt.rcParams['text.usetex'] = True

##############################

def cbar_label(scalar_attr):
    
    """
    Returns label for plot colorbar.
    """
    #add density; fix others
    cbar_label_dict = {'pressure': r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$', \
                      'WVEL': 'Velocity (m/s)', \
                      'Delta_u': r'$|\Delta \vec{u}|_n$', \
                      'ZETA': 'Vorticity (1/s)', \
                      'zetanorm': r'Vorticity per $f_{mean}$', \
                      'OW': r'OW $(1/s^2)$', \
                      's': r'Strain $(1/s^2)$', \
                      'zeta_geos': r'Vorticity per $f_{mean}$', \
                      'OW_geos': r'OW $(1/s^2)$', \
                      'Ro_l': r'$Ro_{\ell}$', \
                      'geos_metric': 'Velocity ratio', \
                      'DIVU': 'Horizontal velocity divergence (1/s)', \
                      'DIVUEk': 'Divergence of Ekman current (1/s)'}
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
 
    depth = - ds_grid.Z[k_plot].values
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

#def get_quiver(ax, ecco_ds_grid, u_plot, v_plot, latmin, latmax, lonmin, lonmax, resolution, quiv_scale):
    
#    """
#    Resamples to lat-lon grid and gets quiver object given an ax.
#    """
"""    
    ds_grid = ecco_ds_grid.copy()
    
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, u_nearest = ecco_resample(ds_grid, u_plot, latmin, latmax, lonmin, lonmax, resolution)
    v_nearest = ecco_resample(ds_grid, v_plot, latmin, latmax, lonmin, lonmax, resolution)[4]
            
    skip = (slice(0, -1, 1), slice(0, -1, 1))
        
    quiv = ax.quiver(new_grid_lon_centers[skip], new_grid_lat_centers[skip], u_nearest[skip], v_nearest[skip], color='k', \
                     transform=ccrs.PlateCarree(), scale=quiv_scale, scale_units='width', regrid_shape=30)
    
    return quiv
"""
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

def ArcCir_pcolormesh(scalar_field_name, date_string, datdir_primary, time_ave_type, plot_plane_type, spatial_bounds, resolutions, outfile, extend='both'):
    
    """
    Creates pcolormesh plot of a scalar variable in a subdomain of the Arctic.
    """
    
    ds_grid = load_grid(datdir_primary) #Load the ECCO grid
    
    #Load the scalar DataSet
    
    if field_is_primary(scalar_field_name): #Call function that loads ECCO data
        scalar_ds = load_ECCO_data_file(scalar_field_name, date_string, datdir_primary, time_ave_type)

    elif not field_is_primary(scalar_field_name): #Call function that loads computed data
        scalar_ds = load_comp_data_file(scalar_field_name, date_string, datdir_secondary, time_ave_type)
        
    if plot_plane_type == 'depth_const': #will add options 'latitude_const' and 'longitude_const'
        depth, latmin, latmax, lonmin, lonmax = spatial_bounds[0], spatial_bounds[1], spatial_bounds[2], spatial_bounds[3], spatial_bounds[4]
        lat_res, lon_res = resolutions[0], resolutions[1]

    lons, lats, lon_edges, lat_edges, scalar_field = ds_to_field(ds_grid, scalar_ds.isel(k=int(depth)), get_field_variable(scalar_field_name), latmin, latmax, lonmin, lonmax, lat_res, lon_res)
    
    vmin, vmax = get_scalar_bounds(scalar_field) #Set scalar bounds

    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-135))    
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
        
    #Create pcolormesh object #tba - make a function for cmap
    ax, color = get_pcolormesh(ax, lons, lats, scalar_field, 'viridis', vmin, vmax) 
        
    ax = plot_geography(ax)
    ax.set_title(pcolormesh_k_title(ds_grid, int(depth), get_field_variable(scalar_field_name), date_string))
    
    fig.colorbar((color), ax=ax, label=cbar_label(scalar_field_name), extend=extend, location='bottom')
    
    plt.savefig(outfile)
    plt.close()
    
##############################
"""
def ArcCir_pcolormesh_quiver(ecco_ds_grid, k_plot, scalars, vecEs, vecNs, \
                           resolution, cmap, datestr, lon_centers, lat_centers, scalar_attr='PHIHYDcR', xvec_attr='UVEL', \
                           scalar_bounds=[1, 1], extend='both', logscale=False, outfile="", \
                           lats_lons=[70.0, 85.0, -175.5, -90.5], quiv_scale=0.3):
    
"""
    #Creates pcolormesh plot of scalar variable in a subdomain of the Arctic,
    #    overlaid with quiver plot of vector variable.
"""

    scalar_mean = comp_temp_mean(scalars)
    scalar = scalar_mean
        
    vecE_mean, vecN_mean = comp_temp_mean(vecEs), comp_temp_mean(vecNs)
    vecE, vecN = vecE_mean, vecN_mean

    latmin, latmax, lonmin, lonmax = lats_lons

    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-135))
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
    
    vmin, vmax = get_scalar_bounds(scalar_bounds, scalar) #Set scalar bounds
    
    #Create pcolormesh object
    ax, color = get_pcolormesh(ax, lon_centers, lat_centers, scalar, cmap, vmin, vmax, logscale=logscale) 
   
    quiv = get_quiver(ax, ecco_ds_grid, vecE, vecN, latmin, latmax, lonmin, lonmax, resolution, quiv_scale)
    
    ax = plot_geography(ax)
    ax.set_title(pcolormesh_quiver_title(ecco_ds_grid, k_plot, datestr, scalar_attr, xvec_attr))
        
    fig.colorbar((color), ax=ax, label=cbar_label(scalar_attr), extend=extend, location='bottom')
    
    plt.savefig(outfile)
    plt.close()
    
    return scalar_mean, vecE_mean, vecN_mean
"""