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

import load_data_files

from functions_ecco_general import load_grid, scalar_to_grid, vector_to_grid
from functions_field_variables import get_field_variable, get_vector_comps, field_is_primary, get_cmap_and_symmetry

##############################

plt.rcParams['font.size'] = 16
plt.rcParams['text.usetex'] = True

##############################

def get_cbar_label(scalar_field_name):
    
    """
    Return label for plot colorbar.
    """

    cbar_label_dict = {'pressure': r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$', \
                       'density': r'Density anomaly $(kg/{m}^3)$', \
                       'vertical_vel': 'Velocity (m/s)', \
                       'vorticity': 'Vorticity (1/s)', \
                       'normal_strain': r'Normal strain $(1/s^2)$', \
                       'shear_strain': r'Shear strain $(1/s^2)$', \
                       '2D_div_vel': 'Horizontal velocity divergence (1/s)'}
    
    label = cbar_label_dict[scalar_field_name]
    
    return label

##############################

def get_plot_title(scalar_field_name, vector_field_name, plot_plane_type, spatial_bounds, ds_grid, date_string):
 
    """
    Return main title for plot.
    """

    if plot_plane_type == 'depth_const':
        k_val = int(spatial_bounds[0])
        depth = -ds_grid.Z[k_val].values
        depth_string = str(depth) + ' m depth'
        
    field_titles = {'RHOAnoma': 'Density Anomaly', \
                    'PHIHYDcR': 'Hydrostatic Pressure Anomaly', \
                    'WVEL': 'Vertical Velocity', \
                    'ZETA': 'Vorticity', \
                    'NORMAL': 'Normal Strain', \
                    'SHEAR': 'Shear Strain', \
                    'DIVU': 'Divergence of Horizontal Velocity', \
                    'UVELVVEL': 'Horizontal Velocity', \
                    'UGVG': 'Geostrophic Velocity', \
                    'EXFtauxEXFtauy': 'Surface Wind-on-Ocean Stress', \
                    'UEkVEk': 'Ekman Velocity'}

    scalar_field_title = field_titles[get_field_variable(scalar_field_name)]
    
    if vector_field_name is None:
        title = '{} in BGR at {}, {} \n'.format(scalar_field_title, depth_string, date_string)

    elif vector_field_name is not None:
        vector_field_title = field_titles[get_field_variable(vector_field_name)]
        title = '{} and {} in BGR \n at {}, {} \n'.format(scalar_field_title, vector_field_title, depth_string, date_string)
        
    return title

##############################

def plot_geography(ax, labels=True):
    
    """
    Add land, coastlines, and grid to an ax.
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

def get_pcolormesh(ax, lon_centers, lat_centers, scalar, field_name, vmin, vmax, logscale=False):
    
    """
    Create pcolormesh object given an axis.
    """
    
    #Get appropriate colormap and symmetry (boolean)
    cmap, symmetry = get_cmap_and_symmetry(field_name)
    
    if symmetry: #If variable is symmetric about zero, reset bounds on colorbar to be symmetric
        abs_max_value = max(abs(vmax), abs(vmin))
        vmin, vmax = -abs_max_value, abs_max_value
        
    if logscale:
        color = ax.pcolormesh(lon_centers, lat_centers, scalar, transform=ccrs.PlateCarree(), cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    elif not logscale:
        color = ax.pcolormesh(lon_centers, lat_centers, scalar, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
    
    return ax, color

##############################

def get_quiver(ax, lon_centers, lat_centers, vec_E_comp, vec_N_comp, quiv_scale=0.3):
    
    """
    Resample vector field to lat-lon grid and get quiver object; add quiver to given ax.
    """

    skip = (slice(0, -1, 1), slice(0, -1, 1))
    quiv = ax.quiver(lon_centers[skip], lat_centers[skip], vec_E_comp[skip], vec_N_comp[skip], color='k', \
                     transform=ccrs.PlateCarree(), scale=quiv_scale, scale_units='width', regrid_shape=30)
    
    return quiv

##############################

def ArcCir_pcolormesh(scalar_field_name, date_string, datdir_primary, datdir_secondary, time_ave_type, plot_plane_type, spatial_bounds, resolutions, outfile, extend='neither', vector_field_name=None):
    
    """
    Create pcolormesh plot of a scalar variable in a subdomain of the Arctic.
    Optional - overlay quiver plot of a vector variable.
    """
    
    ds_grid = load_grid(datdir_primary) #Load the ECCO grid

    #Load scalar DataSet
    scalar_ds = load_data_files.main(field_name=scalar_field_name, date_string=date_string, datdir_primary=datdir_primary, datdir_secondary=datdir_secondary, time_ave_type=time_ave_type)
    
    if plot_plane_type == 'depth_const': #will add options 'latitude_const' and 'longitude_const'
        depth, latmin, latmax, lonmin, lonmax = spatial_bounds[0], spatial_bounds[1], spatial_bounds[2], spatial_bounds[3], spatial_bounds[4]
        lat_res, lon_res = resolutions[0], resolutions[1]
        
    lons, lats, lon_edges, lat_edges, scalar_field = scalar_to_grid(ds_grid, scalar_ds, get_field_variable(scalar_field_name), depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res)
    
    vmin, vmax = get_scalar_bounds(scalar_field) #Set scalar bounds

    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-135))    
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
        
    ax, color = get_pcolormesh(ax, lons, lats, scalar_field, scalar_field_name, vmin, vmax) #Create pcolormesh object 
    
    if vector_field_name is not None:
       
        #Load vector DataSet
        vector_ds = load_data_files.main(field_name=vector_field_name, date_string=date_string, datdir_primary=datdir_primary, datdir_secondary=datdir_secondary, time_ave_type=time_ave_type)
        
        lons, lats, lon_edges, lat_edges, vec_E_comp, vec_N_comp = vector_to_grid(ds_grid, vector_ds, vector_field_name, depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res)
        quiv = get_quiver(ax, lons, lats, vec_E_comp, vec_N_comp) #Create quiver object 
   
    ax = plot_geography(ax)
    ax.set_title(get_plot_title(scalar_field_name, vector_field_name, plot_plane_type, spatial_bounds, ds_grid, date_string))
    
    fig.colorbar((color), ax=ax, label=get_cbar_label(scalar_field_name), extend=extend, location='bottom')
    
    plt.savefig(outfile)
    plt.close()