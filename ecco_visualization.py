import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py as ecco

from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
from xgcm import Grid

from ecco_general import get_month_name, get_scalar_in_xy, ds_to_field, comp_temp_mean, ecco_resample

plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True

def cbar_label(scalar_attr):
    
    """
    Returns label for plot colorbar.
    
    scalar_attr = variable that the colorbar corresponds to
    """
    
    cbar_label_dict = {'PHIHYDcR': r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$', \
                      'Delta_u': r'$|\Delta u|_n$', \
                      'zeta': 'Vorticity (1/s)', \
                      'W': r'W $(1/s^2)$', \
                      's': r'Strain $(1/s^2)$'}
    label = cbar_label_dict[scalar_attr]
    
    return label

def contourf_quiver_title(ecco_ds_grid, k_plot, datestr, scalar_attr, xvec_attr, resid=False):
    
    """
    Returns title for contourf-quiver plot.
    
    ecco_ds_grid = ECCO grid
    k_plot = k-value to plot at
    datestr = string containing date
    """
    
    ds_grid = ecco_ds_grid.copy()
    
    if k_plot == 0:
        depthstr = 'ocean surface'
    
    elif k_plot != 0:
        depth = - ds_grid.Z[k_plot].values
        depthstr = str(depth) + ' m depth'
        
    scalar_dict = {'PHIHYDcR': 'Pressure anomaly'}
    scalar_str = scalar_dict[scalar_attr]
    
    vector_dict = {'UVEL': 'water velocity'}
    vector_str = vector_dict[xvec_attr]
    
    if resid:
        title = scalar_str + ' and ' + vector_str + \
            ' residuals (relative to annual mean) \n in Arctic Circle at {}, {} \n'.format(depthstr, \
                                                                                           datestr)

    else:
        title = scalar_str + ' and ' + vector_str + ' in Arctic Circle \n at {}, {} \n'.format(depthstr, \
                                                                                               datestr)
    
    return title

def pcolormesh_title(ds_grid, k_plot, variable, datestr):
    
    if k_plot == 0:
        depthstr = 'ocean surface'
    
    elif k_plot != 0:
        depth = - ds_grid.Z[k_plot].values
        depthstr = str(depth) + ' m depth'
        
    variable_dict = {'Delta_u': r'$|\Delta u|_n$', \
                    'zeta': 'Vorticity', \
                    'W': 'Okubo-Weiss parameter', \
                    's': 'Strain'}
    variable_name = variable_dict[variable]
    
    title = variable_name + ' in Arctic Circle at {}, {} \n'.format(depthstr, datestr)
    
    return title

def get_contours(ax, new_grid_lon_centers, new_grid_lat_centers, field, vmin, vmax, no_levels, cmap, lines=True):
    
    """
    Gets contour lines and fill objects given an ax.
    """
    
    filled_contours = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, field, \
                      levels=np.linspace(vmin, vmax, no_levels), \
                      transform=ccrs.PlateCarree(), extend='both', cmap=cmap)
    
    if lines:
        
        contour_lines = ax.contour(new_grid_lon_centers, new_grid_lat_centers, field, colors='r', \
                         alpha=0.8, linewidths=0.5, zorder=100, transform=ccrs.PlateCarree(), \
                         levels=np.linspace(vmin, vmax, no_levels))
    
        return filled_contours, contour_lines
    
    else:
        return filled_contours

def get_quiver(ax, ecco_ds_grid, u_plot, v_plot, latmin, latmax, lonmin, lonmax, resolution):
    
    """
    Resamples to lat-lon grid and gets quiver object given an ax.
    """
    
    ds_grid = ecco_ds_grid.copy()
    
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, u_nearest = ecco_resample(ds_grid, u_plot, latmin, latmax, lonmin, lonmax, resolution)
    v_nearest = ecco_resample(ds_grid, v_plot, latmin, latmax, lonmin, lonmax, resolution)[4]
            
    quiv = ax.quiver(new_grid_lon_centers, new_grid_lat_centers, u_nearest, v_nearest, color='k', \
                     transform=ccrs.PlateCarree(), scale=1, regrid_shape=60, zorder=150)
    
    return quiv

def plot_geography(ax):
    
    """
    Adds land, coastlines, and grid to an ax.
    """
    
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    ax.gridlines()
    
    return ax

def ArcCir_pcolormesh(ecco_ds_grid, k_plot, scalars, resolution, cmap, lon_centers, lat_centers, datestr, scalar_attr='Delta_u', scalar_bounds=None, outfile="", lats_lons=[70.0, 85.0, -180.0, -90.0]):
    
    scalar_mean = comp_temp_mean(scalars)
    scalar = scalar_mean
        
    latmin, latmax, lonmin, lonmax = lats_lons

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())    
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
    
    if scalar_bounds is None:
        vmin, vmax = np.nanmin(scalar), np.nanmax(scalar)
    else:
        vmin, vmax = scalar_bounds[0], scalar_bounds[1]
        
    color = ax.pcolormesh(lon_centers, lat_centers, scalar, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
    
    ax = plot_geography(ax)
    ax.set_title(pcolormesh_title(ecco_ds_grid, k_plot, scalar_attr, datestr))
    
    cbar = fig.colorbar((color), label=cbar_label(scalar_attr), extend='both')
    
    plt.savefig(outfile)
    plt.close()
    
    return scalar_mean
    
def ArcCir_contourf_quiver(ecco_ds_grid, k_plot, scalars, vecEs, vecNs, \
                           resolution, cmap, datestr, lon_centers, lat_centers, lon_edges, lat_edges, scalar_attr='PHIHYDcR', xvec_attr='UVEL', scalar_bounds=None, outfile="", \
                           lats_lons=[70.0, 85.0, -180.0, -90.0], no_levels=30, scale_factor=1, \
                           arrow_spacing=10, quiv_scale=1):
    
    """
    Creates contourf plot of scalar variable in a subdomain of the Arctic,
        overlaid with quiver plot of vector variable.
    
    ecco_ds_grid = ECCO grid
    k_plot = depth index to plot at
    scalars = scalar data
    vecEs, vecNs = east/northward vector component data
    resolution = resolution (both lat and lon) in degrees
    cmap = colormap name
    scalar_bounds = min/max of scalar data
    datestr = date string to be used in plot title
    """

    scalar_mean = comp_temp_mean(scalars)
    scalar = scalar_mean
        
    vecE_mean, vecN_mean = comp_temp_mean(vecEs), comp_temp_mean(vecNs)
    vecE, vecN = vecE_mean, vecN_mean

    latmin, latmax, lonmin, lonmax = lats_lons

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
    
    if scalar_bounds is None:
        vmin, vmax = np.nanmin(scalar), np.nanmax(scalar)
    else:
        vmin, vmax = scalar_bounds[0], scalar_bounds[1]
        
    filled_contours, contour_lines = get_contours(ax, lon_centers, lat_centers, scalar, vmin, vmax, no_levels, cmap)
    
    quiv = get_quiver(ax, ecco_ds_grid, vecE, vecN, latmin, latmax, lonmin, lonmax, resolution)
    
    plot_geography(ax)
    ax.set_title(contourf_quiver_title(ecco_ds_grid, k_plot, datestr, scalar_attr, xvec_attr))
        
    cbar = fig.colorbar(filled_contours, label=cbar_label(scalar_attr))
    
    plt.savefig(outfile)
    plt.close()
    
    return scalar_mean, vecE_mean, vecN_mean

def ArcCir_contourf_quiver_grid(ecco_ds_grid, k_plot, scalars, vecEs, vecNs, \
                                resolution, cmap, monthstrs, \
                                yearstrs, lon_centers, lat_centers, lon_edges, lat_edges, scalar_attr='PHIHYDcR', xvec_attr='UVEL', scalar_bounds=None, outfile="", lats_lons=[70.0, 85.0, -180.0, -90.0], \
                                no_levels=15, scale_factor=1, arrow_spacing=10, quiv_scale=0.3, \
                                nrows=3, ncols=4, resid=False):
    
    """
    Creates array of contourf plots of scalar variable in a subdomain of the Arctic,
        overlaid with quiver plot of vector variable.
    
    ecco_ds_grid = ECCO grid
    k_plot = depth index to plot at
    scalars = scalar DataSets
    vecEs, vecNs = east/northward vector component data
    scalar_bounds = bounds on scalar attribute
    resolution = resolution (both lat and lon) in degrees
    cmap = colormap name
    month/yearstrs = strings of months/years to plot
    outfile = output file name
    resid = whether this is a residual plot
    """
    
    plt.rcParams['font.size'] = 40
    
    mainfig = plt.figure(figsize=(48, 40))
    
    nplots = nrows * ncols   
    row, col = -1, 0
        
    for i in range(nplots):
        
        if i % ncols == 0:
            row += 1
            col = 0
        else:
            col += 1
            
        scalar, vecE, vecN = scalars[i], vecEs[i], vecNs[i]
        monthstr, yearstr = monthstrs[i], yearstrs[i]
        
        latmin, latmax, lonmin, lonmax = lats_lons

        ax = mainfig.add_subplot(nrows, ncols, i + 1, projection=ccrs.NorthPolarStereo())
        ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
        
        if scalar_bounds is None:
            vmin, vmax = np.nanmin(scalar), np.nanmax(scalar)
        else:
            vmin, vmax = scalar_bounds[0], scalar_bounds[1]
            
        filled_contours, contour_lines = get_contours(ax, lon_centers, lat_centers, scalar, vmin, vmax, no_levels, cmap)

        quiv = get_quiver(ax, ecco_ds_grid, vecE, vecN, latmin, latmax, lonmin, lonmax, resolution)
    
        plot_geography(ax)
        ax.set_title('\n {} {}'.format(get_month_name(monthstr), yearstr))
        
    mainfig.suptitle(contourf_quiver_title(ecco_ds_grid, k_plot, yearstrs[0], scalar_attr, xvec_attr, resid=resid), \
                     size=80)
    mainfig.tight_layout()
    cbar = mainfig.colorbar(filled_contours, ax=mainfig.get_axes(), aspect=40, pad=0.05, \
                            label=cbar_label(scalar_attr), location='bottom')
    mainfig.savefig(outfile)
    plt.close()
    
    plt.rcdefaults()