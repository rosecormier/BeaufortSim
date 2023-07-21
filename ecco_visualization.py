"""
Rosalie Cormier, 2023
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py as ecco

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from ecco_general import get_month_name, get_scalar_in_xy, ds_to_field, comp_temp_mean, ecco_resample

plt.rcParams['font.size'] = 16
plt.rcParams['text.usetex'] = True

def cbar_label(scalar_attr):
    
    """
    Returns label for plot colorbar.
    
    scalar_attr = variable that the colorbar corresponds to
    """
    
    cbar_label_dict = {'PHIHYDcR': r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$', \
                      'WVEL': 'Velocity (m/s)', \
                      'Delta_u': r'$|\Delta \vec{u}|_n$', \
                      'ZETA': 'Vorticity (1/s)', \
                      'zetanorm': r'Vorticity per $f_{mean}$', \
                      'OW': r'OW $(1/s^2)$', \
                      's': r'Strain $(1/s^2)$', \
                      'zeta_geos': r'Vorticity per $f_{mean}$', \
                      'OW_geos': r'OW $(1/s^2)$', \
                      'Ro_l': r'$Ro_{\ell}$', \
                      'geos_metric': 'Velocity ratio'}
    label = cbar_label_dict[scalar_attr]
    
    return label

def pcolormesh_quiver_title(ecco_ds_grid, k_plot, datestr, scalar_attr, xvec_attr, resid=False):
    
    """
    Returns title for contourf-quiver plot.
    
    ecco_ds_grid = ECCO grid
    k_plot = k-value to plot at
    datestr = string containing date
    """
    
    ds_grid = ecco_ds_grid.copy()
   
    depth = - ds_grid.Z[k_plot].values
    depthstr = str(depth) + ' m depth'
        
    scalar_dict = {'PHIHYDcR': 'Pressure anomaly', \
                  'OW': 'Okubo-Weiss parameter', \
                  'OW_geos': r'Okubo-Weiss parameter (computed from $\vec{u}_g$)', \
                  'Ro_l': 'Local Rossby number', \
                  'ZETA': 'Vorticity'}
    scalar_str = scalar_dict[scalar_attr]
    
    vector_dict = {'UVEL': 'water velocity', \
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
                    'PHIHYDcR': 'Hydrostatic pressure anomaly'}
    variable_name = variable_dict[variable]
    
    title = variable_name + ' in BGR at {}, {} \n'.format(depthstr, datestr)
    
    return title

def get_quiver(ax, ecco_ds_grid, u_plot, v_plot, latmin, latmax, lonmin, lonmax, resolution, quiv_scale):
    
    """
    Resamples to lat-lon grid and gets quiver object given an ax.
    """
    
    ds_grid = ecco_ds_grid.copy()
    
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, u_nearest = ecco_resample(ds_grid, u_plot, latmin, latmax, lonmin, lonmax, resolution)
    v_nearest = ecco_resample(ds_grid, v_plot, latmin, latmax, lonmin, lonmax, resolution)[4]
            
    skip = (slice(0, -1, 1), slice(0, -1, 1))
        
    quiv = ax.quiver(new_grid_lon_centers[skip], new_grid_lat_centers[skip], u_nearest[skip], v_nearest[skip], color='k', \
                     transform=ccrs.PlateCarree(), scale=quiv_scale, scale_units='width', regrid_shape=30)
    
    return quiv

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

def ArcCir_pcolormesh(ecco_ds_grid, scalars, resolution, cmap, lon_centers, lat_centers, k_centers, datestr, scalar_attr='Delta_u', scalar_bounds=[1, 1], k_plot=None, extend='both', logscale=False, outfile="", lats_lons=[70.0, 85.0, -180.0, -90.0]):
    
    scalar = comp_temp_mean(scalars)
    
    #Set bounds on data

    if scalar_bounds[0] == scalar_bounds[1]:
        vmin, vmax = np.nanmin(scalar), np.nanmax(scalar)
    else:
        vmin, vmax = scalar_bounds[0], scalar_bounds[1]
        
    latmin, latmax, lonmin, lonmax = lats_lons #Spatial bounds
    
    if latmin != latmax and lonmin != lonmax: #If plotting a horizontal plane

        fig = plt.figure(figsize=(10, 8))
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-135))    
        ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
        
        if logscale:
            color = ax.pcolormesh(lon_centers, lat_centers, scalar, transform=ccrs.PlateCarree(), cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
        
        elif not logscale:
            color = ax.pcolormesh(lon_centers, lat_centers, scalar, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
    
        ax = plot_geography(ax)
        ax.set_title(pcolormesh_k_title(ecco_ds_grid, k_plot, scalar_attr, datestr))
        
    elif latmin != latmax and lonmin == lonmax: #If plotting along a line of constant longitude
        
        fig = plt.figure(figsize=(10, 8))
        ax = plt.axes()
        
        color = ax.pcolormesh(lat_centers, k_centers, scalar, cmap=cmap)
    
    fig.colorbar((color), ax=ax, label=cbar_label(scalar_attr), extend=extend, location='bottom')
    
    plt.savefig(outfile)
    plt.close()
    
    return scalar
    
def ArcCir_pcolormesh_quiver(ecco_ds_grid, k_plot, scalars, vecEs, vecNs, \
                           resolution, cmap, datestr, lon_centers, lat_centers, scalar_attr='PHIHYDcR', xvec_attr='UVEL', \
                           scalar_bounds=[1, 1], extend='both', logscale=False, outfile="", \
                           lats_lons=[70.0, 85.0, -175.5, -90.5], quiv_scale=0.3):
    
    """
    Creates pcolormesh plot of scalar variable in a subdomain of the Arctic,
        overlaid with quiver plot of vector variable.
    """

    scalar_mean = comp_temp_mean(scalars)
    scalar = scalar_mean
        
    vecE_mean, vecN_mean = comp_temp_mean(vecEs), comp_temp_mean(vecNs)
    vecE, vecN = vecE_mean, vecN_mean

    latmin, latmax, lonmin, lonmax = lats_lons

    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-135))
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
    
    if scalar_bounds[0] == scalar_bounds[1]:
        vmin, vmax = np.nanmin(scalar), np.nanmax(scalar)
    else:
        vmin, vmax = scalar_bounds[0], scalar_bounds[1]

    if logscale:
        color = ax.pcolormesh(lon_centers, lat_centers, scalar, transform=ccrs.PlateCarree(), cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
        
    elif not logscale:
        color = ax.pcolormesh(lon_centers, lat_centers, scalar, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
    
    quiv = get_quiver(ax, ecco_ds_grid, vecE, vecN, latmin, latmax, lonmin, lonmax, resolution, quiv_scale)
    
    ax = plot_geography(ax)
    ax.set_title(pcolormesh_quiver_title(ecco_ds_grid, k_plot, datestr, scalar_attr, xvec_attr))
        
    fig.colorbar((color), ax=ax, label=cbar_label(scalar_attr), extend=extend, location='bottom')
    
    plt.savefig(outfile)
    plt.close()
    
    return scalar_mean, vecE_mean, vecN_mean
"""
def ArcCir_contourf_quiver_grid(ecco_ds_grid, k_plot, scalars, vecEs, vecNs, \
                                resolution, cmap, yearstr, lon_centers, lat_centers, scalar_attr='PHIHYDcR', xvec_attr='UVEL', \
                                scalar_bounds=None, outfile="", lats_lons=[70.0, 85.0, -180.0, -90.0], \
                                no_levels=15, scale_factor=1, arrow_spacing=10, quiv_scale=10, resid=False):
    
    
    #Creates array of contourf plots of scalar variable in a subdomain of the Arctic,
    #    overlaid with quiver plot of vector variable.
    
    ecco_ds_grid = ECCO grid
    k_plot = depth index to plot at
    scalars = scalar DataSets
    vecEs, vecNs = east/northward vector component data
    scalar_bounds = bounds on scalar attribute
    resolution = resolution (both lat and lon) in degrees
    cmap = colormap name
    yearstr = string of year to plot
    outfile = output file name
    resid = whether this is a residual plot
    
    
    plt.rcParams['font.size'] = 40
    
    mainfig = plt.figure(figsize=(48, 40))
    
    months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    
    nrows, ncols = 3, 4
    nplots = nrows * ncols   
    row, col = -1, 0
        
    for i in range(nplots):
        
        monthstr = months[i]
        
        if i % ncols == 0:
            row += 1
            col = 0
        else:
            col += 1
            
        scalar, vecE, vecN = scalars[i], vecEs[i], vecNs[i]
        
        latmin, latmax, lonmin, lonmax = lats_lons

        ax = mainfig.add_subplot(nrows, ncols, i + 1, projection=ccrs.NorthPolarStereo(central_longitude=-135))
        ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
        
        if scalar_bounds is None:
            vmin, vmax = np.nanmin(scalar), np.nanmax(scalar)
        else:
            vmin, vmax = scalar_bounds[0], scalar_bounds[1]
            
        filled_contours, contour_lines = get_contours(ax, lon_centers, lat_centers, scalar, vmin, vmax, no_levels, cmap)

        quiv = get_quiver(ax, ecco_ds_grid, vecE, vecN, latmin, latmax, lonmin, lonmax, resolution, quiv_scale)
    
        plot_geography(ax, labels=False)
        ax.set_title('\n {} {}'.format(get_month_name(monthstr), yearstr))
        
    mainfig.suptitle(contourf_quiver_title(ecco_ds_grid, k_plot, yearstr, scalar_attr, xvec_attr, resid=resid), size=80)
    mainfig.tight_layout()
    cbar = mainfig.colorbar(filled_contours, ax=mainfig.get_axes(), aspect=40, pad=0.05, \
                            label=cbar_label(scalar_attr), location='bottom')
    mainfig.savefig(outfile)
    plt.close()
    
    plt.rcdefaults()
"""