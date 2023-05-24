import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py as ecco

from matplotlib.gridspec import GridSpec
from xgcm import Grid

from ecco_general import comp_temp_mean_scalar, comp_temp_mean_vector

plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True
vir_nanmasked = plt.get_cmap('viridis_r').copy()
vir_nanmasked.set_bad('black')

def cbar_label(scalar_attr):
    
    """
    Returns label for plot colorbar.
    
    scalar_attr = variable that the colorbar corresponds to
    """
    
    cbar_label_dict = {'PHIHYDcR': r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$'}
    label = cbar_label_dict[scalar_attr]
    
    return label

def contourf_quiver_title(ecco_ds_grid, k_plot, datestr, scalar_attr, xvec_attr, resid=False):
    
    """
    Returns title for contourf-quiver plot.
    
    ecco_ds_grid = ECCO grid
    k_plot = k-value to plot at
    datestr = string containing date
    """
    
    if k_plot == 0:
        depthstr = 'ocean surface'
    
    elif k_plot != 0:
        depth = - ecco_ds_grid.Z[k_plot].values
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

def ArcCir_contourf_quiver(ecco_ds_grid, k_val, ecco_ds_scalars, ecco_ds_vectors, scalar_attr, \
                           xvec_attr, yvec_attr, resolution, cmap, scalar_bounds, datestr, outfile="", \
                           lats_lons=[70.0, 85.0, -180.0, -90.0], no_levels=30, scale_factor=1, \
                           arrow_spacing=10, quiv_scale=1):
    
    """
    Creates contourf plot of scalar variable in a subdomain of the Arctic,
        overlaid with quiver plot of vector variable.
    
    ecco_ds_grid = ECCO grid
    k_plot = depth index to plot at
    ecco_ds_scalar = scalar DataSet(s)
    ecco_ds_vector = vector DataSet(s)
    scalar_attr = string corresponding to scalar attribute to plot
    xvec_attr = string corresponding to x-comp of vector to plot
    yvec_attr = string corresponding to y-comp of vector to plot
    resolution = resolution (both lat and lon) in degrees
    cmap = colormap name
    scalar_bounds = min/max of scalar data
    datestr = date string to be used in plot title
    outfile = output file name
    lats_lons = latmin, latmax, lonmin, lonmax
    no_levels = number of contour levels
    scale_factor = colorbar multiplier
    arrow_spacing = quiver arrow spacing in gridpoints
    quiv_scale = quiver plot scale
    """
    
    skip_k_scalar, skip_k_vector = False, False
    
    if len(ecco_ds_scalars) == 1:
        ecco_ds_scalar = ecco_ds_scalars[0]
        
    elif len(ecco_ds_scalars) > 1:
        ecco_ds_scalar, skip_k_scalar = comp_temp_mean_scalar(k_val, ecco_ds_scalars, scalar_attr)
        scalar_mean = ecco_ds_scalar
        
    if len(ecco_ds_vectors) == 1:
        ecco_ds_vector = ecco_ds_vectors[0]
        
    elif len(ecco_ds_vectors) > 1:
        ecco_ds_vector, skip_k_vector = comp_temp_mean_vector(k_val, ecco_ds_vectors, xvec_attr, yvec_attr)
        vector_mean = ecco_ds_vector
    
    ds_grid = get_scalar_in_xy(ecco_ds_grid, k_val, ecco_ds_scalar, scalar_attr, skip_k=skip_k_scalar)
    curr_field = (ds_grid[scalar_attr]).squeeze()
    
    ds_grid = ecco_ds_grid.copy()

    velc = get_vector_in_xy(ecco_ds_grid, k_val, ecco_ds_vector, xvec_attr, yvec_attr, skip_k=skip_k_vector)
    velE = velc['X'] * ds_grid['CS'] - velc['Y'] * ds_grid['SN']
    velN = velc['X'] * ds_grid['SN'] + velc['Y'] * ds_grid['CS']
    
    delta_lat, delta_lon = resolution, resolution
    latmin, latmax, lonmin, lonmax = lats_lons

    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, \
    field_nearest = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, curr_field, latmin, latmax, delta_lat, \
                                            lonmin, lonmax, delta_lon, fill_value=np.NaN, \
                                            mapping_method='nearest_neighbor', radius_of_influence=120000)
    
    vmin, vmax = scalar_bounds[0], scalar_bounds[1]

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
    
    cs1 = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, field_nearest, \
                      levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)), no_levels), \
                      transform=ccrs.PlateCarree(), extend='both', cmap=cmap)
    cs2 = ax.contour(new_grid_lon_centers, new_grid_lat_centers, field_nearest, colors='r', \
                     alpha=0.8, linewidths=0.5, zorder=100, transform=ccrs.PlateCarree(), \
                     levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)), no_levels))
    
    u_plot, v_plot = (velE).squeeze(), (velN).squeeze()

    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, u_nearest = \
    ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, u_plot, latmin, latmax, delta_lat, lonmin, \
                            lonmax, delta_lon, fill_value=np.NaN, mapping_method='nearest_neighbor', \
                            radius_of_influence=120000)
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, v_nearest = \
    ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, v_plot, latmin, latmax, delta_lat, lonmin, \
                            lonmax, delta_lon, fill_value=np.NaN, mapping_method='nearest_neighbor', \
                            radius_of_influence=120000)
            
    quiv = ax.quiver(new_grid_lon_centers, new_grid_lat_centers, u_nearest, v_nearest, color='k', \
                     transform=ccrs.PlateCarree(), scale=1, regrid_shape=60, zorder=150)
    
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    ax.gridlines()
    ax.set_title(contourf_quiver_title(ds_grid, k_val, datestr, scalar_attr, xvec_attr))
        
    cbar = fig.colorbar(cs1, ticks=range(int(np.floor(vmin)), int(np.ceil(vmax)), 1), \
                        label=cbar_label(scalar_attr))
    
    plt.savefig(outfile)
    plt.close()
    
    if len(ecco_ds_scalars) > 1 and len(ecco_ds_scalars) > 1:
        return scalar_mean, vector_mean

def ArcCir_contourf_quiver_grid(ecco_ds_grid, k_plot, ecco_ds_scalars, ecco_ds_vectors, scalar_attr, \
                                scalar_bounds, xvec_attr, yvec_attr, resolution, cmap, monthstrs, \
                                yearstrs, outfile="", lats_lons=[70.0, 85.0, -180.0, -90.0], \
                                no_levels=15, scale_factor=1, arrow_spacing=10, quiv_scale=0.3, \
                                nrows=3, ncols=4, resid=False):
    
    """
    Creates array of contourf plots of scalar variable in a subdomain of the Arctic,
        overlaid with quiver plots of vector variable.
    
    ecco_ds_grid = ECCO grid
    k_plot = depth index to plot at
    ecco_ds_scalars = scalar DataSets
    ecco_ds_vectors = vector DataSets
    scalar_attr = string corresponding to scalar attribute to plot
    scalar_bounds = bounds on scalar attribute
    xvec/yvec_attr = string corresponding to x/y-comp of vector to plot
    resolution = resolution (both lat and lon) in degrees
    cmap = colormap name
    month/yearstrs = strings of months/years to plot
    outfile = output file name
    lats_lons = latmin, latmax, lonmin, lonmax
    no_levels = number of contour levels
    scale_factor = colorbar multiplier
    arrow_spacing = quiver arrow spacing in gridpoints
    quiv_scale = quiver plot scale
    nrows, ncols = number of rows and columns (resp.) in grid
    resid = whether this is a residual plot
    """
    
    plt.rcParams['font.size'] = 40
    
    skip_k_scalar, skip_k_vector = False, False
    
    monthnames = {"01": "January", "02": "February", "03": "March", "04": "April", "05": "May", \
                  "06": "June", "07": "July", "08": "August", "09": "September", "10": "October", \
                  "11": "November", "12": "December"}
    
    mainfig = plt.figure(figsize=(48, 40))
    
    nplots = nrows * ncols   
    row, col = -1, 0
        
    for i in range(nplots):
        
        if i % ncols == 0:
            row += 1
            col = 0
        else:
            col += 1
            
        ecco_ds_scalar = ecco_ds_scalars[i]
        ecco_ds_vector = ecco_ds_vectors[i]
        monthstr = monthstrs[i]
        yearstr = yearstrs[i]
        
        ds_grid = get_scalar_in_xy(ecco_ds_grid, k_plot, ecco_ds_scalar, scalar_attr, skip_k=skip_k_scalar)
        curr_field = (ds_grid[scalar_attr]).squeeze()
        
        velc = get_vector_in_xy(ecco_ds_grid, k_plot, ecco_ds_vector, xvec_attr, yvec_attr)
        velE = velc['X'] * ds_grid['CS'] - velc['Y'] * ds_grid['SN']
        velN = velc['X'] * ds_grid['SN'] + velc['Y'] * ds_grid['CS']

        delta_lat, delta_lon = resolution, resolution
        latmin, latmax, lonmin, lonmax = lats_lons

        new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, \
        field_nearest = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, curr_field, latmin, \
                                                latmax, delta_lat, lonmin, lonmax, delta_lon, \
                                                fill_value=np.NaN, mapping_method='nearest_neighbor', \
                                                radius_of_influence=120000)

        ax = mainfig.add_subplot(nrows, ncols, i + 1, projection=ccrs.NorthPolarStereo())
        ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
        
        vmin, vmax = scalar_bounds[0], scalar_bounds[1]
        
        cs1 = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, field_nearest, \
                          levels=np.linspace(vmin, vmax, no_levels), transform=ccrs.PlateCarree(), \
                          extend='both', cmap=cmap)
        cs2 = ax.contour(new_grid_lon_centers, new_grid_lat_centers, field_nearest, colors='r', alpha=0.8, \
                         linewidths=1.0, zorder=100, transform=ccrs.PlateCarree(), \
                         levels=np.linspace(vmin, vmax, no_levels))

        u_plot, v_plot = (velE).squeeze(), (velN).squeeze()

        new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, u_nearest = \
        ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, u_plot, latmin, latmax, delta_lat, lonmin, lonmax, \
                                delta_lon, fill_value=np.NaN, mapping_method='nearest_neighbor', \
                                radius_of_influence=120000)
        new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, v_nearest = \
        ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, v_plot, latmin, latmax, delta_lat, lonmin, lonmax, \
                                delta_lon, fill_value=np.NaN, mapping_method='nearest_neighbor', \
                                radius_of_influence=120000)

        quiv = ax.quiver(new_grid_lon_centers, new_grid_lat_centers, u_nearest, v_nearest, color='k', \
                         transform=ccrs.PlateCarree(), scale=1, regrid_shape=30, zorder=150)

        ax.add_feature(cfeature.LAND)
        ax.coastlines()
        ax.gridlines()
        ax.set_title('\n {} {}'.format(monthnames[monthstr], yearstr))
        
    mainfig.suptitle(contourf_quiver_title(ds_grid, k_plot, yearstrs[0], scalar_attr, xvec_attr, resid=True), \
                     size=80)
    mainfig.tight_layout()
    cbar = mainfig.colorbar(cs1, ax=mainfig.get_axes(), aspect=40, pad=0.05, ticks=range(vmin, vmax, 1), \
                            label=cbar_label(scalar_attr), location='bottom')
    mainfig.savefig(outfile)
    plt.close()
    
    plt.rcdefaults()