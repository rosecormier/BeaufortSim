import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py as ecco

from matplotlib.gridspec import GridSpec
from xgcm import Grid

plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True
vir_nanmasked = plt.get_cmap('viridis_r').copy()
vir_nanmasked.set_bad('black')

def get_scalar_in_xy(ecco_ds_grid, k_val, ecco_ds_scalar, scalar_attr, skip_k=False):
    
    """
    Loads scalar field in xy-grid.
    
    ecco_ds_grid = ECCO grid
    k_val = depth index of interest
    ecco_ds_scalar = DataSet containing field
    scalar_attr = name of field
    skip_k = whether to skip isolating a value of k (True if value has already been isolated)
    """
    
    ds_grid = ecco_ds_grid.copy()

    if skip_k:
        ds_grid[scalar_attr] = ecco_ds_scalar[scalar_attr]
    
    else:
        ecco_ds_scalar_k = ecco_ds_scalar.isel(k=k_val)
        ds_grid[scalar_attr] = ecco_ds_scalar_k[scalar_attr]
        
    ds_grid = ds_grid.load()
    
    return ds_grid
    
def get_vector_in_xy(ecco_ds_grid, k_val, ecco_ds_vector, xvec_attr, yvec_attr, skip_k=False):
    
    """
    Loads vector field in xy-grid.
    
    ecco_ds_grid = ECCO grid
    k_val = depth index of interest
    ecco_ds_vector = DataSet containing vector field
    xvec_attr = name of x-comp of vector field
    yvec_attr = name of y-comp of vector field
    skip_k = whether to skip isolating a value of k (True if already isolated)
    """

    ds_grid = ecco_ds_grid.copy()
    ds_grid = ds_grid.load()
    XGCM_grid = ecco.get_llc_grid(ds_grid)
    
    if skip_k:
        velc = XGCM_grid.interp_2d_vector({'X': ecco_ds_vector[xvec_attr], \
                                           'Y': ecco_ds_vector[yvec_attr]}, \
                                           boundary='fill')
                                           
    else:
        velc = XGCM_grid.interp_2d_vector({'X': (ecco_ds_vector[xvec_attr]).isel(k=k_val), \
                                                 'Y': (ecco_ds_vector[yvec_attr]).isel(k=k_val)}, \
                                                 boundary='fill')

    return velc

def comp_temp_mean_scalar(k_val, ecco_ds_scalars, scalar_attr):
    
    """
    Computes temporal mean of a scalar field.

    k_val = depth value of interest
    ecco_ds_scalars = scalar Datasets
    scalar_attr = string corresponding to scalar attribute of interest
    """ 
    
    mean_field = ((ecco_ds_scalars[0]).copy()) * 0 
    concat_field = ((ecco_ds_scalars[0]).copy()).isel(k=k_val) * 0
    
    for dataset in ecco_ds_scalars:
        curr_field = dataset.isel(k=k_val)
        concat_field = xr.concat((concat_field, curr_field), dim='time')

    mean_scalar_field = (concat_field[scalar_attr]).sum(dim=['time']) / len(ecco_ds_scalars)
    mean_scalar_field = mean_scalar_field.where(mean_scalar_field != 0)
    mean_field[scalar_attr] = mean_scalar_field
    (mean_field[scalar_attr]).data = (mean_field[scalar_attr]).values

    skip_k = True
    
    return mean_field, skip_k

def comp_temp_mean_vector(k_val, ecco_ds_vectors, xvec_attr, yvec_attr):
    
    """
    Computes temporal mean of a vector field.
    
    k_val = depth value of interest
    ecco_ds_vectors = vector Datasets
    xvec_attr = string corresponding to x-comp of attribute of interest
    yvec_attr = string corresponding to y-comp of attribute of interest
    """
    
    mean_fields = ((ecco_ds_vectors[0]).copy()) * 0
    concat_x_field = ((ecco_ds_vectors[0]).copy()).isel(k=k_val) * 0
    concat_y_field = ((ecco_ds_vectors[0]).copy()).isel(k=k_val) * 0
    
    for dataset in ecco_ds_vectors:
        
        curr_x_field = dataset.isel(k=k_val)
        curr_y_field = dataset.isel(k=k_val)
        
        concat_x_field = xr.concat((concat_x_field, curr_x_field), dim='time')
        concat_y_field = xr.concat((concat_y_field, curr_y_field), dim='time')
        
    mean_x_field = (concat_x_field[xvec_attr]).sum(dim=['time']) / len(ecco_ds_vectors)
    mean_x_field = mean_x_field.where(mean_x_field != 0)
    mean_y_field = (concat_y_field[yvec_attr]).sum(dim=['time']) / len(ecco_ds_vectors)
    mean_y_field = mean_y_field.where(mean_y_field != 0)
    
    mean_fields[xvec_attr] = mean_x_field
    (mean_fields[xvec_attr]).data = (mean_fields[xvec_attr]).values
    mean_fields[yvec_attr] = mean_y_field
    (mean_fields[yvec_attr]).data = (mean_fields[yvec_attr]).values
    
    skip_k = True
    
    return mean_fields, skip_k

def comp_residuals(attributes, datasets, mean_dataset):
    
    """
    Computes residuals relative to a mean.
    
    attributes = variables to compute residuals of
    datasets = datasets to compute residuals for
    mean_dataset = mean of all attributes
    """
        
    residual_list = []
        
    for dataset in datasets:
            
        residual = dataset.copy() * 0 
            
        concat_attr = xr.concat((dataset, -1*mean_dataset), dim='time')
        
        for attribute in attributes:
            residual_field = (concat_attr[attribute]).sum(dim=['time'])
            residual_field = residual_field.where(residual_field != 0)
            residual[attribute] = residual_field
            (residual[attribute]).data = (residual[attribute]).values
            
        residual_list.append(residual)
    
    return residual_list

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
        ecco_ds_scalar, skip_k_scalar = comp_temp_mean_scalar(1, ecco_ds_scalars, scalar_attr)
        scalar_mean = ecco_ds_scalar
        
    if len(ecco_ds_vectors) == 1:
        ecco_ds_vector = ecco_ds_vectors[0]
        
    elif len(ecco_ds_vectors) > 1:
        ecco_ds_vector, skip_k_vector = comp_temp_mean_vector(1, ecco_ds_vectors, xvec_attr, yvec_attr)
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