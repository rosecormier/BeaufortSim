import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py as ecco

from matplotlib.gridspec import GridSpec

plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True

def get_scalar_in_xy(ecco_ds_grid, ecco_ds_scalar, scalar_attr, k_val=-1):
    
    """
    Loads scalar field in xy-grid.
    
    ecco_ds_grid = ECCO grid
    k_plot = depth index of interest
    ecco_ds_scalar = DataSet containing field
    scalar_attr = name of field
    """
    
    ds_grid = ecco_ds_grid.copy()

    if k_val < 0:
        
        ds_grid[scalar_attr] = ecco_ds_scalar[scalar_attr]
    
    elif k_val >= 0:
        
        ecco_ds_scalar_k = ecco_ds_scalar.isel(k=k_val)
        ds_grid[scalar_attr] = ecco_ds_scalar_k[scalar_attr]
        
    ds_grid = ds_grid.load()
    
    return ds_grid
    
def get_vector_in_xy(ecco_ds_grid, ecco_ds_vector, xvec_attr, yvec_attr, k_val=-1):
    
    """
    Loads vector field in xy-grid.
    
    ecco_ds_grid = ECCO grid
    k_plot = depth index of interest
    ecco_ds_vector = DataSet containing vector field
    xvec_attr = name of x-comp of vector field
    yvec_attr = name of y-comp of vector field
    """

    ds_grid = ecco_ds_grid.copy()
    ds_grid = ds_grid.load()
    XGCM_grid = ecco.get_llc_grid(ds_grid)
    
    if k_val < 0:
        velc = XGCM_grid.interp_2d_vector({'X': ecco_ds_vector[xvec_attr], \
                                           'Y': ecco_ds_vector[yvec_attr]}, boundary='fill')
    
    elif k_val >= 0:
        velc = XGCM_grid.interp_2d_vector({'X': (ecco_ds_vector[xvec_attr]).isel(k=k_val), \
                                       'Y': (ecco_ds_vector[yvec_attr]).isel(k=k_val)}, \
                                      boundary='fill')
    
    return velc
        
def ArcCir_contourf(k_plot, ecco_ds, attribute, ecco_ds_grid, resolution, cmap, no_levels, \
                    vis_dir, filename, scale_factor=1):
    
    """
    Creates contourf plot of scalar variable in a subdomain of the Arctic.
    
    k_plot = depth index to plot at
    ecco_ds = DataSet
    attribute = string corresponding to attribute to plot
    ecco_ds_grid = ECCO grid
    resolution = resolution (both lat and lon) in degrees
    cmap = colormap name
    no_levels = number of contour levels
    vis_dir = visualization directory
    filename = output file name
    scale_factor = colorbar multiplier
    """
    
    ds_grid = get_scalar_in_xy(ecco_ds_grid, ecco_ds, attribute, k_val=k_plot)
    
    field = ds_grid[attribute].squeeze()
    
    new_grid_delta_lat, new_grid_delta_lon = resolution, resolution
    new_grid_min_lat, new_grid_max_lat = -90, 90
    new_grid_min_lon, new_grid_max_lon = -180, 180

    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, \
    new_grid_lat_edges, field_nearest = \
    ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, field, new_grid_min_lat, \
                            new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, \
                            new_grid_max_lon, new_grid_delta_lon, fill_value = np.NaN, \
                            mapping_method = 'nearest_neighbor', \
                            radius_of_influence = 120000)
    
    field_copy = field.isel(tile=6).squeeze().values.copy()
    field_copy = field_copy[~np.isnan(field_copy)]
    
    vmax, vmin = scale_factor * np.max(field_copy), scale_factor * np.min(field_copy)
    
    fig = plt.figure(figsize=(6, 8), dpi=90)
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    
    cs1 = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, field_nearest, \
                      levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)), no_levels), \
                      transform=ccrs.PlateCarree(), extend='both', cmap=cmap)
    cs2 = ax.contour(new_grid_lon_centers, new_grid_lat_centers, field_nearest, \
                     colors='r', alpha=0.8, linewidths=0.5, zorder=100, \
                     transform=ccrs.PlateCarree(), levels=np.linspace(int(np.floor(vmin)), \
                                                                      int(np.ceil(vmax)), \
                                                                      no_levels))
    
    ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    ax.gridlines()
    
    cbar = fig.colorbar(cs1, ticks=range(int(np.floor(vmin)), int(np.ceil(vmax)), 1))
    
    plt.savefig(vis_dir + filename + '.pdf')
    plt.close()

def comp_temporal_mean(ecco_ds_grid, k_val, ecco_ds_scalars, scalar_attr):
    
    """
    Computes temporal mean of a scalar field.
    
    ecco_ds_grid = ECCO grid
    k_val = depth value of interest
    ecco_ds_scalars = scalar DataSets
    scalar_attr = string corresponding to scalar attribute of interest
    """ 

    if k_val == 0:
            depthstr = 'ocean surface'
    
    elif k_val != 0:
        depth = - ((ecco_ds_scalars[0])[scalar_attr]).Z[k_val].values
        depthstr = str(depth) + ' m depth'
    
    concat_field = ((ecco_ds_scalars[0]).copy()).isel(k=k_val)* 0
    
    for dataset in ecco_ds_scalars:
        
        curr_field = (dataset.isel(k=k_val))
        concat_field = xr.concat((concat_field, curr_field), dim='time')

    mean_field = (concat_field[scalar_attr]).sum(dim=['time']) / len(ecco_ds_scalars)
    
    return mean_field, depthstr

def ArcCir_contourf_quiver(ecco_ds_grid, k_plot, ecco_ds_scalar, ecco_ds_vector, scalar_attr, \
                           xvec_attr, yvec_attr, resolution, cmap, scalar_bounds, \
                           outfile="", latmin=70.0, latmax=85.0, lonmin=-180.0, lonmax=-90.0, \
                           no_levels=30, scale_factor=1, arrow_spacing=10, quiv_scale=1, title="", \
                           zplane=False, zplanedepth=0):
    
    """
    Creates contourf plot of scalar variable in a subdomain of the Arctic,
        overlaid with quiver plot of vector variable.
    
    ecco_ds_grid = ECCO grid
    k_plot = depth index to plot at
    ecco_ds_scalar = scalar DataSet
    ecco_ds_vector = vector DataSet
    scalar_attr = string corresponding to scalar attribute to plot
    xvec_attr = string corresponding to x-comp of vector to plot
    yvec_attr = string corresponding to y-comp of vector to plot
    resolution = resolution (both lat and lon) in degrees
    cmap = colormap name
    outfile = output file name
    vmin = minimum of scalar data
    vmax = maximum of scalar data
    lat/lonmin/max = latitude/longitude bounds
    no_levels = number of contour levels
    scale_factor = colorbar multiplier
    arrow_spacing = quiver arrow spacing in gridpoints
    quiv_scale = quiver plot scale
    title = plot title
    """
    
    if zplane:
        k_val = -1
        depthstr = zplanedepth
        
    else:
        k_val = k_plot
    
    ds_grid = get_scalar_in_xy(ecco_ds_grid, ecco_ds_scalar, scalar_attr, k_val=k_val)

    velc = get_vector_in_xy(ecco_ds_grid, ecco_ds_vector, xvec_attr, yvec_attr, k_val=k_val)
    velE = velc['X'] * ds_grid['CS'] - velc['Y'] * ds_grid['SN']
    velN = velc['X'] * ds_grid['SN'] + velc['Y'] * ds_grid['CS']
    
    curr_field = (ds_grid[scalar_attr]).squeeze()
    
    new_grid_delta_lat, new_grid_delta_lon = resolution, resolution
    new_grid_min_lat, new_grid_max_lat = latmin, latmax
    new_grid_min_lon, new_grid_max_lon = lonmin, lonmax

    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, \
    field_nearest = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, curr_field, new_grid_min_lat, \
                                            new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, \
                                            new_grid_max_lon, new_grid_delta_lon, fill_value = np.NaN, \
                                            mapping_method = 'nearest_neighbor', \
                                            radius_of_influence = 120000)
    
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
    ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, u_plot, new_grid_min_lat, new_grid_max_lat, \
                            new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon, \
                            fill_value = np.NaN, mapping_method = 'nearest_neighbor', \
                            radius_of_influence = 120000)
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, v_nearest = \
    ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, v_plot, new_grid_min_lat, new_grid_max_lat, \
                            new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon, \
                            fill_value = np.NaN, mapping_method = 'nearest_neighbor', \
                            radius_of_influence = 120000)
    
    quiv = ax.quiver(new_grid_lon_centers, new_grid_lat_centers, u_nearest, v_nearest, color='k', \
                     transform=ccrs.PlateCarree(), scale=1, regrid_shape=60, zorder=150)
    
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    ax.gridlines()
    
    if not zplane:
        
        if k_plot == 0:
            depthstr = 'ocean surface'
        
        elif k_plot != 0:
            depth = - (ecco_ds_scalar[scalar_attr]).Z[k_plot].values
            depthstr = str(depth) + ' m depth'
    
    ax.set_title(title)
        
    cbar = fig.colorbar(cs1, ticks=range(int(np.floor(vmin)), int(np.ceil(vmax)), 1), \
                        label=r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$')
    
    plt.savefig(outfile)
    plt.close()
    
def ArcCir_contourf_quiver_grid(ecco_ds_grid, k_plot, ecco_ds_scalars, ecco_ds_vectors, scalar_attr, \
                                scalar_bounds, xvec_attr, yvec_attr, resolution, \
                                cmap, monthstrs, yearstrs, outfile="", latmin=70.0, latmax=85.0, \
                                lonmin=-180.0, lonmax=-90.0, no_levels=15, scale_factor=1, \
                                arrow_spacing=10, quiv_scale=0.3, nrows=3, ncols=4, title=""):
    
    """
    Creates array of contourf plots of scalar variable in a subdomain of the Arctic,
        overlaid with quiver plots of vector variable.
    
    ecco_ds_grid = ECCO grid
    k_plot = depth index to plot at
    ecco_ds_scalars = scalar DataSets
    ecco_ds_vectors = vector DataSets
    scalar_attr = string corresponding to scalar attribute to plot
    scalar_bounds = bounds on scalar attribute
    xvec_attr = string corresponding to x-comp of vector to plot
    yvec_attr = string corresponding to y-comp of vector to plot
    resolution = resolution (both lat and lon) in degrees
    cmap = colormap name
    monthstrs = strings of months to plot
    yearstrs = strings of years to plot
    outfile = output file name
    lat/lonmin/max = latitude/longitude bounds
    no_levels = number of contour levels
    scale_factor = colorbar multiplier
    arrow_spacing = quiver arrow spacing in gridpoints
    quiv_scale = quiver plot scale
    nrows, ncols = number of rows and columns (resp.) in grid
    title = plot title
    """
    
    plt.rcParams['font.size'] = 40
    
    monthnames = {"01": "January", "02": "February", "03": "March", "04": "April", "05": "May", \
                  "06": "June", "07": "July", "08": "August", "09": "September", "10": "October", \
                  "11": "November", "12": "December"}
    
    if k_plot == 0:
        depthstr = 'ocean surface'

    elif k_plot != 0:
        depth = - (ecco_ds_scalars[0][scalar_attr]).Z[k_plot].values
        depthstr = str(depth) + ' m depth'
    
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
        
        ds_grid = get_scalar_in_xy(ecco_ds_grid, ecco_ds_scalar, scalar_attr, k_val=k_plot)
        
        velc = get_vector_in_xy(ecco_ds_grid, ecco_ds_vector, xvec_attr, yvec_attr, k_val=k_plot)
        velE = velc['X'] * ds_grid['CS'] - velc['Y'] * ds_grid['SN']
        velN = velc['X'] * ds_grid['SN'] + velc['Y'] * ds_grid['CS']

        curr_field = ds_grid[scalar_attr].squeeze()

        new_grid_delta_lat, new_grid_delta_lon = resolution, resolution
        new_grid_min_lat, new_grid_max_lat = latmin, latmax
        new_grid_min_lon, new_grid_max_lon = lonmin, lonmax

        new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, \
        field_nearest = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, curr_field, new_grid_min_lat, \
                                                new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, \
                                                new_grid_max_lon, new_grid_delta_lon, \
                                                fill_value = np.NaN, mapping_method = 'nearest_neighbor', \
                                                radius_of_influence = 120000)

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
        ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, u_plot, new_grid_min_lat, new_grid_max_lat, \
                                new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon, \
                                fill_value = np.NaN, mapping_method = 'nearest_neighbor', \
                                radius_of_influence = 120000)
        new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, v_nearest = \
        ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, v_plot, new_grid_min_lat, new_grid_max_lat, \
                                new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, new_grid_delta_lon, \
                                fill_value = np.NaN, mapping_method = 'nearest_neighbor', \
                                radius_of_influence = 120000)

        quiv = ax.quiver(new_grid_lon_centers, new_grid_lat_centers, u_nearest, v_nearest, color='k', \
                         transform=ccrs.PlateCarree(), scale=1, regrid_shape=30, zorder=150)

        ax.add_feature(cfeature.LAND)
        ax.coastlines()
        ax.gridlines()
        
        ax.set_title('\n {} {}'.format(monthnames[monthstr], yearstr))

    if title == "":
        title = "Monthly mean pressure anomaly and water velocity in Arctic Circle at {} \n".format(depthstr)
    
    mainfig.suptitle(title, size=80)
    mainfig.tight_layout()
    cbar = mainfig.colorbar(cs1, ax=mainfig.get_axes(), aspect=40, pad=0.05, ticks=range(vmin, vmax, 1), \
                            label=r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$', location='bottom')
    mainfig.savefig(outfile)
    plt.close()
    """
    if resid:
        
        absmax = 1.8
    
        for i in range(nplots):

            if i % ncols == 0:
                row += 1
                col = 0
            else:
                col += 1

            new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, resid = \
            ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, tmp_plots[i] - mean_plot, new_grid_min_lat, \
                                    new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, \
                                    new_grid_max_lon, new_grid_delta_lon, fill_value = np.NaN, \
                                    mapping_method = 'nearest_neighbor', radius_of_influence = 120000)

            ax = resid_fig.add_subplot(nrows, ncols, i + 1, projection=ccrs.NorthPolarStereo())
            ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())

            cs1 = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, resid, \
                              transform=ccrs.PlateCarree(), levels=np.linspace(-absmax, absmax, no_levels), \
                              extend='both', cmap='seismic')
            cs2 = ax.contour(new_grid_lon_centers, new_grid_lat_centers, resid, colors='r', alpha=0.8, \
                             linewidths=1.0, zorder=100, transform=ccrs.PlateCarree(), \
                             levels=np.linspace(-absmax, absmax, no_levels))

            ax.add_feature(cfeature.LAND)
            ax.coastlines()
            ax.gridlines()
             
            ax.set_title('\n {} {}'.format(monthnames[monthstrs[i]], yearstrs[i]))

        resid_title = "Monthly residuals of mean pressure anomaly and water velocity, relative to annual \
        mean, \n in Arctic Circle at {} \n".format(depthstr)
            
        resid_fig.suptitle(resid_title, size=80)
        resid_fig.tight_layout()
        cbar = mainfig.colorbar(cs1, ax=resid_fig.get_axes(), aspect=40, pad=0.05, \
                                ticks=range(int(np.ceil(-absmax)), int(np.floor(absmax))), \
                                label=r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$', location='bottom')
        resid_fig.savefig('test.png')
        plt.close()
    """
    
    plt.rcdefaults()