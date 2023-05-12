import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py as ecco

def cmap_zerocent_scale(plot, scale_factor):
    
    """
    Centre colormap at zero and scale relative to existing absolute maximum value, given plot object and scale_factor (float)
    Return new colormap limits
    """
    
    curr_clim = plot.get_clim()
    new_clim = (scale_factor * np.max(np.abs(curr_clim))) * np.array([-1, 1])
    plot.set_clim(new_clim)
    
    return new_clim
        
def ArcCir_contourf(k_plot, ecco_ds, attribute, ecco_ds_grid, resolution, cmap, no_levels, vis_dir, filename, scale_factor=1):
    
    """
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
    
    ecco_ds_k = ecco_ds.isel(k=k_plot)

    ds_grid = ecco_ds_grid.copy()
    ds_grid[attribute] = ecco_ds_k[attribute]
    ds_grid = ds_grid.load()

    tmp_plot = ds_grid[attribute].squeeze()
    
    new_grid_delta_lat, new_grid_delta_lon = resolution, resolution
    
    new_grid_min_lat = -90
    new_grid_max_lat = 90

    new_grid_min_lon = -180
    new_grid_max_lon = 180

    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, field_nearest_quarter_deg = ecco.resample_to_latlon(ds_grid.XC, \
                                ds_grid.YC, tmp_plot, new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, \
                                new_grid_delta_lon, fill_value = np.NaN, mapping_method = 'nearest_neighbor', radius_of_influence = 120000)
    
    field_copy = tmp_plot.isel(tile=6).squeeze().values.copy()
    field_copy = field_copy[~np.isnan(field_copy)]
    
    vmax = scale_factor * np.max(field_copy)
    vmin = scale_factor * np.min(field_copy)
    
    fig = plt.figure(figsize=(12, 6), dpi=90)
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    #curr_obj = ecco.plot_proj_to_latlon_grid(ds_grid.XC, ds_grid.YC, tmp_plot, projection_type='stereo', plot_type='contourf', cmin=cmin, cmax=cmax, cmap=cmap, show_colorbar=True, lat_lim=65,zorder=50)
    cs1 = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, field_nearest_quarter_deg, levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)),no_levels), transform=ccrs.PlateCarree(), extend='both', cmap='viridis')
    cs2 = ax.contour(new_grid_lon_centers, new_grid_lat_centers, field_nearest_quarter_deg, colors='r', alpha=0.8, linewidths=0.5, zorder=100, transform=ccrs.PlateCarree(), levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)),no_levels))
    
    ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    ax.gridlines()
    
    cbar = fig.colorbar(cs1, ticks=range(int(np.floor(vmin)), int(np.ceil(vmax)), 1))
    
    plt.savefig(vis_dir + filename + '.pdf')
    plt.close()
    
def ArcCir_contourf_quiver(ecco_ds_grid, k_plot, ecco_ds_scalar, ecco_ds_vector, scalar_attr, xvec_attr, yvec_attr, resolution, cmap, vis_dir, filename, no_levels=30, scale_factor=1, arrow_spacing=10, quiv_scale=1):
    
    """
    ecco_ds_grid = ECCO grid
    k_plot = depth index to plot at
    ecco_ds_scalar = scalar DataSet
    ecco_ds_vector = vector DataSet
    scalar_attr = string corresponding to scalar attribute to plot
    xvec_attr = string corresponding to x-comp of vector to plot
    yvec_attr = string corresponding to y-comp of vector to plot
    resolution = resolution (both lat and lon) in degrees
    cmap = colormap name
    vis_dir = visualization directory
    filename = output file name, no extension
    no_levels = number of contour levels
    scale_factor = colorbar multiplier
    arrow_spacing = quiver arrow spacing in gridpoints
    quiv_scale = quiver plot scale
    """
    
    ecco_ds_scalar_k = ecco_ds_scalar.isel(k=k_plot)

    ds_grid = ecco_ds_grid.copy()
    ds_grid[scalar_attr] = ecco_ds_scalar_k[scalar_attr]
    ds_grid = ds_grid.load()
    
    XGCM_grid = ecco.get_llc_grid(ds_grid)
    velc = XGCM_grid.interp_2d_vector({'X': ecco_ds_vector[xvec_attr].isel(k=k_plot), 'Y': ecco_ds_vector[yvec_attr].isel(k=k_plot)}, boundary='fill')
    
    velE = velc['X'] * ds_grid['CS'] - velc['Y'] * ds_grid['SN']
    velN = velc['X'] * ds_grid['SN'] + velc['Y'] * ds_grid['CS']
    
    tmp_plot = ds_grid[scalar_attr].squeeze()
    
    new_grid_delta_lat, new_grid_delta_lon = resolution, resolution
    
    new_grid_min_lat = -89
    new_grid_max_lat = 89

    new_grid_min_lon = -180
    new_grid_max_lon = 180

    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, field_nearest = ecco.resample_to_latlon(ds_grid.XC, \
                                ds_grid.YC, tmp_plot, new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, \
                                new_grid_delta_lon, fill_value = np.NaN, mapping_method = 'nearest_neighbor', radius_of_influence = 120000)
    
    field_copy = tmp_plot.isel(tile=6).squeeze().values.copy()
    field_copy = field_copy[~np.isnan(field_copy)]
    
    vmax = scale_factor * np.max(field_copy)
    vmin = scale_factor * np.min(field_copy)
    
    fig = plt.figure(figsize=(12, 6), dpi=90)
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    
    cs1 = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, field_nearest, levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)),no_levels), transform=ccrs.PlateCarree(), extend='both', cmap='viridis')
    cs2 = ax.contour(new_grid_lon_centers, new_grid_lat_centers, field_nearest, colors='r', alpha=0.8, linewidths=0.5, zorder=100, transform=ccrs.PlateCarree(), levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)),no_levels))
    
    u_plot = (velE).squeeze()
    v_plot = (velN).squeeze()
    
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, u_nearest = ecco.resample_to_latlon(ds_grid.XC, \
                                ds_grid.YC, u_plot, new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, \
                                new_grid_delta_lon, fill_value = np.NaN, mapping_method = 'nearest_neighbor', radius_of_influence = 120000)
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, v_nearest = ecco.resample_to_latlon(ds_grid.XC, \
                                ds_grid.YC, v_plot, new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, \
                                new_grid_delta_lon, fill_value = np.NaN, mapping_method = 'nearest_neighbor', radius_of_influence = 120000)
    
    #print(new_grid_lon_centers, new_grid_lat_centers, u_nearest, v_nearest)
    quiv = ax.quiver(new_grid_lon_centers, new_grid_lat_centers, u_nearest, v_nearest, color='k', transform=ccrs.PlateCarree(), scale=1, regrid_shape=1000)
    
    ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    ax.gridlines()
    
    cbar = fig.colorbar(cs1, ticks=range(int(np.floor(vmin)), int(np.ceil(vmax)), 1))
    
    plt.savefig(vis_dir + filename + '.pdf')
    plt.close()
