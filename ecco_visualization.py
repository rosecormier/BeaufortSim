import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py as ecco

from matplotlib.gridspec import GridSpec

plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True

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
    
    fig = plt.figure(figsize=(6, 8), dpi=90)
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    #curr_obj = ecco.plot_proj_to_latlon_grid(ds_grid.XC, ds_grid.YC, tmp_plot, projection_type='stereo', plot_type='contourf', cmin=cmin, cmax=cmax, cmap=cmap, show_colorbar=True, lat_lim=65,zorder=50)
    cs1 = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, field_nearest_quarter_deg, levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)), no_levels), transform=ccrs.PlateCarree(), extend='both', cmap='viridis')
    cs2 = ax.contour(new_grid_lon_centers, new_grid_lat_centers, field_nearest_quarter_deg, colors='r', alpha=0.8, linewidths=0.5, zorder=100, transform=ccrs.PlateCarree(), levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)), no_levels))
    
    #ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    ax.gridlines()
    
    cbar = fig.colorbar(cs1, ticks=range(int(np.floor(vmin)), int(np.ceil(vmax)), 1))
    
    plt.savefig(vis_dir + filename + '.pdf')
    plt.close()
    
def ArcCir_contourf_quiver(ecco_ds_grid, k_plot, ecco_ds_scalar, ecco_ds_vector, scalar_attr, xvec_attr, yvec_attr, resolution, cmap, monthstr, yearstr, outfile="", vmin=0, vmax=0, latmin=70.0, latmax=85.0, lonmin=-180.0, lonmax=-90.0, no_levels=30, scale_factor=1, arrow_spacing=10, quiv_scale=1):
    
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
    monthstr = month to plot
    yearstr = year to plot
    outfile = output file name
    vmin = minimum of scalar data
    vmax = maximum of scalar data
    lat/lonmin/max = latitude/longitude bounds
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
    velc = XGCM_grid.interp_2d_vector({'X': (ecco_ds_vector[xvec_attr]).isel(k=k_plot), 'Y': (ecco_ds_vector[yvec_attr]).isel(k=k_plot)}, boundary='fill')
    
    velE = velc['X'] * ds_grid['CS'] - velc['Y'] * ds_grid['SN']
    velN = velc['X'] * ds_grid['SN'] + velc['Y'] * ds_grid['CS']
    
    tmp_plot = ds_grid[scalar_attr].squeeze()
    
    new_grid_delta_lat, new_grid_delta_lon = resolution, resolution
    
    new_grid_min_lat = latmin
    new_grid_max_lat = latmax

    new_grid_min_lon = lonmin
    new_grid_max_lon = lonmax

    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, field_nearest = ecco.resample_to_latlon(ds_grid.XC, \
                                ds_grid.YC, tmp_plot, new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, \
                                new_grid_delta_lon, fill_value = np.NaN, mapping_method = 'nearest_neighbor', radius_of_influence = 120000)
    
    field_copy = field_nearest.copy()
    field_copy = field_copy[~np.isnan(field_copy)]
    
    if vmin == vmax:
        vmax = scale_factor * np.max(field_copy)
        vmin = scale_factor * np.min(field_copy)
    
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
    
    cs1 = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, field_nearest, levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)), no_levels), transform=ccrs.PlateCarree(), extend='both', cmap='viridis')
    cs2 = ax.contour(new_grid_lon_centers, new_grid_lat_centers, field_nearest, colors='r', alpha=0.8, linewidths=0.5, zorder=100, transform=ccrs.PlateCarree(), levels=np.linspace(int(np.floor(vmin)), int(np.ceil(vmax)), no_levels))
    
    u_plot = (velE).squeeze()
    v_plot = (velN).squeeze()
    
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, u_nearest = ecco.resample_to_latlon(ds_grid.XC, \
                                ds_grid.YC, u_plot, new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, \
                                new_grid_delta_lon, fill_value = np.NaN, mapping_method = 'nearest_neighbor', radius_of_influence = 120000)
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, v_nearest = ecco.resample_to_latlon(ds_grid.XC, \
                                ds_grid.YC, v_plot, new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, \
                                new_grid_delta_lon, fill_value = np.NaN, mapping_method = 'nearest_neighbor', radius_of_influence = 120000)
    
    quiv = ax.quiver(new_grid_lon_centers, new_grid_lat_centers, u_nearest, v_nearest, color='k', transform=ccrs.PlateCarree(), scale=1, regrid_shape=60, zorder=150)
    
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    ax.gridlines()
    
    if k_plot == 0:
        depthstr = 'ocean surface'
        
    elif k_plot != 0:
        depth = - (ecco_ds_scalar[scalar_attr]).Z[k_plot].values
        depthstr = str(depth) + ' m depth'
       
    ax.set_title('Pressure anomaly and water velocity in Arctic Circle \n at {} ({}-{})'.format(depthstr, yearstr, monthstr))
        
    cbar = fig.colorbar(cs1, ticks=range(int(np.floor(vmin)), int(np.ceil(vmax)), 1), label=r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$')
    
    plt.savefig(outfile)
    plt.close()
    
def ArcCir_contourf_quiver_grid(ecco_ds_grid, k_plot, ecco_ds_scalars, ecco_ds_vectors, scalar_attr, vmin, vmax, xvec_attr, yvec_attr, resolution, cmap, monthstrs, yearstrs, outfile="", latmin=70.0, latmax=85.0, lonmin=-180.0, lonmax=-90.0, no_levels=30, scale_factor=1, arrow_spacing=10, quiv_scale=0.5, nrows=3, ncols=4):
    
    """
    ecco_ds_grid = ECCO grid
    k_plot = depth index to plot at
    ecco_ds_scalars = scalar DataSets
    ecco_ds_vectors = vector DataSets
    scalar_attr = string corresponding to scalar attribute to plot
    vmin = minimum of scalar attribute
    vmax = maximum of scalar attribute
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
    """
    
    plt.rcParams['font.size'] = 40
    
    monthnames = {"01": "January", "02": "February", "03": "March", "04": "April", "05": "May", "06": "June",
                 "07": "July", "08": "August", "09": "September", "10": "October", "11": "November", "12": "December"}
    
    if k_plot == 0:
        depthstr = 'ocean surface'

    elif k_plot != 0:
        depth = - (ecco_ds_scalars[0][scalar_attr]).Z[k_plot].values
        depthstr = str(depth) + ' m depth'
    
    mainfig = plt.figure(figsize=(44, 36))
    
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
        
        ecco_ds_scalar_k = ecco_ds_scalar.isel(k=k_plot)

        ds_grid = ecco_ds_grid.copy()
        ds_grid[scalar_attr] = ecco_ds_scalar_k[scalar_attr]
        ds_grid = ds_grid.load()

        XGCM_grid = ecco.get_llc_grid(ds_grid)
        velc = XGCM_grid.interp_2d_vector({'X': (ecco_ds_vector[xvec_attr]).isel(k=k_plot), 'Y': (ecco_ds_vector[yvec_attr]).isel(k=k_plot)}, boundary='fill')

        velE = velc['X'] * ds_grid['CS'] - velc['Y'] * ds_grid['SN']
        velN = velc['X'] * ds_grid['SN'] + velc['Y'] * ds_grid['CS']

        tmp_plot = ds_grid[scalar_attr].squeeze()

        new_grid_delta_lat, new_grid_delta_lon = resolution, resolution

        new_grid_min_lat = latmin
        new_grid_max_lat = latmax

        new_grid_min_lon = lonmin
        new_grid_max_lon = lonmax

        new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, field_nearest = ecco.resample_to_latlon(ds_grid.XC, \
                                    ds_grid.YC, tmp_plot, new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, \
                                    new_grid_delta_lon, fill_value = np.NaN, mapping_method = 'nearest_neighbor', radius_of_influence = 120000)

        ax = mainfig.add_subplot(nrows, ncols, i + 1, projection=ccrs.NorthPolarStereo())
        ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
        
        cs1 = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, field_nearest, levels=np.linspace(vmin, vmax, no_levels), transform=ccrs.PlateCarree(), extend='both', cmap='viridis')
        cs2 = ax.contour(new_grid_lon_centers, new_grid_lat_centers, field_nearest, colors='r', alpha=0.8, linewidths=0.5, zorder=100, transform=ccrs.PlateCarree(), levels=np.linspace(vmin, vmax, no_levels))

        u_plot = (velE).squeeze()
        v_plot = (velN).squeeze()

        new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, u_nearest = ecco.resample_to_latlon(ds_grid.XC, \
                                    ds_grid.YC, u_plot, new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, \
                                    new_grid_delta_lon, fill_value = np.NaN, mapping_method = 'nearest_neighbor', radius_of_influence = 120000)
        new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, v_nearest = ecco.resample_to_latlon(ds_grid.XC, \
                                    ds_grid.YC, v_plot, new_grid_min_lat, new_grid_max_lat, new_grid_delta_lat, new_grid_min_lon, new_grid_max_lon, \
                                    new_grid_delta_lon, fill_value = np.NaN, mapping_method = 'nearest_neighbor', radius_of_influence = 120000)

        quiv = ax.quiver(new_grid_lon_centers, new_grid_lat_centers, u_nearest, v_nearest, color='k', transform=ccrs.PlateCarree(), scale=1, regrid_shape=30, zorder=150)

        ax.add_feature(cfeature.LAND)
        ax.coastlines()
        ax.gridlines()
        
        ax.set_title('{} {}'.format(monthnames[monthstr], yearstr))

    mainfig.suptitle('Pressure anomaly and water velocity in Arctic Circle at {}'.format(depthstr), size=80)
    mainfig.savefig(outfile)
    plt.close()
    
    plt.rcdefaults()