"""
Auxiliary visualization file for ECCO data.

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

from ecco_general import get_month_name, get_scalar_in_xy, ds_to_field, comp_temp_mean, ecco_resample, load_dataset
from ecco_field_variables import get_field_vars

from vorticity_functions import comp_local_Ro, comp_vorticity, comp_normal_strain, comp_shear_strain, comp_OkuboWeiss

import save_seasonal_avgs

##############################

plt.rcParams['font.size'] = 16
plt.rcParams['text.usetex'] = True

##############################

def cbar_label(scalar_attr):
    
    """
    Returns label for plot colorbar.
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

##############################

def pcolormesh_quiver_title(ecco_ds_grid, k_plot, datestr, scalar_attr, xvec_attr, resid=False):
    
    """
    Returns title for contourf-quiver plot.
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
                    'PHIHYDcR': 'Hydrostatic pressure anomaly'}
    variable_name = variable_dict[variable]
    
    title = variable_name + ' in BGR at {}, {} \n'.format(depthstr, datestr)
    
    return title

##############################

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
    
##############################
    
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

##############################

def plot_Ro_l(Ro_l_list, zeta_field, lon_centers, lat_centers, seasonal, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, scalar_bounds=[1e-4, 1e-2], monthstr=None, seas_monthstr=None, seas_yearstr=None):

    """
    Computes and plots local Rossby number corresponding to monthly vorticity field.
    """

    Ro_l = comp_local_Ro(zeta_field, lat_centers) #Compute Ro_l
            
    #Define output file name
        
    if not seasonal:
        Ro_l_outfile = join(outdir, 'monthly', 'localRo_k{}_{}{}.pdf'.format(str(k), monthstr, yearstr))
    elif seasonal:
        Ro_l_outfile = join(outdir, 'seasonal', 'localRo_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr))

    #Plot Ro_l
    ArcCir_pcolormesh(ds_grid, [Ro_l], resolution, 'Reds', lon_centers, lat_centers, None, datestr, 'Ro_l', scalar_bounds=scalar_bounds, k_plot=k, extend='both', logscale=True, outfile=Ro_l_outfile, lats_lons=lats_lons)

    #Save Ro_l data and return it

    Ro_l_list.append(Ro_l)
    return Ro_l_list
    
##############################

def plot_OW(OW_list, zeta_field, seasonal, yearstr, year, outdir, k, datdirname, ds_grid, lon_centers, lat_centers, latmin, latmax, lonmin, lonmax, resolution, datestr, lats_lons, monthstr=None, datdir=None, season_start=None, season_end=None, endyearstr=None, seas_monthstr=None, seas_yearstr=None, seasonaldatdir=None, scalar_bounds=[-1e-14, 1e-14]):
    
    """
    Computes and plots Okubo-Weiss parameter corresponding to monthly velocity and vorticity profiles.
    """
    
    if not seasonal:
        
        #Get monthly velocity data
        
        vel_monthly_shortname, vel_monthly_nc_str = get_field_vars('UVELVVEL')
        vel_file = join(datdir, vel_monthly_shortname, vel_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        ds_vel = load_dataset(vel_file)
        (ds_vel['UVEL']).data, (ds_vel['VVEL']).data = (ds_vel['UVEL']).values, (ds_vel['VVEL']).values
        
        #Define output file name
        OW_outfile = join(outdir, 'monthly', 'OW_k{}_{}{}.pdf'.format(str(k), monthstr, yearstr))
        
    elif seasonal:
        
        #Get seasonal velocity data
        
        vel_seas_file = join(seasonaldatdir, "avg_UVELVVEL_"+season_start+yearstr+"-"+season_end+endyearstr+".nc")
        ds_vel = xr.open_mfdataset(vel_seas_file, engine="scipy")
        ds_vel.load()
        
        #Define output file name
        OW_outfile = join(outdir, 'seasonal', 'OW_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr))
        
    #Compute strain terms
    
    xgcm_grid = ecco.get_llc_grid(ds_grid)
    normal_strain = comp_normal_strain(xgcm_grid, ds_vel['UVEL'], ds_vel['VVEL'], ds_grid.dxG, ds_grid.dyG, ds_grid.rA).isel(k=k).squeeze()
    shear_strain = comp_shear_strain(xgcm_grid, ds_vel['UVEL'], ds_vel['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz).isel(k=k).squeeze()
    normal_strain= ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    shear_strain = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    
    OW = comp_OkuboWeiss(zeta_field, normal_strain, shear_strain) #Compute OW 
    
    #Plot OW
    ArcCir_pcolormesh(ds_grid, [OW], resolution, 'seismic', lon_centers, lat_centers, None, datestr, 'OW', scalar_bounds=scalar_bounds, k_plot=k, extend='both', outfile=OW_outfile, lats_lons=lats_lons)

    OW_list.append(OW)
    return OW_list

##############################

def plot_pcolormesh_k_plane(ds_grid, ds_scalar_list, k, scalar_attr, latmin, latmax, lonmin, lonmax, resolution, cmap, datestr, vmin, vmax, outfile, lats_lons, datdir, year, Ro_l_list, OW_list, yearstr, outdir=None, monthstr=None, seas_monthstr=None, seas_yearstr=None, logscale=True, seasonal=False, multiple_seas=False, annual=False, season_start=None, season_end=None, endyearstr=None, season_years=None, years=None, startyr=None, datdirname=None, seasonaldatdir=None, data_seasons=None, lon_centers=None, lat_centers=None):
    
    """
    Creates pcolormesh plot on plane of constant k.
    """
    
    #Take temporal average, if needed
 
    ds_scalar_mean = comp_temp_mean(ds_scalar_list)
    
    if type(ds_scalar_list[0]) == xr.Dataset: 
        
        ds_scalar_mean = ds_scalar_mean.isel(k=k) #Isolate k-plane
    
        #Convert scalar DataSet to useful field
        lon_centers, lat_centers, lon_edges, lat_edges, scalar = ds_to_field(ds_grid, ds_scalar_mean, scalar_attr, latmin, latmax, lonmin, lonmax, resolution)
        
    else: 
        scalar = ds_scalar_mean #lat/lon_centers are to be input as kwargs in this case
    
    if seasonal:
        seas_yearstr = yearstr
        datestr = '{}, {}'.format(seas_monthstr, seas_yearstr)
    
    #Plot scalar data
    ArcCir_pcolormesh(ds_grid, [scalar], resolution, cmap, lon_centers, lat_centers, None, datestr, scalar_attr, scalar_bounds=[vmin, vmax], k_plot=k, extend='both', outfile=outfile, lats_lons=lats_lons)    
  
    if scalar_attr == 'ZETA': #If vorticity, also compute and plot Ro_l, OW

        if annual: #If plotting an annual average

            #Plot Ro_l
            ArcCir_pcolormesh(ds_grid, Ro_l_list, resolution, 'Reds', lon_centers, lat_centers, None, yearstr, 'Ro_l', scalar_bounds=[1e-4, 1e-2], k_plot=k, extend='both', logscale=True, outfile=join(outdir, 'yearly', 'localRo_k{}_{}.pdf'.format(str(k), yearstr)), lats_lons=lats_lons)

            #Plot OW
            ArcCir_pcolormesh(ds_grid, OW_list, resolution, 'seismic', lon_centers, lat_centers, None, yearstr, 'OW', scalar_bounds=[-0.1e-13, 0.1e-13], k_plot=k, extend='both', outfile=join(outdir, 'yearly', 'OW_k{}_{}.pdf'.format(str(k), yearstr)), lats_lons=lats_lons)
        
        elif multiple_seas: #If plotting interannual seasonal average
            
            seas_yearstr = str(startyr) + "-" + str(startyr + (years-1) + season_years[-1]) #For titles
     
            ArcCir_pcolormesh(ds_grid, Ro_l_list, resolution, 'Reds', lon_centers, lat_centers, None, '{}, {}'.format(seas_monthstr, seas_yearstr), 'Ro_l', scalar_bounds=[1e-4, 1e-2], k_plot=k, extend='both', logscale=True, outfile=join(outdir, 'interannual', 'localRo_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)

            ArcCir_pcolormesh(ds_grid, OW_list, resolution, 'seismic', lon_centers, lat_centers, None, '{}, {}'.format(seas_monthstr, seas_yearstr), 'OW', scalar_bounds=[-0.1e-13, 0.1e-13], k_plot=k, extend='both', outfile=join(outdir, 'interannual', 'OW_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)
        
        elif not seasonal:
                
            #Compute and plot local Rossby number for the month
            Ro_l_list = plot_Ro_l(Ro_l_list, scalar, lon_centers, lat_centers, False, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, monthstr=monthstr)
                
            #Compute and plot OW for the month
            OW_list = plot_OW(OW_list, scalar, False, yearstr, year, outdir, k, datdirname, ds_grid, lon_centers, lat_centers, latmin, latmax, lonmin, lonmax, resolution, datestr, lats_lons, monthstr=monthstr, datdir=datdir)
          
            return Ro_l_list, OW_list
            
        elif seasonal:    
            
            #Compute and plot local Rossby number for the season
            Ro_l_list = plot_Ro_l(Ro_l_list, scalar, lon_centers, lat_centers, True, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, seas_monthstr=seas_monthstr, seas_yearstr=seas_yearstr)
                
            #Make sure the seasonal velocity file exists already
                
            vel_seas_file = join(seasonaldatdir, "avg_UVELVVEL_"+season_start+yearstr+"-"+season_end+endyearstr+".nc") #Define filename
                
            if not os.path.exists(vel_seas_file): #If it doesn't exist, compute it
                save_seasonal_avgs.main(field='UVELVVEL', years=[year], start_month=season_start, end_month=season_end, usecompdata=False, datdir=datdirname, outdir=seasonaldatdir)
                
            #Compute and plot OW for the season
            OW_list = plot_OW(OW_list, scalar, True, yearstr, year, outdir, k, datdirname, ds_grid, lon_centers, lat_centers, latmin, latmax, lonmin, lonmax, resolution, datestr, lats_lons, season_start=season_start, season_end=season_end, endyearstr=endyearstr, seas_monthstr=seas_monthstr, seas_yearstr=seas_yearstr, seasonaldatdir=seasonaldatdir)
                
            data_seasons.append(scalar)
    
    if seasonal:
        return Ro_l_list, OW_list, data_seasons, lon_centers, lat_centers
            
##############################

def plot_pcm_quiver_k_plane(ds_grid, ds_scalar, k, scalar_attr, latmin, latmax, lonmin, lonmax, resolution, vector_dir, vector_monthly_nc_str, yearstr, monthstr, year, xvec_attr, yvec_attr, datdirname, compdatdir, lats_lons, vectorECCO, Delta_u_outfile=None):
    
    #Convert scalar DataSet to useful field
    lon_centers, lat_centers, lon_edges, lat_edges, scalar = ds_to_field(ds_grid, ds_scalar.isel(k=k), scalar_attr, latmin, latmax, lonmin, lonmax, resolution)

    if vectorECCO:
                            
        ds_vector = load_ECCO_dataset.main(variable_dir=vector_dir, variable_monthly_nc_str=vector_monthly_nc_str, yearstr=yearstr, monthstr=monthstr, year=year, scalar_attr=None, xvec_attr=xvec_attr, datdir=datdirname)#config['datdir'])
                            
        #Interpolate and rotate vector
                            
        (ds_vector[xvec_attr]).data, (ds_vector[yvec_attr]).data = (ds_vector[xvec_attr]).values, (ds_vector[yvec_attr]).values  
        vecE, vecN = rotate_vector(ds_grid, ds_vector, xvec_attr, yvec_attr)
        vecE, vecN = vecE.isel(k=k).squeeze(), vecN.isel(k=k).squeeze()
                            
    elif not vectorECCO:
                            
        curr_vector_file = join(vector_dir, vector_monthly_nc_str+yearstr+"-"+monthstr+".nc")

        if os.path.exists(curr_vector_file): #Look for the file
            ds_vector = xr.open_mfdataset(curr_vector_file, engine="scipy") #Load monthly vector file into workspace

        else: #If it doesn't exist, compute it
            
            compute_monthly_avgs.main(latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax, startyr=year, years=1, datdir=datdirname, outdir=compdatdir)#config['datdir'], outdir=compdatdir)
            ds_vector = load_dataset(curr_vector_file)
            vecE, vecN = rotate_u_g(ds_grid, ds_vector[xvec_attr], ds_vector[yvec_attr], k)
                            
            if xvec_attr == 'UG': #If u_g, also plot geostrophy metric
                                
                vel_monthly_shortname, vel_monthly_nc_str = get_field_vars('UVELVVEL')
 
                ds_vel = load_ECCO_dataset.main(variable_dir=join(datdir, vel_monthly_shortname), variable_monthly_nc_str=vel_monthly_nc_str, yearstr=yearstr, monthstr=monthstr, year=year, scalar_attr='UVEL', xvec_attr=None, datdir=datdirname)#config['datdir'])
                                
                #Interpolate velocities to centres
                (ds_vel['UVEL']).data, (ds_vel['VVEL']).data = (ds_vel['UVEL']).values, (ds_vel['VVEL']).values
                velocity_interp = get_vector_in_xy(ds_grid, ds_vel, 'UVEL', 'VVEL') 
                u, v = velocity_interp['X'].isel(k=k), velocity_interp['Y'].isel(k=k)

                #Compute our geostrophy metric
                Delta_u = comp_geos_metric(u.squeeze(), v.squeeze(), vecE, vecN)
                                
                lon_centers, lat_centers, lon_edges, lat_edges, Delta_u_plot = ecco_resample(ds_grid, Delta_u, latmin, latmax, lonmin, lonmax, resolution)
                                
                #ArcCir_pcolormesh(ds_grid, k, [Delta_u_plot], resolution, 'Reds', lon_centers, lat_centers, None, monthstr+"-"+yearstr, 'Delta_u', scalar_bounds=[0, 1], k_plot=k, extend='max', outfile=join(outdir, 'monthly', '{}_k{}_{}{}.pdf'.format('Delta_u', str(k), monthstr, yearstr)), lats_lons=lats_lons) 
                ArcCir_pcolormesh(ds_grid, k, [Delta_u_plot], resolution, 'Reds', lon_centers, lat_centers, None, monthstr+"-"+yearstr, 'Delta_u', scalar_bounds=[0, 1], k_plot=k, extend='max', outfile=Delta_u_outfile, lats_lons=lats_lons)
    
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
##############################