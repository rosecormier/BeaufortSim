"""
Auxiliary functions for visualizing ECCO data.

R. Cormier, 2024
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
from matplotlib.lines import Line2D
from os.path import join

import load_data_files

from functions_ecco_general import load_grid, scalar_to_grid, vector_to_grid
from functions_field_variables import get_field_variable, get_variable_str, \
get_cmap_and_symmetry, get_cbar_label, get_field_title

##############################

plt.rcParams['font.size'] = 16
plt.rcParams['text.usetex'] = True

##############################

def pcolormesh_quiver_title(ecco_ds_grid, k_plot, datestr, scalar_attr, 
                            xvec_attr, resid=False):
    """
    Returns title for contourf-quiver plot.
    """
    
    ds_grid = ecco_ds_grid.copy()
   
    depth = - ds_grid.Z[k_plot].values
    depthstr = str(depth) + ' m depth'
        
    scalar_dict = {'PHIHYDcR': 'Pressure anomaly', \
                  'OW': 'Okubo-Weiss parameter', \
                  'OW_geos': \
                       r'Okubo-Weiss parameter (computed from $\vec{u}_g$)', \
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
        title = '{} and {}'.format(scalar_str, vector_str) \
                + ' residuals (relative to annual mean) \n' \
                + 'in Arctic Circle at {}, {} \n'.format(depthstr, datestr)
    else:
        title = '{} and {}'.format(scalar_str, vector_str) \
                + ' in BGR \n at {}, {} \n'.format(depthstr, datestr)
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
                    'zeta_geos': r'Vorticity (computed from $\vec{u}_g$),' \
                                + r' normalized by $f_{mean}$', \
                    'OW_geos': r'Okubo-Weiss parameter' \
                                + r' (computed from $\vec{u}_g$)', \
                    'Ro_l': 'Local Rossby number', \
                    'geos_metric': r'Metric for geostrophy' \
        + r' $\frac{||\vec{u} - \vec{u}_g||}{|\vec{u}|| + ||\vec{u}_g||}$', \
                    'PHIHYDcR': 'Hydrostatic pressure anomaly', \
                    'DIVU': 'Divergence of horizontal water velocity', \
                    'DIVUEk': 'Divergence of Ekman current', \
                    'SIheff': 'Sea ice thickness'}
    variable_name = variable_dict[variable]
    
    title = variable_name + ' in Beaufort Gyre region' \ 
                          + ' at {}, {} \n'.format(depthstr, datestr)
    return title

##############################

def get_plot_title(scalar_field_name, vector_field_name, plot_plane_type, 
                   spatial_bounds, ds_grid, date_string):
    """
    Return main title for plot.
    """
    
    scalar_field_title, multiple_depths = get_field_title(scalar_field_name)

    if multiple_depths:
    
        if plot_plane_type == 'depth_index_const':
            k_val = int(spatial_bounds[0])
            depth = -ds_grid.Z[k_val].values
            depth_string = str(round(depth)) + ' m depth'
    
        if vector_field_name is None:
            title = '{} in BGR at {}, {} \n'.format(scalar_field_title, 
                                                    depth_string, date_string)

        elif vector_field_name is not None:
            vector_field_title, vector_multiple_depths = get_field_title(
                vector_field_name)
            if vector_multiple_depths:
                title = '{} and {} \n at {}, {} \n'.format(
                    scalar_field_title, vector_field_title, depth_string, 
                    date_string)
            elif not vector_multiple_depths:
                title = '{} at {} and {} in BGR, \n {} \n'.format(
                    scalar_field_title, depth_string, vector_field_title, 
                    date_string)
            
    elif not multiple_depths:
        
        if vector_field_name is None:
            title = '{}, {} \n'.format(scalar_field_title, date_string)

        elif vector_field_name is not None:
            vector_field_title, vector_multiple_depths = get_field_title(
                vector_field_name)
            
            if vector_multiple_depths:
                if plot_plane_type == 'depth_index_const':
                    k_val = int(spatial_bounds[0])
                    depth = -ds_grid.Z[k_val].values
                    depth_string = str(depth) + ' m depth'
                title = '{} and {} at {}, \n {} \n'.format(
                    scalar_field_title, vector_field_title, depth_string, 
                    date_string)
                
            elif not vector_multiple_depths:
                title = '{} \n and {}, {}'.format(scalar_field_title, 
                                                  vector_field_title, 
                                                  date_string)
        
    return title

##############################

def plot_geography(ax, labels=True, latmin=70.5, latmax=80.0, lonmin=-155.0,
                   lonmax=-120.0):
    """
    Add land, coastlines, and grid to an ax.
    """
    
    ax.add_feature(cfeature.LAND)
    ax.coastlines()
    lines = ax.gridlines(draw_labels = labels)
    lines2 = ax.gridlines(draw_labels = False, 
                          xlim = (lonmin, lonmax), 
                          ylim = (latmin, latmax), 
                          color = 'k')
    
    if labels:
        
        lines.xlocator = mticker.FixedLocator(np.arange(-180,-90,5))
        lines.xformatter = LONGITUDE_FORMATTER
        lines.xpadding = 10
        lines.xlabel_style = {'size': 12, 'color': 'grey', 'rotation': 0}
        lines.ylocator = mticker.FixedLocator(np.arange(70,86,2))
        lines.yformatter = LATITUDE_FORMATTER
        lines.ylabel_style = {'size': 12, 'color': 'grey', 'rotation': 0}
        lines.x_inline = False
        lines.y_inline = False
        lines.xlabels_top = None
        lines.ylabels_right = None
        
        lines2.xlocator = mticker.FixedLocator(np.arange(-155, -119.9, 35))
        lines2.ylocator = mticker.FixedLocator(np.arange(70.5, 80.1, 9.5))
        
        plt.draw()
        
    return ax

##############################

def get_scalar_bounds(scalar_field, scalar_bounds=None):
    """
    Set vmin, vmax for a scalar field.
    """
    
    if scalar_bounds == None: #Default to actual bounds on data
        vmin, vmax = np.nanmin(scalar_field), np.nanmax(scalar_field)
    else: #Case where scalar bounds are explicitly given
        vmin, vmax = scalar_bounds[0], scalar_bounds[1]
    return vmin, vmax

##############################

def get_pcolormesh(ax, lon_centers, lat_centers, scalar, field_name, vmin, vmax,
                   logscale=False):
    """
    Create pcolormesh object given an axis.
    """

    cmap, symmetry = get_cmap_and_symmetry(field_name)
    
    if symmetry: 
        #Reset bounds on colorbar to be symmetric
        abs_max_value = max(abs(vmax), abs(vmin))
        vmin, vmax = -abs_max_value, abs_max_value
        
    cmap = plt.get_cmap(cmap).copy()
    cmap.set_bad('grey')
        
    if logscale:
        color = ax.pcolormesh(lon_centers, lat_centers, scalar, 
                              transform = ccrs.PlateCarree(), 
                              cmap = cmap, 
                              norm = colors.LogNorm(vmin=vmin, vmax=vmax))
    elif not logscale:
        color = ax.pcolormesh(lon_centers, lat_centers, scalar, 
                              transform = ccrs.PlateCarree(), 
                              cmap = cmap, 
                              vmin = vmin, vmax = vmax)
    return ax, color

##############################

def get_quiver(ax, lon_centers, lat_centers, vec_E_comp, vec_N_comp, 
               quiv_scale=0.35):
    """
    Resample vector to lat-lon grid and get quiver object; add quiver to ax.
    """

    skip = (slice(0, -1, 2), slice(0, -1, 2))
    quiv = ax.quiver(lon_centers[skip], lat_centers[skip], 
                     vec_E_comp[skip], vec_N_comp[skip], 
                     color = 'k', 
                     transform = ccrs.PlateCarree(),
                     scale = quiv_scale, 
                     scale_units = 'width', 
                     regrid_shape = 30)
    return quiv

##############################

def ArcCir_pcolormesh(scalar_field_name, date_string, datdir_primary, 
                      datdir_secondary, time_ave_type, plot_plane_type, 
                      spatial_bounds, resolutions, outfile, extend='both', 
                      vector_field_name=None):
    """
    Create pcolormesh plot of a scalar variable in a subdomain of the Arctic.
    Optional - overlay quiver plot of a vector variable.
    """
    
    ds_grid = load_grid(datdir_primary)

    scalar_ds = load_data_files.main(field_name=scalar_field_name, 
                                     date_string=date_string, 
                                     datdir_primary=datdir_primary, 
                                     datdir_secondary=datdir_secondary, 
                                     time_ave_type=time_ave_type)

    if plot_plane_type == 'depth_index_const':
        depth = spatial_bounds[0]
        latmin, latmax = spatial_bounds[1], spatial_bounds[2] 
        lonmin, lonmax = spatial_bounds[3], spatial_bounds[4]
        lat_res, lon_res = resolutions[0], resolutions[1]
    #will add options 'latitude_const' and 'longitude_const'
        
    lons, lats, lon_edges, lat_edges, scalar_field = scalar_to_grid(ds_grid, 
                                        scalar_ds, 
                                        get_variable_str(scalar_field_name),
                                        depth, 
                                        latmin, latmax, lonmin, lonmax, 
                                        lat_res, lon_res)
    
    vmin, vmax = get_scalar_bounds(scalar_field) #Set scalar bounds

    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-135))    
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
    
    ax, color = get_pcolormesh(ax, lons, lats, 
                               scalar_field, scalar_field_name,
                               vmin, vmax) 
    
    if vector_field_name is not None:
        vector_ds = load_data_files.main(field_name=vector_field_name, 
                                         date_string=date_string, 
                                         datdir_primary=datdir_primary, 
                                         datdir_secondary=datdir_secondary, 
                                         time_ave_type=time_ave_type)
        vector_to_grid_output = vector_to_grid(ds_grid, vector_ds, 
                                               vector_field_name, depth, latmin,
                                               latmax, lonmin, lonmax, lat_res,
                                               lon_res)
        lons, lats = vector_to_grid_output[0], vector_to_grid_output[1]
        vec_E_comp = vector_to_grid_output[4]
        vec_N_comp = vector_to_grid_output[5]
        quiv = get_quiver(ax, lons, lats, vec_E_comp, vec_N_comp) 
   
    ax = plot_geography(ax)
    ax.set_title(get_plot_title(scalar_field_name, vector_field_name, 
                                plot_plane_type, spatial_bounds, ds_grid, 
                                date_string))
    
    fig.colorbar((color), ax = ax, 
                 label = get_cbar_label(scalar_field_name), 
                 extend = extend, 
                 location = 'bottom', 
                 shrink = 0.85, 
                 ticks = np.linspace(vmin, vmax, 5))
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()
    
##############################

def plot_Ro_l(Ro_l_list, zeta_field, lon_centers, lat_centers, seasonal, outdir, 
              k, yearstr, ds_grid, resolution, datestr, lats_lons, 
              scalar_bounds=[1e-4, 1e-2], monthstr=None, seas_monthstr=None, 
              seas_yearstr=None, quiver=False, vecE=None, vecN=None, 
              xvec_attr=None):
    """
    Computes and plots local Rossby number corresponding to monthly vorticity field.
    """

    Ro_l = comp_local_Ro(zeta_field, lat_centers)
            
    #Define output file name
    if not seasonal:
        Ro_l_outfile = join(outdir, 'monthly', 
                            'localRo_k{}_{}{}.pdf'.format(str(k), 
                                                          monthstr, yearstr))
    elif seasonal:
        Ro_l_outfile = join(outdir, 'seasonal', 
                            'localRo_k{}_{}_{}.pdf'.format(str(k), 
                                                           seas_monthstr, 
                                                           seas_yearstr))


    if not quiver:
        ArcCir_pcolormesh(ds_grid, [Ro_l], resolution, 'Reds', lon_centers, 
                          lat_centers, None, datestr, 'Ro_l', 
                          scalar_bounds=scalar_bounds, k_plot=k, extend='both',
                          logscale=True, outfile=Ro_l_outfile, 
                          lats_lons=lats_lons)
        
    elif quiver:  #need to fix outfile for this case
        ArcCir_pcolormesh_quiver(ds_grid, k, [Ro_l], [vecE], [vecN], resolution,
                                 'Reds', datestr, lon_centers, lat_centers, 
                                 scalar_attr='Ro_l', xvec_attr=xvec_attr, 
                                 scalar_bounds=scalar_bounds, logscale=True, 
                                 outfile=Ro_l_outfile, lats_lons=lats_lons)

    Ro_l_list.append(Ro_l)
    return Ro_l_list
    
##############################

def plot_OW(OW_list, zeta_field, lon_centers, lat_centers, seasonal, outdir, k,
            yearstr, ds_grid, resolution, datestr, lats_lons, monthstr=None, 
            datdir=None, season_start=None, season_end=None, endyearstr=None, 
            seas_monthstr=None, seas_yearstr=None, seasonaldatdir=None,
            scalar_bounds=[-1e-14, 1e-14], quiver=False, vecE=None, vecN=None, 
            xvec_attr=None):
    """
    Computes and plots Okubo-Weiss parameter 
    corresponding to monthly velocity and vorticity profiles.
    """
    
    if not seasonal:
        
        #Get monthly velocity data
        vel_monthly_shortname, vel_monthly_nc_str = get_field_vars('UVELVVEL')
        vel_file = join(datdir, vel_monthly_shortname, 
                        "{}{}-{}_ECCO_V4r4_native_llc0090.nc".format(
                                        vel_monthly_nc_str, yearstr, monthstr))
        ds_vel = load_dataset(vel_file)
        ds_vel.UVEL.data = ds_vel.UVEL.values
        ds_vel.VVEL.data = ds_vel.VVEL.values
        
        #Define output file name
        OW_outfile = join(outdir, 'monthly', 
                          'OW_k{}_{}{}.pdf'.format(str(k), monthstr, yearstr))
        
    elif seasonal:
        
        #Get seasonal velocity data
        vel_seas_file = join(seasonaldatdir, 
                             "avg_UVELVVEL_{}{}-{}{}.nc".format(
                                 season_start, yearstr, season_end, endyearstr))
        ds_vel = xr.open_mfdataset(vel_seas_file, engine="scipy")
        ds_vel.load()
        
        #Define output file name
        OW_outfile = join(outdir, 'seasonal', 
                          'OW_k{}_{}_{}.pdf'.format(
                              str(k), seas_monthstr, seas_yearstr))
 
    OW = get_OW_field(ds_grid, ds_vel, k, lats_lons, resolution, zeta_field)

    if not quiver:
        ArcCir_pcolormesh(ds_grid, [OW], resolution, 'seismic', lon_centers, 
                          lat_centers, None, datestr, 'OW', 
                          scalar_bounds=scalar_bounds, k_plot=k, extend='both',
                          outfile=OW_outfile, lats_lons=lats_lons)
        
    elif quiver: #need to fix outfile for this case
        ArcCir_pcolormesh_quiver(ds_grid, k, [OW], [vecE], [vecN], resolution, 
                                 'seismic', datestr, lon_centers, lat_centers,
                                 scalar_attr='OW', xvec_attr=xvec_attr, 
                                 scalar_bounds=scalar_bounds, 
                                 outfile=OW_outfile, lats_lons=lats_lons)

    OW_list.append(OW)
    return OW_list

##############################

def plot_Ek_vel_divergence(div_u_Ek_list, vecEs, vecNs, lon_centers, 
                           lat_centers, seasonal, outdir, k, yearstr, ds_grid, 
                           resolution, datestr, lats_lons, monthstr=None, 
                           datdir=None, season_start=None, season_end=None, 
                           endyearstr=None, seas_monthstr=None, 
                           seas_yearstr=None, seasonaldatdir=None, 
                           scalar_bounds=[-1, 1], quiver=False, vecE=None, 
                           vecN=None, xvec_attr=None):
    
    #Define output file name
    if not seasonal:
        div_u_Ek_outfile = join(outdir, 'monthly', 
                                'divuEk_uEk_k{}_{}{}.pdf'.format(
                                    str(k), monthstr, yearstr))
    elif seasonal:
        div_u_Ek_outfile = join(outdir, 'seasonal', 
                                'divuEk_uEk_k{}_{}_{}.pdf'.format(
                                    str(k), seas_monthstr, seas_yearstr))
    
    u_Ek, v_Ek = comp_temp_mean(vecEs), comp_temp_mean(vecNs)
    xgcm_grid = ecco.get_llc_grid(ds_grid)
    
    #Compute divergence
    lon_centers, lat_centers, lon_edges, lat_edges, div_u_Ek = \
    comp_2d_Ek_divergence(xgcm_grid, u_Ek, v_Ek, ds_grid.dxC, ds_grid.dyC, 
                          ds_grid.rA, ds_grid, lats_lons, resolution)

    if not quiver:
        ArcCir_pcolormesh(ds_grid, [div_u_Ek], resolution, 'PuOr', lon_centers, 
                          lat_centers, None, datestr, 'DIVUEk', 
                          scalar_bounds=scalar_bounds, k_plot=k, extend='both',
                          outfile=div_u_Ek_outfile, lats_lons=lats_lons)
    elif quiver: 
        ArcCir_pcolormesh_quiver(ds_grid, k, [div_u_Ek], [vecE], [vecN], 
                                 resolution, 'PuOr', datestr, lon_centers, 
                                 lat_centers, scalar_attr='DIVUEk', 
                                 xvec_attr='UEk', scalar_bounds=scalar_bounds, 
                                 outfile=div_u_Ek_outfile, lats_lons=lats_lons)
    div_u_Ek_list.append(div_u_Ek)
    return div_u_Ek_list

##############################

def plot_pcolormesh_k_plane(ds_grid, ds_scalar_list, k, scalar_attr, resolution,
                            cmap, datestr, vmin, vmax, outfile, lats_lons, 
                            datdir, Ro_l_list, OW_list, yearstr, year=None, 
                            outdir=None, monthstr=None, seas_monthstr=None, 
                            seas_yearstr=None, seasonal=False, 
                            multiple_seas=False, annual=False, 
                            season_start=None, season_end=None, endyearstr=None, 
                            season_years=None, years=None, startyr=None, 
                            datdirname=None, seasonaldatdir=None, 
                            data_seasons=None, lon_centers=None, 
                            lat_centers=None):
    """
    Creates pcolormesh plot on plane of constant k.
    """
    
    latmin, latmax, lonmin, lonmax = lats_lons
    
    #Take temporal average, if needed
    ds_scalar_mean = comp_temp_mean(ds_scalar_list)
    
    if type(ds_scalar_list[0]) == xr.Dataset: 
        
        if scalar_attr != 'DIVUEk' and scalar_attr != 'SIheff':
            ds_scalar_mean = ds_scalar_mean.isel(k=k) #Isolate k-plane
            
        ds_grid[scalar_attr] = ds_scalar_mean[scalar_attr]
        
        #Convert scalar DataSet to useful field
        lon_centers, lat_centers, lon_edges, lat_edges, scalar = ds_to_field(
            ds_grid, ds_scalar_mean, scalar_attr, latmin, latmax, lonmin, 
            lonmax, resolution)
    
    else: #Typically data are numpy arrays in this case
        scalar = ds_scalar_mean
    
    if seasonal:
        seas_yearstr = yearstr
        datestr = '{}, {}'.format(seas_monthstr, seas_yearstr)
        data_seasons.append(scalar)
    
    ArcCir_pcolormesh(ds_grid, [scalar], resolution, cmap, lon_centers, lat_centers, None, datestr, scalar_attr, scalar_bounds=[vmin, vmax], k_plot=k, extend='both', outfile=outfile, lats_lons=lats_lons)    
  
    if scalar_attr == 'ZETA': #If vorticity, also compute and plot Ro_l, OW

        if annual:

            #Plot Ro_l
            ArcCir_pcolormesh(ds_grid, Ro_l_list, resolution, 'Reds', lon_centers, lat_centers, None, yearstr, 'Ro_l', k_plot=k, extend='both', logscale=True, outfile=join(outdir, 'yearly', 'localRo_k{}_{}.pdf'.format(str(k), yearstr)), lats_lons=lats_lons)

            #Plot OW
            ArcCir_pcolormesh(ds_grid, OW_list, resolution, 'seismic', lon_centers, lat_centers, None, yearstr, 'OW', k_plot=k, extend='both', outfile=join(outdir, 'yearly', 'OW_k{}_{}.pdf'.format(str(k), yearstr)), lats_lons=lats_lons)
        
        elif multiple_seas: #If plotting interannual seasonal average
            
            seas_yearstr = str(startyr) + "-" + str(startyr + (years-1) + season_years[-1]) #For titles
     
            ArcCir_pcolormesh(ds_grid, Ro_l_list, resolution, 'Reds', lon_centers, lat_centers, None, '{}, {}'.format(seas_monthstr, seas_yearstr), 'Ro_l', k_plot=k, extend='both', logscale=True, outfile=join(outdir, 'interannual', 'localRo_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)

            ArcCir_pcolormesh(ds_grid, OW_list, resolution, 'seismic', lon_centers, lat_centers, None, '{}, {}'.format(seas_monthstr, seas_yearstr), 'OW', k_plot=k, extend='both', outfile=join(outdir, 'interannual', 'OW_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)
        
        elif not seasonal:
                
            #Compute and plot local Rossby number for the month
            Ro_l_list = plot_Ro_l(Ro_l_list, scalar, lon_centers, lat_centers, False, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, monthstr=monthstr)
                
            #Compute and plot OW for the month
            OW_list = plot_OW(OW_list, scalar, lon_centers, lat_centers, False, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, monthstr=monthstr, datdir=datdir)
        
        elif seasonal:    
            
            #Compute and plot local Rossby number for the season
            Ro_l_list = plot_Ro_l(Ro_l_list, scalar, lon_centers, lat_centers, True, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, seas_monthstr=seas_monthstr, seas_yearstr=seas_yearstr)
                
            #Make sure the seasonal velocity file exists already
                
            vel_seas_file = join(seasonaldatdir, "avg_UVELVVEL_"+season_start+yearstr+"-"+season_end+endyearstr+".nc") #Define filename
                
            if not os.path.exists(vel_seas_file): #If it doesn't exist, compute it
                save_seasonal_avgs.main(field='UVELVVEL', years=[year], start_month=season_start, end_month=season_end, usecompdata=False, datdir=datdirname, outdir=seasonaldatdir)
                
            #Compute and plot OW for the season
            OW_list = plot_OW(OW_list, scalar, True, yearstr, year, outdir, k, datdirname, ds_grid, lon_centers, lat_centers, latmin, latmax, lonmin, lonmax, resolution, datestr, lats_lons, season_start=season_start, season_end=season_end, endyearstr=endyearstr, seas_monthstr=seas_monthstr, seas_yearstr=seas_yearstr, seasonaldatdir=seasonaldatdir)
    
    if seasonal:
        return Ro_l_list, OW_list, data_seasons, lon_centers, lat_centers
    elif not seasonal and not annual and not multiple_seas:
        return Ro_l_list, OW_list
            
##############################

def plot_pcm_quiver_k_plane(ds_grid, ds_scalar_list, k, scalar_attr, xvec_attr, vecE, vecN, resolution, cmap, datestr, vmin, vmax, outfile, lats_lons, datdir, Ro_l_list, OW_list, div_u_Ek_list, yearstr, outdir=None, monthstr=None, seas_monthstr=None, seas_yearstr=None, seasonal=False, multiple_seas=False, annual=False, season_start=None, season_end=None, endyearstr=None, season_years=None, years=None, startyr=None, datdirname=None, seasonaldatdir=None, scalar_data_seasons=None, vecE_data_seasons=None, vecN_data_seasons=None, lon_centers=None, lat_centers=None):
    
    """
    Creates pcolormesh + quiver plot on plane of constant k.
    """
    
    latmin, latmax, lonmin, lonmax = lats_lons #Set spatial bounds 
    
    ds_scalar_mean = comp_temp_mean(ds_scalar_list) #Take temporal average, if needed
    
    if type(ds_scalar_list[0]) == xr.Dataset: 
        
        if scalar_attr != 'DIVUEk':
            ds_scalar_mean = ds_scalar_mean.isel(k=k) #Isolate k-plane
        
        ds_grid[scalar_attr] = ds_scalar_mean[scalar_attr]
        
        #Convert scalar DataSet to useful field
        lon_centers, lat_centers, lon_edges, lat_edges, scalar = ds_to_field(ds_grid, ds_scalar_mean, scalar_attr, latmin, latmax, lonmin, lonmax, resolution)
   
    else: #Typically data are numpy arrays in this case
        scalar = ds_scalar_mean #lat/lon_centers are to be input as kwargs in this case
    
    if seasonal:
        
        seas_yearstr = yearstr
        datestr = '{}, {}'.format(seas_monthstr, seas_yearstr)
        
        scalar_data_seasons.append(scalar)
        vecE_data_seasons.append(vecE)
        vecN_data_seasons.append(vecN)

    #Create main plot
    ArcCir_pcolormesh_quiver(ds_grid, k, [scalar], [vecE], [vecN], resolution, cmap, datestr, lon_centers, lat_centers, scalar_attr=scalar_attr, xvec_attr=xvec_attr, scalar_bounds=[vmin, vmax], extend='both', logscale=False, outfile=outfile, lats_lons=lats_lons, quiv_scale=0.3)
 
    if scalar_attr == 'ZETA': #If vorticity, also compute and plot Ro_l, OW overlaid with quiver
        
        yvec_attr = get_vector_partner(xvec_attr)
        vec = xvec_attr + yvec_attr

        if annual: #If plotting an annual average

            #Plot Ro_l
            ArcCir_pcolormesh_quiver(ds_grid, k, Ro_l_list, [vecE], [vecN], resolution, 'Reds', yearstr, lon_centers, lat_centers, scalar_attr='Ro_l', xvec_attr=xvec_attr, logscale=True, outfile=join(outdir, 'yearly', 'localRo_{}_k{}_{}.pdf'.format(get_variable_str(vec), str(k), yearstr)), lats_lons=lats_lons)
    
            #Plot OW
            ArcCir_pcolormesh_quiver(ds_grid, k, OW_list, [vecE], [vecN], resolution, 'seismic', yearstr, lon_centers, lat_centers, scalar_attr='OW', xvec_attr=xvec_attr, outfile=join(outdir, 'yearly', 'OW_{}_k{}_{}.pdf'.format(get_variable_str(vec), str(k), yearstr)), lats_lons=lats_lons)
        
        elif multiple_seas: #If plotting interannual seasonal average
            
            seas_yearstr = str(startyr) + "-" + str(startyr + (years-1) + season_years[-1]) #For titles
            
            #Plot Ro_l
            ArcCir_pcolormesh_quiver(ds_grid, k, Ro_l_list, [vecE], [vecN], resolution, 'Reds', '{}, {}'.format(seas_monthstr, seas_yearstr), lon_centers, lat_centers, scalar_attr='Ro_l', xvec_attr=xvec_attr, logscale=True, outfile=join(outdir, 'interannual', 'localRo_{}_k{}_{}_{}.pdf'.format(get_variable_str(vec), str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)
                 
            #Plot OW
            ArcCir_pcolormesh_quiver(ds_grid, k, OW_list, [vecE], [vecN], resolution, 'seismic', '{}, {}'.format(seas_monthstr, seas_yearstr), lon_centers, lat_centers, scalar_attr='OW', xvec_attr=xvec_attr, outfile=join(outdir, 'interannual', 'OW_{}_k{}_{}_{}.pdf'.format(get_variable_str(vec), str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)
             
        elif not seasonal: #Plot each month
            
            #Plot Ro_l
            Ro_l_list = plot_Ro_l(Ro_l_list, scalar, lon_centers, lat_centers, False, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, monthstr=monthstr, quiver=True, vecE=vecE, vecN=vecN, xvec_attr=xvec_attr)
                
            #Plot OW
            OW_list = plot_OW(OW_list, scalar, lon_centers, lat_centers, False, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, monthstr=monthstr, datdir=datdir, quiver=True, vecE=vecE, vecN=vecN, xvec_attr=xvec_attr)
    
        elif seasonal: #Plot each season
            
            #Compute and plot local Rossby number for the season
            Ro_l_list = plot_Ro_l(Ro_l_list, scalar, lon_centers, lat_centers, True, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, seas_monthstr=seas_monthstr, seas_yearstr=seas_yearstr)
                
            #Make sure the seasonal velocity file exists already
                
            vel_seas_file = join(seasonaldatdir, "avg_UVELVVEL_"+season_start+yearstr+"-"+season_end+endyearstr+".nc") #Define filename
                
            if not os.path.exists(vel_seas_file): #If it doesn't exist, compute it
                save_seasonal_avgs.main(field='UVELVVEL', years=[year], start_month=season_start, end_month=season_end, usecompdata=False, datdir=datdirname, outdir=seasonaldatdir)
                
            #Compute and plot OW for the season
            OW_list = plot_OW(OW_list, scalar, lon_centers, lat_centers, True, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, season_start=season_start, season_end=season_end, endyearstr=endyearstr, seas_monthstr=seas_monthstr, seas_yearstr=seas_yearstr, seasonaldatdir=seasonaldatdir)
            
    elif xvec_attr == 'UEk': #If Ekman current, also compute and plot its divergence

        if annual: #If plotting an annual average

            #Plot div_u_Ek
            ArcCir_pcolormesh_quiver(ds_grid, k, div_u_Ek_list, [vecE], [vecN], resolution, 'PuOr', yearstr, lon_centers, lat_centers, scalar_attr='DIVUEk', xvec_attr='UEk', outfile=join(outdir, 'yearly', 'divuEk_uEk_k{}_{}.pdf'.format(str(k), yearstr)), lats_lons=lats_lons)
        
        elif multiple_seas: #If plotting interannual seasonal average
            
            seas_yearstr = str(startyr) + "-" + str(startyr + (years-1) + season_years[-1]) #For titles
            
            #Plot div_u_Ek
            ArcCir_pcolormesh_quiver(ds_grid, k, div_u_Ek_list, [vecE], [vecN], resolution, 'PuOr', '{}, {}'.format(seas_monthstr, seas_yearstr), lon_centers, lat_centers, scalar_attr='DIVUEk', xvec_attr='UEk', outfile=join(outdir, 'interannual', 'divuEk_uEk_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)
          
        elif not seasonal: #Plot each month
            
            #Plot div_u_Ek
            div_u_Ek_list = plot_Ek_vel_divergence(div_u_Ek_list, [vecE], [vecN], lon_centers, lat_centers, False, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, monthstr=monthstr, quiver=True, vecE=vecE, vecN=vecN, xvec_attr='UEk')
           
        elif seasonal: #Plot each season
            
            #Compute and plot div_u_Ek for the season
            div_u_Ek_list = plot_Ek_vel_divergence(div_u_Ek_list, [vecE], [vecN], lon_centers, lat_centers, True, outdir, k, yearstr, ds_grid, resolution, datestr, lats_lons, quiver=True, vecE=vecE, vecN=vecN, xvec_attr='UEk', seas_monthstr=seas_monthstr, seas_yearstr=seas_yearstr)
    
    if seasonal:
        return Ro_l_list, OW_list, div_u_Ek_list, scalar_data_seasons, vecE_data_seasons, vecN_data_seasons, lon_centers, lat_centers
    elif not seasonal and not annual and not multiple_seas:
        return Ro_l_list, OW_list, div_u_Ek_list
    
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