"""
General functions for use with ECCO data.

Rosalie Cormier, 2024
"""

import os
import numpy as np
import xarray as xr
import ecco_v4_py as ecco

from os.path import join

from functions_ecco_download import ecco_podaac_download
from functions_field_variables import get_field_variable, get_vector_comps, field_is_primary, get_monthly_shortname, get_monthly_nc_string, get_seasonal_shortname, get_seasonal_nc_string

##############################

def compute_temporal_mean(timeseries):
    
    """
    Compute temporal mean of a field.
    """ 
    
    mean = (timeseries[0]).copy() / len(timeseries)
    
    if len(timeseries) > 1:
        for i in range(1, len(timeseries)):
            mean = mean + (timeseries[i]).copy() / len(timeseries)
        
    return mean

##############################

def get_monthstr(i):
    
    """
    Returns a string corresponding to month i.
    """
    
    month_dict = {1: "01", 2: "02", 3: "03", 4: "04", 5: "05", 6: "06", 7: "07", 8: "08", 9: "09", 10: "10", 11: "11", 12: "12"}
    
    return month_dict[((i-1) % 12) + 1]

##############################

def get_month_end(monthstr, yearstr):
    
    """
    Return string representing last day of specified month.
    """
    
    month_end_dict = {"01": "31", "03": "31", "04": "30", "05": "31", "06": "30", "07": "31", "08": "31", "09": "30", "10": "31", "11": "30", "12": "31"}
    
    if monthstr != "02":
        endmonth = month_end_dict[monthstr]
    elif monthstr == "02":
        if int(yearstr) % 4 == 0:
            endmonth = "29"
        else: 
            endmonth = "28"
    
    return endmonth

##############################

def get_args_from_date_string(date_string, time_ave_type, time_kwargs):
    
    if time_ave_type == 'monthly':
        
        month, year = date_string[5:9], date_string[0:4]
        #initial_month, initial_year = month, year
        #final_month, final_year = month, year
        
        return [month, year, month, year]
        
    elif time_ave_type == 'seasonal':
        
        season_start, season_end = time_kwargs[0], time_kwargs[1]
        
        if int(season_start) < int(season_end):
            year = date_string[6:10]
            initial_month, initial_year = date_string[0:2], year
            final_month, final_year = date_string[3:5], year
        elif int(season_end) < int(season_start):
            initial_month, initial_year = date_string[0:2], date_string[6:10]
            final_month, final_year = date_string[3:5], date_string[11:15]
            
        return [initial_month, initial_year, final_month, final_year]

##############################

def load_grid(datdir_primary):
    
    """
    Loads ECCO grid.
    """
    
    grid_params_shortname = "ECCO_L4_GEOMETRY_LLC0090GRID_V4R4"
    grid_params_file = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
    grid_params_directory = join(datdir_primary, grid_params_shortname)

    if not os.path.exists(grid_params_directory): 
        
        os.makedirs(grid_params_directory)
    
        #Download ECCO grid parameters (date is arbitrary)
        ecco_podaac_download(ShortName=grid_params_shortname, StartDate="2000-01-01", EndDate="2000-01-02", download_root_dir=datdir_primary, n_workers=6, force_redownload=False)

    ds_grid = xr.open_dataset(join(grid_params_directory, grid_params_file)) #Load grid parameters
    
    return ds_grid

##############################

def load_primary_data_file(field_name, date_string, datdir_primary, time_ave_type):

    """
    Load a specified ECCO (primary) DataSet, or a time-averaged primary-data DataSet.
    """
    
    if time_ave_type == 'monthly':
        field_shortname, field_nc_string = get_monthly_shortname(get_field_variable(field_name)), get_monthly_nc_string(get_field_variable(field_name))
        file_suffix = '_ECCO_V4r4_native_llc0090.nc'
        
    elif time_ave_type == 'seasonal':
        field_shortname, field_nc_string = get_seasonal_shortname(get_field_variable(field_name)), get_seasonal_nc_string(get_field_variable(field_name))
        file_suffix = '.nc'
    
    data_file_path = join(datdir_primary, field_shortname, field_nc_string+date_string+file_suffix)
    print(data_file_path)
    try: #This option should work for files (e.g. monthly averages) that come directly from ECCO
        dataset = xr.open_mfdataset(data_file_path, parallel=True, data_vars='minimal', coords='minimal', compat='override')
    except: #This case should catch the computed (e.g. seasonal) averages
        dataset = xr.open_mfdataset(data_file_path, engine="scipy")
    
    return dataset

##############################

def rotate_vector(curr_ds_grid, vector_ds, vector_field_name, vector_comps):
    
    """
    Get eastward and northward components of vector in x-y coordinates.
    """

    if field_is_primary(vector_field_name):
        
        vector_ds[vector_comps[0]].data = vector_ds[vector_comps[0]].values
        vector_ds[vector_comps[1]].data = vector_ds[vector_comps[1]].values
        
        xgcm_grid = ecco.get_llc_grid(curr_ds_grid)
        interp_vector = xgcm_grid.interp_2d_vector({'X': vector_ds[vector_comps[0]], 'Y': vector_ds[vector_comps[1]]}, boundary='fill')
        vel_E_comp = interp_vector['X'] * curr_ds_grid['CS'] - interp_vector['Y'] * curr_ds_grid['SN']
        vel_N_comp = interp_vector['X'] * curr_ds_grid['SN'] + interp_vector['Y'] * curr_ds_grid['CS']
        
    elif not field_is_primary(vector_field_name):
        
        vector_ds[vector_comps[0]].data = vector_ds[vector_comps[0]].values
        vector_ds[vector_comps[1]].data = vector_ds[vector_comps[1]].values
        
        vel_E_comp = vector_ds[vector_comps[0]] * curr_ds_grid['CS'] - vector_ds[vector_comps[1]] * curr_ds_grid['SN']
        vel_N_comp = vector_ds[vector_comps[0]] * curr_ds_grid['SN'] + vector_ds[vector_comps[1]] * curr_ds_grid['CS']
    
    return vel_E_comp, vel_N_comp

##############################

def scalar_to_grid(ds_grid, scalar_ds, field_variable, depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res):
    
    """
    Resample scalar DataSet attribute at specific k-value (depth) to lat-lon grid.
    If field is vertical velocity (w), interpolate along z-axis (helpful for visualization).
    """
    
    curr_ds_grid = ds_grid.copy()
    
    if field_variable == 'WVEL': #If w, interpolate vertically
        xgcm_grid = ecco.get_llc_grid(curr_ds_grid)
        scalar_ds['WVEL'] = xgcm_grid.interp(scalar_ds.WVEL, axis='Z')
    
    scalar_ds = scalar_ds.isel(k=int(depth)) #Isolate plane at specified depth
    
    curr_ds_grid[field_variable] = scalar_ds[field_variable].squeeze()
    curr_ds_grid.load()
    
    latmin, latmax, lonmin, lonmax = float(latmin), float(latmax), float(lonmin), float(lonmax)
    lat_res, lon_res = float(lat_res), float(lon_res)
    
    lon_centers, lat_centers, lon_edges, lat_edges, field = ecco.resample_to_latlon(curr_ds_grid.XC, \
                                            curr_ds_grid.YC, curr_ds_grid[field_variable], latmin, latmax, \
                                            lat_res, lonmin, lonmax, lon_res, fill_value=np.NaN, \
                                            mapping_method='nearest_neighbor', radius_of_influence=120000)
    
    return lon_centers, lat_centers, lon_edges, lat_edges, field

##############################

def vector_to_grid(ds_grid, vector_ds, vector_field_name, depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res):
    
    """
    Interpolate/rotate 2D vector from DataSet (x-y coordinates) onto east-north axes.
    Resample vector components at specific k-value (depth) to lat-lon grid.
    """
    
    curr_ds_grid = ds_grid.copy()
    
    vector_comps = get_vector_comps(vector_field_name)
    vec_E_comp, vec_N_comp = rotate_vector(curr_ds_grid, vector_ds.isel(k=int(depth)), vector_field_name, vector_comps) #Rotate vector
    vec_E_comp, vec_N_comp = vec_E_comp.squeeze(), vec_N_comp.squeeze() #Squeeze along time axis
    
    curr_ds_grid[vector_comps[0]] = vec_E_comp
    curr_ds_grid[vector_comps[1]] = vec_N_comp
    curr_ds_grid.load()
    
    latmin, latmax, lonmin, lonmax = float(latmin), float(latmax), float(lonmin), float(lonmax)
    lat_res, lon_res = float(lat_res), float(lon_res)
    
    lon_centers, lat_centers, lon_edges, lat_edges, field_E_comp = ecco.resample_to_latlon(curr_ds_grid.XG, \
                                            curr_ds_grid.YG, curr_ds_grid[vector_comps[0]], latmin, latmax, \
                                            lat_res, lonmin, lonmax, lon_res, fill_value=np.NaN, \
                                            mapping_method='nearest_neighbor', radius_of_influence=120000)
    field_N_comp = ecco.resample_to_latlon(ds_grid.XG, ds_grid.YG, curr_ds_grid[vector_comps[1]], latmin, latmax, \
                                            lat_res, lonmin, lonmax, lon_res, fill_value=np.NaN, \
                                            mapping_method='nearest_neighbor', radius_of_influence=120000)[4]
    
    return lon_centers, lat_centers, lon_edges, lat_edges, field_E_comp, field_N_comp