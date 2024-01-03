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
from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string

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
        ecco_podaac_download(ShortName=grid_params_shortname, StartDate="2000-01-01", \
                     EndDate="2000-01-02", download_root_dir=datdir, n_workers=6, \
                     force_redownload=False)

    ds_grid = xr.open_dataset(join(grid_params_directory, grid_params_file)) #Load grid parameters
    
    return ds_grid

##############################

def get_monthstr(i):
    
    """
    Returns a string corresponding to month i.
    """
    
    month_dict = {1: "01", 2: "02", 3: "03", 4: "04", 5: "05", 6: "06",
              7: "07", 8: "08", 9: "09", 10: "10", 11: "11", 12: "12"}
    
    return month_dict[((i-1) % 12) + 1]

##############################

def get_month_name(monthstr):
    
    monthnames = {"01": [0, "January"], "02": [1, "February"], "03": [2, "March"], "04": [3, "April"], "05": [4, "May"], \
                  "06": [5, "June"], "07": [6, "July"], "08": [7, "August"], "09": [8, "September"], "10": [9, "October"], \
                  "11": [10, "November"], "12": [11, "December"]}

    month_i, monthname = monthnames[monthstr]
    
    return month_i, monthname

##############################

def get_month_end(monthstr, yearstr):
    
    """
    Return string representing last day of specified month.
    """
    
    month_end_dict = {"01": "31", "03": "31", "04": "30",
                 "05": "31", "06": "30", "07": "31", "08": "31",
                 "09": "30", "10": "31", "11": "30", "12": "31"}
    
    if monthstr != "02":
        endmonth = month_end_dict[monthstr]
        
    elif monthstr == "02":
        if int(yearstr) % 4 == 0:
            endmonth = "29"
        else: 
            endmonth = "28"
    
    return endmonth

##############################

def load_ECCO_data_file(field_name, date_string, datdir_primary, time_ave_type):

    """
    Load a specified ECCO (primary) DataSet.
    """
    
    if time_ave_type == 'monthly': #will update to include other options
        field_shortname, field_nc_string = get_monthly_shortname(get_field_variable(field_name)), get_monthly_nc_string(get_field_variable(field_name))
    
    data_file = join(datdir_primary, field_shortname, field_nc_string+date_string+"_ECCO_V4r4_native_llc0090.nc")
    
    dataset = xr.open_mfdataset(data_file, parallel=True, data_vars='minimal', coords='minimal', compat='override')
    dataset.load()
    
    return dataset

##############################

def rotate_vector(curr_ds_grid, vector_ds, vector_comps):
    
    """
    Get eastward and northward components of vector in x-y coordinates.
    """
    
    #xgcm_grid = ecco.get_llc_grid(curr_ds_grid)
    #interp_vector = xgcm_grid.interp_2d_vector({'X': vector_ds[vector_comps[0]], 'Y': vector_ds[vector_comps[1]]}, boundary='fill')

    vel_E_comp = vector_ds[vector_comps[0]] * curr_ds_grid['CS'] - vector_ds[vector_comps[1]] * curr_ds_grid['SN']
    vel_N_comp = vector_ds[vector_comps[0]] * curr_ds_grid['SN'] + vector_ds[vector_comps[1]] * curr_ds_grid['CS']
    #vel_N_comp = interp_vector['X'] * curr_ds_grid['SN'] + interp_vector['Y'] * curr_ds_grid['CS']
    
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

def vector_to_grid(ds_grid, vector_ds, vector_comps, depth, latmin, latmax, lonmin, lonmax, lat_res, lon_res):
    
    """
    Interpolate/rotate 2D vector from DataSet (x-y coordinates) onto east-north axes.
    Resample vector components at specific k-value (depth) to lat-lon grid.
    """
    
    curr_ds_grid = ds_grid.copy()
    
    vec_E_comp, vec_N_comp = rotate_vector(curr_ds_grid, vector_ds, vector_comps) #Rotate vector
    
    #Isolate plane at specified depth and squeeze along time axis
    vec_E_comp, vec_N_comp = vec_E_comp.isel(k=int(depth)).squeeze(), vec_N_comp.isel(k=int(depth)).squeeze() 
    
    curr_ds_grid[vector_comps[0]] = vec_E_comp
    curr_ds_grid[vector_comps[1]] = vec_N_comp
    curr_ds_grid.load()
    
    latmin, latmax, lonmin, lonmax = float(latmin), float(latmax), float(lonmin), float(lonmax)
    lat_res, lon_res = float(lat_res), float(lon_res)
    
    lon_centers, lat_centers, lon_edges, lat_edges, field_E_comp = ecco.resample_to_latlon(ds_grid.XG, \
                                            ds_grid.YG, curr_ds_grid[vector_comps[0]], latmin, latmax, \
                                            lat_res, lonmin, lonmax, lon_res, fill_value=np.NaN, \
                                            mapping_method='nearest_neighbor', radius_of_influence=120000)
    lon_centers, lat_centers, lon_edges, lat_edges, field_N_comp = ecco.resample_to_latlon(ds_grid.XG, \
                                            ds_grid.YG, curr_ds_grid[vector_comps[1]], latmin, latmax, \
                                            lat_res, lonmin, lonmax, lon_res, fill_value=np.NaN, \
                                            mapping_method='nearest_neighbor', radius_of_influence=120000)
    
    return lon_centers, lat_centers, lon_edges, lat_edges, field_E_comp, field_N_comp

##############################

def get_season_months_and_years(start_month, end_month):
    
    months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

    season_start_i, season_end_i = months.index(start_month), months.index(end_month)
    
    if season_end_i >= season_start_i:
    
        season_months = months[season_start_i:season_end_i+1]
        season_years = []

        for month in season_months:
            season_years.append(0)

    elif season_end_i < season_start_i:

        season_1 = months[season_start_i:]
        season_2 = months[0:season_end_i+1]

        season_years = []

        for month in season_1:
            season_years.append(0)

        for month in season_2:
            season_years.append(1)

        season_months = season_1 + season_2
        
    return season_months, season_years