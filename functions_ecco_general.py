"""
General functions for use with ECCO data.

Rosalie Cormier, 2023
"""

import os
import numpy as np
import xarray as xr
import ecco_v4_py as ecco

from os.path import join

from functions_ecco_download import ecco_podaac_download

def load_grid(datdir):
    
    """
    Loads ECCO grid.
    """
    
    grid_params_shortname = "ECCO_L4_GEOMETRY_LLC0090GRID_V4R4"
    grid_params_file = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
    grid_params_directory = join(datdir, grid_params_shortname)

    if not os.path.exists(grid_params_directory): 
        
        os.makedirs(grid_params_directory)
    
        #Download ECCO grid parameters (date is arbitrary)
        ecco_podaac_download(ShortName=grid_params_shortname, StartDate="2000-01-01", \
                     EndDate="2000-01-02", download_root_dir=datdir, n_workers=6, \
                     force_redownload=False)

    ds_grid = xr.open_dataset(join(grid_params_directory, grid_params_file)) #Load grid parameters
    
    return ds_grid

def get_vector_partner(x_comp):
    
    y_comps = {'UVEL': 'VVEL', \
              'UG': 'VG', \
              'EXFtaux': 'EXFtauy', \
              'UEk': 'VEk'}
    y_comp = y_comps[x_comp]
    
    return y_comp

def get_monthstr(i):
    
    """
    Returns a string corresponding to month i.
    """
    
    month_dict = {0: "01", 1: "02", 2: "03", 3: "04", 4: "05", 5: "06",
              6: "07", 7: "08", 8: "09", 9: "10", 10: "11", 11: "12"}
    
    return month_dict[i % 12]

def get_month_name(monthstr):
    
    monthnames = {"01": [0, "January"], "02": [1, "February"], "03": [2, "March"], "04": [3, "April"], "05": [4, "May"], \
                  "06": [5, "June"], "07": [6, "July"], "08": [7, "August"], "09": [8, "September"], "10": [9, "October"], \
                  "11": [10, "November"], "12": [11, "December"]}

    month_i, monthname = monthnames[monthstr]
    
    return month_i, monthname

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

def load_dataset(curr_file):

    """
    Opens ECCO dataset given a file.
    """
    
    dataset = xr.open_mfdataset(curr_file, parallel=True, data_vars='minimal', coords='minimal', compat='override')
    dataset.load()
    
    return dataset

def get_scalar_in_xy(ecco_ds_grid, ecco_ds_scalar, scalar_attr):
    
    """
    Loads scalar field in xy-grid.
    """
    
    ds_grid = ecco_ds_grid.copy()
    ds_grid[scalar_attr] = ecco_ds_scalar[scalar_attr]
    ds_grid = ds_grid.load()
    
    return ds_grid
    
def get_vector_in_xy(ecco_ds_grid, ecco_ds_vector, xvec_attr, yvec_attr): 
    
    """
    Loads vector field in xy-grid.
    """

    ds_grid = ecco_ds_grid.copy()
    ds_grid = ds_grid.load()
    XGCM_grid = ecco.get_llc_grid(ds_grid)

    velc = XGCM_grid.interp_2d_vector({'X': ecco_ds_vector[xvec_attr], \
                                           'Y': ecco_ds_vector[yvec_attr]}, \
                                           boundary='fill')

    return velc

def rotate_vector(ecco_ds_grid, ecco_ds_vector, xvec_attr, yvec_attr):
    
    """
    Gets eastward and northward components of xy-vector.
    """
    
    velc = get_vector_in_xy(ecco_ds_grid, ecco_ds_vector, xvec_attr, yvec_attr)
    velE = velc['X'] * ecco_ds_grid['CS'] - velc['Y'] * ecco_ds_grid['SN']
    velN = velc['X'] * ecco_ds_grid['SN'] + velc['Y'] * ecco_ds_grid['CS']
    
    return velE, velN

def ds_to_field(ecco_ds_grid, ecco_ds_scalar, scalar_attr, latmin, latmax, lonmin, lonmax, resolution):
    
    """
    Resamples scalar DataSet attribute to lat-lon grid.
    """
    
    ds_grid = get_scalar_in_xy(ecco_ds_grid, ecco_ds_scalar, scalar_attr)
    curr_field = (ds_grid[scalar_attr]).squeeze()
    
    #ds_grid = ecco_ds_grid.copy()
    
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, \
    field_nearest = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, curr_field, latmin, latmax, resolution, \
                                            lonmin, lonmax, resolution, fill_value=np.NaN, \
                                            mapping_method='nearest_neighbor', radius_of_influence=120000)
    
    return new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, field_nearest

def comp_temp_mean(timeseries):
    
    """
    Computes temporal mean of a field.
    """ 
    
    mean = (timeseries[0]).copy() / len(timeseries)
    
    if len(timeseries) > 1:
        for i in range(1, len(timeseries)):
            mean = mean + (timeseries[i]).copy() / len(timeseries)
        
    return mean

def comp_residuals(fields, mean):
    
    """
    Computes residuals relative to a mean.
    """
        
    residual_list = []
        
    for field in fields:
            
        residual = field.copy() - mean.copy()
        residual_list.append(residual)
        
    return residual_list

def ecco_resample(ds_grid, curr_field, latmin, latmax, lonmin, lonmax, resolution):
    
    """
    Resamples field to lat-lon grid.
    """
    
    lon_centers, lat_centers, lon_edges, lat_edges, field = ecco.resample_to_latlon(ds_grid.XG, ds_grid.YG, curr_field, latmin, latmax, resolution, lonmin, lonmax, resolution, fill_value=np.NaN, mapping_method='nearest_neighbor', radius_of_influence=120000)
    
    return lon_centers, lat_centers, lon_edges, lat_edges, field

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