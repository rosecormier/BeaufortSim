"""
Rosalie Cormier, 2023
"""

import numpy as np
import xarray as xr
import ecco_v4_py as ecco

from os.path import join
from ecco_download import ecco_podaac_download

def load_grid(datdir):
    
    """
    Loads ECCO grid.
    
    datdir = directory where data is stored
    """
    
    grid_params_shortname = "ECCO_L4_GEOMETRY_LLC0090GRID_V4R4"

    #Download ECCO grid parameters (date is arbitrary)
    ecco_podaac_download(ShortName=grid_params_shortname, StartDate="2000-01-01", \
                     EndDate="2000-01-02", download_root_dir=datdir, n_workers=6, \
                     force_redownload=False)

    grid_params_file = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
    grid_params_file_path = join(datdir, grid_params_shortname, grid_params_file)

    ds_grid = xr.open_dataset(grid_params_file_path) #Load grid parameters
    
    return ds_grid

def get_starting_i(startmo):

    """
    Get index representing start month.
    
    startmo = starting month
    """
    
    month_dict = {0: "01", 1: "02", 2: "03", 3: "04", 4: "05", 5: "06",
              6: "07", 7: "08", 8: "09", 9: "10", 10: "11", 11: "12"}
    month_key_list, month_val_list = list(month_dict.keys()), list(month_dict.values())
    
    month_index = month_val_list.index(startmo)
    i = month_key_list[month_index]
    
    return i

def get_monthstr(i):
    
    """
    Returns a string corresponding to month i.
    """
    
    month_dict = {0: "01", 1: "02", 2: "03", 3: "04", 4: "05", 5: "06",
              6: "07", 7: "08", 8: "09", 9: "10", 10: "11", 11: "12"}
    
    return month_dict[i % 12]

def get_month_name(monthstr):
    
    monthnames = {"01": "January", "02": "February", "03": "March", "04": "April", "05": "May", \
                  "06": "June", "07": "July", "08": "August", "09": "September", "10": "October", \
                  "11": "November", "12": "December"}
    monthname = monthnames[monthstr]
    return monthname

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
    
    return dataset

def get_scalar_in_xy(ecco_ds_grid, k_val, ecco_ds_scalar, scalar_attr):
    
    """
    Loads scalar field in xy-grid.
    
    ecco_ds_grid = ECCO grid
    k_val = depth index of interest
    ecco_ds_scalar = DataSet containing field
    scalar_attr = name of field
    """
    
    ds_grid = ecco_ds_grid.copy()
    ds_grid[scalar_attr] = ecco_ds_scalar[scalar_attr]
    ds_grid = ds_grid.load()
    
    return ds_grid
    
def get_vector_in_xy(ecco_ds_grid, k_val, ecco_ds_vector, xvec_attr, yvec_attr):
    
    """
    Loads vector field in xy-grid.
    
    ecco_ds_grid = ECCO grid
    k_val = depth index of interest
    ecco_ds_vector = DataSet containing vector field
    xvec_attr = name of x-comp of vector field
    yvec_attr = name of y-comp of vector field
    """

    ds_grid = ecco_ds_grid.copy()
    ds_grid = ds_grid.load()
    XGCM_grid = ecco.get_llc_grid(ds_grid)

    velc = XGCM_grid.interp_2d_vector({'X': ecco_ds_vector[xvec_attr], \
                                           'Y': ecco_ds_vector[yvec_attr]}, \
                                           boundary='fill')

    return velc

def rotate_vector(ecco_ds_grid, k_val, ecco_ds_vector, xvec_attr, yvec_attr):
    
    """
    Gets eastward and northward components of xy-vector.
    
    ecco_ds_grid = grid DataSet
    k_val = depth value of index
    ecco_ds_vector = DataSet containing vector
    x/yvec_attr = attributes corresponding to vector components
    """
    
    velc = get_vector_in_xy(ecco_ds_grid, k_val, ecco_ds_vector, xvec_attr, yvec_attr)
    velE = velc['X'] * ecco_ds_grid['CS'] - velc['Y'] * ecco_ds_grid['SN']
    velN = velc['X'] * ecco_ds_grid['SN'] + velc['Y'] * ecco_ds_grid['CS']
    
    return velE, velN

def ds_to_field(ecco_ds_grid, ecco_ds_scalar, scalar_attr, k_val, latmin, latmax, lonmin, lonmax, resolution):
    
    """
    Resamples scalar DataSet attribute to lat-lon grid
    """
    
    ds_grid = get_scalar_in_xy(ecco_ds_grid, k_val, ecco_ds_scalar, scalar_attr)
    curr_field = (ds_grid[scalar_attr]).squeeze()
    
    ds_grid = ecco_ds_grid.copy()
    
    new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, \
    field_nearest = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, curr_field, latmin, latmax, resolution, \
                                            lonmin, lonmax, resolution, fill_value=np.NaN, \
                                            mapping_method='nearest_neighbor', radius_of_influence=120000)
    
    return new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, field_nearest

def comp_temp_mean(timeseries):
    
    """
    Computes temporal mean of a field.
    """ 
    
    mean = (timeseries[0]).copy()
    
    for i in range(len(timeseries)):
        mean = mean + (timeseries[i]).copy() / len(timeseries)
        
    return mean

def comp_residuals(fields, mean):
    
    """
    Computes residuals relative to a mean.
    
    fields = data to compute residuals for
    mean = mean of all in fields
    """
        
    residual_list = []
        
    for field in fields:
            
        residual = field.copy() - mean.copy()
        residual_list.append(residual)
        
    return residual_list