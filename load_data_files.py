"""
Rosalie Cormier, 2024
"""

import xarray as xr

from os.path import join

from functions_ecco_general import get_monthstr, load_grid
from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string, get_seasonal_shortname, get_seasonal_nc_string

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
    
    dataset.load()
        
    return dataset

##############################

def load_secondary_data_file(field_name, date_string, datdir_secondary, time_ave_type):
    
    """
    Loads DataSet for a given computed field at a given time.
    """
    
    if time_ave_type == 'monthly':
        field_shortname, field_nc_string = get_monthly_shortname(get_field_variable(field_name)), get_monthly_nc_string(get_field_variable(field_name))
    elif time_ave_type == 'seasonal':
        field_shortname, field_nc_string = get_seasonal_shortname(get_field_variable(field_name)), get_seasonal_nc_string(get_field_variable(field_name))
        
    filename = field_nc_string + date_string + '.nc'
    path_to_file = join(datdir_secondary, field_shortname, filename)
    
    try: #Try to load DataSet
        computed_ds = xr.open_mfdataset(path_to_file, engine="scipy")
        print("Loaded DataSet.")
        return computed_ds
    
    except:
        print("DataSet does not exist.")