"""
Contains functions to load primary (ECCO) and secondary (computed) datafiles.
The main function receives a field name and calls the correct function to load 
the data.

Rosalie Cormier, 2024
"""

import xarray as xr

from os.path import join

from functions_ecco_general import get_monthstr, load_grid
from functions_field_variables import get_field_variable, field_is_primary, \
get_monthly_shortname, get_monthly_nc_string, get_seasonal_shortname, \
get_seasonal_nc_string

##############################

def load_primary_data_file(field_name, date_string, datdir_primary, 
                           time_ave_type):

    """
    Load a specified ECCO (primary) DataSet, or a time-averaged primary-data 
    DataSet.
    """
    
    if time_ave_type == 'monthly':
        field_shortname = get_monthly_shortname(get_field_variable(field_name))
        field_nc_string = get_monthly_nc_string(get_field_variable(field_name))
        file_suffix = '_ECCO_V4r4_native_llc0090.nc'
    elif time_ave_type == 'seasonal':
        field_shortname = get_seasonal_shortname(get_field_variable(field_name))
        field_nc_string = get_seasonal_nc_string(get_field_variable(field_name))
        file_suffix = '.nc'
    
    data_file_path = join(datdir_primary, field_shortname, 
                          field_nc_string+date_string+file_suffix)

    try: 
        #This option should work for files (e.g. monthly averages) that come 
        #directly from ECCO
        dataset = xr.open_mfdataset(data_file_path, parallel=True, 
                                    data_vars='minimal', coords='minimal', 
                                    compat='override')
    except: 
        #This case should catch the computed (e.g. seasonal) averages
        dataset = xr.open_mfdataset(data_file_path, engine="scipy")
    
    dataset.load()
        
    return dataset

##############################

def load_secondary_data_file(field_name, date_string, datdir_secondary, 
                             time_ave_type):
    
    """
    Load a specified computed DataSet consisting of secondary field variables.
    """
    
    if time_ave_type == 'monthly':
        field_shortname = get_monthly_shortname(get_field_variable(field_name))
        field_nc_string = get_monthly_nc_string(get_field_variable(field_name))
    elif time_ave_type == 'seasonal':
        field_shortname = get_seasonal_shortname(get_field_variable(field_name))
        field_nc_string = get_seasonal_nc_string(get_field_variable(field_name))
        
    filename = field_nc_string + date_string + '.nc'
    path_to_file = join(datdir_secondary, field_shortname, filename)
    
    try: #Try to load DataSet
        computed_ds = xr.open_mfdataset(path_to_file, engine="scipy")
        return computed_ds
    
    except:
        print("DataSet does not exist.")
        
##############################

def main(**kwargs):
   
    if kwargs:
        
        field_name = kwargs.get('field_name')
        
        time_ave_type = kwargs.get('time_ave_type')
        date_string = kwargs.get('date_string')
        
        datdir_primary = kwargs.get('datdir_primary')
        datdir_secondary = kwargs.get('datdir_secondary')
    
    if field_is_primary(field_name):
        dataset = load_primary_data_file(field_name, date_string, 
                                         datdir_primary, time_ave_type)
        
    elif not field_is_primary(field_name):
        dataset = load_secondary_data_file(field_name, date_string, 
                                           datdir_secondary, time_ave_type)

    return dataset
        
##############################

if __name__ == "__main__":
    main()