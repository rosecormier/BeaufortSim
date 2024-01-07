"""
Contains functions (associated with computed data) to:
    -Create a file containing computed data (if it doesn't already exist);
    -Load a computed DataSet.

Rosalie Cormier, 2023
"""

import os
import xarray as xr

from os.path import join

import comp_secondary

from functions_ecco_general import get_monthstr, load_grid
from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string, get_seasonal_shortname, get_seasonal_nc_string

##############################

def create_secondary_data_file(field_name, date_string, datdir_primary, datdir_secondary, rho_ref, nu_E, time_ave_type, time_kwargs):
    
    """
    Iterates over time, checks that computed files for a given field exist, and creates them if they don't.
    """
    
    if time_ave_type == 'monthly':
        field_shortname, field_nc_string = get_monthly_shortname(get_field_variable(field_name)), get_monthly_nc_string(get_field_variable(field_name))
    elif time_ave_type == 'seasonal':
        field_shortname, field_nc_string = get_seasonal_shortname(get_field_variable(field_name)), get_seasonal_nc_string(get_field_variable(field_name))
        
    filename = field_nc_string + date_string + '.nc'
    path_to_file = os.path.join(datdir_secondary, field_shortname, filename)

    if not os.path.exists(path_to_file): #Create file if it doesn't exist
        comp_secondary.main(datdir_primary=datdir_primary, datdir_secondary=datdir_secondary, date_string=date_string, time_ave_type=time_ave_type, time_kwargs=time_kwargs, field_name=field_name, rho_ref=rho_ref, nu_E=nu_E)

    print("Done computing secondary data.")
        
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