"""
Contains functions (associated with computed data) to:
    -Create a file containing computed data (if it doesn't already exist);
    -Load a computed DataSet.

Rosalie Cormier, 2023
"""

import os
import xarray as xr

from os.path import join

from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string

##############################

def create_comp_data_file(field_name, date_string, datdir_secondary, time_ave_type='monthly'):
    
    """
    Checks that a computed file exists, and creates it if it doesn't.
    Note - only does monthly avgs at the moment - will fix in future (i.e., will make time_ave_type variable)
    """
    
    if time_ave_type == 'monthly':
        date_string = yearstr + '-' + monthstr
        field_shortname = get_monthly_shortname(get_field_variable(field_name))
        field_nc_string = get_monthly_nc_string(get_field_variable(field_name))

    filename = field_nc_string + date_string + '.nc'
    path_to_file = os.path.join(datdir_secondary, field_shortname, filename)
    
    #if not os.path.exists(path_to_file): #Execute only if the file doesn't already exist
        #create the file - to be added
        
##############################

#def load_comp_data_file():
    
#    """
#    Loads computed DataSet.
#    """
    
#    #define variable 'filename' here
    
#    try:
#        comp_ds = xr.open_mfdataset(filename, engine="scipy") #Load DataSet
#        return comp_ds
    
#    except:
#        print("DataSet does not exist.")