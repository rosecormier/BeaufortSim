"""
Contains functions (associated with computed data) to:
    -Create a file containing computed data (if it doesn't already exist);
    -Load a computed DataSet.

Rosalie Cormier, 2023
"""

import os
import xarray as xr

from os.path import join

import NEWcomp_data

from functions_ecco_general import get_monthstr
from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string

##############################

def create_comp_data_file(field_name, initial_month, initial_year, final_month, final_year, datdir_primary, datdir_secondary, time_ave_type='monthly'):
    
    """
    Checks that a computed file exists, and creates it if it doesn't.
    Note - only does monthly avgs at the moment - will fix in future (i.e., will make time_ave_type variable)
    """
    
    month, year = int(initial_month), int(initial_year)
    
    if time_ave_type == 'monthly':
    
        field_shortname = get_monthly_shortname(get_field_variable(field_name))
        field_nc_string = get_monthly_nc_string(get_field_variable(field_name))

        while year < int(final_year):
            
            yearstr = str(year)
            
            while month <= 12:
                
                monthstr = get_monthstr(month-1) #The indexing is weird for months
                date_string = yearstr + '-' + monthstr

                filename = field_nc_string + date_string + '.nc'
                path_to_file = os.path.join(datdir_secondary, field_shortname, filename)

                if not os.path.exists(path_to_file): #Execute only if the file doesn't already exist
                    
                    #Create the file
                    NEWcomp_data.main(datdir_primary=datdir_primary, datdir_secondary=datdir_secondary, monthstr=monthstr, yearstr=yearstr, field_name=field_name)

                month += 1
                
            year += 1
            month = 1
        
        if year == int(final_year):

            while month <= int(final_month):

                monthstr = get_monthstr(month-1) #The indexing is weird for months
                date_string = final_year + '-' + monthstr

                filename = field_nc_string + date_string + '.nc'
                path_to_file = os.path.join(datdir_secondary, field_shortname, filename)

                if not os.path.exists(path_to_file): #Execute only if the file doesn't already exist
                    
                    #Create the file
                    NEWcomp_data.main(datdir_primary=datdir_primary, datdir_secondary=datdir_secondary, monthstr=monthstr, yearstr=final_year, field_name=field_name)
                
                month += 1
                
    print("Done computing secondary data.")
        
##############################

def load_comp_data_file():
    
    """
    Loads computed DataSet.
    """
    
#    #define variable 'filename' here
    
#    try:
#        comp_ds = xr.open_mfdataset(filename, engine="scipy") #Load DataSet
#        return comp_ds
    
#    except:
#        print("DataSet does not exist.")