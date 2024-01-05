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

from functions_ecco_general import get_monthstr
from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string

##############################

def create_secondary_data_file(field_name, initial_month, initial_year, final_month, final_year, datdir_primary, datdir_secondary, rho_ref, nu_E, time_ave_type):
    
    """
    Iterates over time, checks that computed files for a given field exist, and creates them if they don't.
    Note - only does monthly avgs at the moment - will fix in future (i.e., will add conditions for other time_ave_type)
    """
    
    month, year = int(initial_month), int(initial_year)
    
    if time_ave_type == 'monthly':
    
        field_shortname = get_monthly_shortname(get_field_variable(field_name))
        field_nc_string = get_monthly_nc_string(get_field_variable(field_name))

        while year < int(final_year):
            
            yearstr = str(year)
            
            while month <= 12:
                
                monthstr = get_monthstr(month)
                date_string = yearstr + '-' + monthstr

                filename = field_nc_string + date_string + '.nc'
                path_to_file = os.path.join(datdir_secondary, field_shortname, filename)

                if not os.path.exists(path_to_file): #Create file if it doesn't exist
                    comp_secondary.main(datdir_primary=datdir_primary, datdir_secondary=datdir_secondary, monthstr=monthstr, yearstr=yearstr, field_name=field_name, rho_ref=rho_ref, nu_E=nu_E, time_ave_type=time_ave_type)

                month += 1
                
            year += 1
            month = 1
        
        if year == int(final_year):

            while month <= int(final_month):

                monthstr = get_monthstr(month)
                date_string = final_year + '-' + monthstr

                filename = field_nc_string + date_string + '.nc'
                path_to_file = os.path.join(datdir_secondary, field_shortname, filename)

                if not os.path.exists(path_to_file): #Create file if it doesn't exist
                    comp_secondary.main(datdir_primary=datdir_primary, datdir_secondary=datdir_secondary, monthstr=monthstr, yearstr=final_year, field_name=field_name, rho_ref=rho_ref, nu_E=nu_E, time_ave_type=time_ave_type)
                
                month += 1
                
    print("Done computing secondary data.")
        
##############################

def load_secondary_data_file(field_name, date_string, datdir_secondary, time_ave_type):
    
    """
    Loads DataSet for a given computed field at a given time.
    Note - only does monthly avgs at the moment - will fix in future (i.e., will add conditions for other time_ave_type)
    """
    
    if time_ave_type == 'monthly':
    
        field_shortname = get_monthly_shortname(get_field_variable(field_name))
        field_nc_string = get_monthly_nc_string(get_field_variable(field_name))

    filename = field_nc_string + date_string + '.nc'
    path_to_file = join(datdir_secondary, field_shortname, filename)
    
    try: #Try to load DataSet
        computed_ds = xr.open_mfdataset(path_to_file, engine="scipy")
        print("Loaded DataSet.")
        return computed_ds
    
    except:
        print("DataSet does not exist.")