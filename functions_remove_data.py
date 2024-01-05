"""
Function to remove primary data.

Rosalie Cormier, 2024
"""

import os
import glob

from os.path import join

from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string

##############################

def remove_primary_files(field_name, datdir_primary, time_ave_type, date_strings):

    if time_ave_type == 'monthly': #will add other options
        field_shortname = get_monthly_shortname(get_field_variable(field_name))
        field_nc_string = get_monthly_nc_string(get_field_variable(field_name))
        
    if os.path.exists(join(datdir_primary, field_shortname)): #To avoid errors, only remove files after confirming directory exists
        for date_string in date_strings: #Iterate over times
            pattern = join(datdir_primary, field_shortname, field_nc_string+date_string+r"*")
            for item in glob.iglob(pattern, recursive=True): #Delete the files
                os.remove(item)  
        try:
            os.rmdir(join(datdir_primary, field_shortname)) #Delete the directory, if empty
        except:
            return join(datdir_primary, field_shortname)