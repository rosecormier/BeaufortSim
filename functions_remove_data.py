"""
Functions to:
    -Assemble a list of date strings to iterate over, particularly for file 
    removal;
    -Remove primary datafiles.

Rosalie Cormier, 2024
"""

import os
import glob

from os.path import join

from functions_ecco_general import get_monthstr
from functions_field_variables import get_field_variable, \
get_monthly_shortname, get_monthly_nc_string, get_seasonal_shortname, \
get_seasonal_nc_string

##############################

def get_date_strings(initial_month, initial_year, final_month, final_year, 
                     time_ave_type, time_kwargs):
    
    date_strings = []
    month, year = int(initial_month), int(initial_year)

    #Append every eligible date string to list 'date_strings'

    if time_ave_type == 'monthly':
    
        while year < int(final_year):
            while month <= 12:
                date_string = str(year) + '-' + get_monthstr(month)
                date_strings.append(date_string)
                month += 1
            year += 1
            month = 1

        if year == int(final_year):
            while month <= int(final_month):
                date_string = final_year + '-' + get_monthstr(month)
                date_strings.append(date_string)
                month += 1
            
    elif time_ave_type == 'seasonal':
        
        season_start, season_end = time_kwargs[0], time_kwargs[1]
        
        if int(season_start) < int(season_end):
            while year <= int(final_year):
                date_string = '{}-{}_{}'.format(season_start, season_end, 
                                                str(year))
                date_strings.append(date_string) 
                year += 1

        elif int(season_end) < int(season_start):
            while year < int(final_year):
                date_string = '{}-{}_{}-{}'.format(season_start, season_end, 
                                                   str(year), str(year+1))
                date_strings.append(date_string)
                year += 1
                
    return date_strings

##############################

def remove_primary_files(field_name, datdir_primary, date_strings):

    field_shortname = get_monthly_shortname(get_field_variable(field_name))
    field_nc_string = get_monthly_nc_string(get_field_variable(field_name))

    #To avoid errors, only remove files after confirming directory exists
    if os.path.exists(join(datdir_primary, field_shortname)):
        for date_string in date_strings: #Iterate over times
            pattern = join(datdir_primary, field_shortname, 
                           field_nc_string+date_string+r"*")
            for item in glob.iglob(pattern, recursive=True): 
                os.remove(item) #Delete the files
        try:
            #Delete the directory, if empty
            os.rmdir(join(datdir_primary, field_shortname)) 
        except:
            return join(datdir_primary, field_shortname)