"""
For a specified ECCO field and date range, iterates over time and:
    -Checks whether ECCO-field data is already downloaded; and
        -If not, downloads the data.

Rosalie Cormier, 2023
"""

import os

from os.path import join

from functions_ecco_download import ecco_podaac_download
from functions_ecco_general import get_monthstr, get_month_end
from NEW_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string

##############################

def main(**kwargs):
    
    if kwargs:
        
        field_name = kwargs.get('field_name')
        initial_month, initial_year = kwargs.get('initial_month'), kwargs.get('initial_year')
        final_month, final_year = kwargs.get('final_month'), kwargs.get('final_year')
        time_ave_type = kwargs.get('time_ave_type')
        datdir_primary = kwargs.get('datdir_primary')
    
        month, year = int(initial_month), int(initial_year)
    
        if time_ave_type == "monthly": #This is the only case supported now; will modify to add others
        
            field_shortname = get_monthly_shortname(get_field_variable(field_name))
            field_nc_string = get_monthly_nc_string(get_field_variable(field_name))
    
            while year < int(final_year):
            
                yearstr = str(year)
            
                while month <= 12:
                
                    monthstr = get_monthstr(month-1) #The indexing is weird for months
                   
                    endmonth = get_month_end(monthstr, yearstr)
                    StartDate, EndDate = yearstr+"-"+monthstr+"-02", yearstr+"-"+monthstr+"-"+endmonth
                    ecco_podaac_download(ShortName=field_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir_primary, n_workers=6, force_redownload=False)
                        
                    month += 1
                
                year += 1
                month = 1
        
            if year == int(final_year):

                while month <= int(final_month):

                    monthstr = get_monthstr(month-1) #The indexing is weird for months
                   
                    endmonth = get_month_end(monthstr, final_year)
                    StartDate, EndDate = final_year+"-"+monthstr+"-02", final_year+"-"+monthstr+"-"+endmonth
                    ecco_podaac_download(ShortName=field_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir_primary, n_workers=6, force_redownload=False)

                    month += 1
    
    print("Done downloading ECCO data.")
    
##############################

if __name__ == "__main__":
    main()