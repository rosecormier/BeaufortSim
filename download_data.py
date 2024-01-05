"""
For a specified ECCO field (scalar or vector) and date range, iterates over time and downloads the field data.
If "seasonal" time-average type is specified, saves seasonal averages of primary data.

Rosalie Cormier, 2023
"""

import os

from os.path import join

from functions_ecco_download import ecco_podaac_download
from functions_ecco_general import get_monthstr, get_month_end, load_ECCO_data_file
from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string, get_seasonal_shortname, get_seasonal_nc_string

##############################

def compute_temporal_mean(timeseries):
    
    """
    Compute temporal mean of a field.
    """ 
    
    mean = (timeseries[0]).copy() / len(timeseries)
    
    if len(timeseries) > 1:
        for i in range(1, len(timeseries)):
            mean = mean + (timeseries[i]).copy() / len(timeseries)
        
    return mean

##############################

def compute_seasonal_average(monthly_fields, datdir_primary, field_name, season_start_string, season_end_string, yearstr):
    
    """
    Compute seasonal average of field from monthly averages; save output to NetCDF file.
    """ 

    seasonal_avg_field = compute_temporal_mean(monthly_fields)
    
    outdir = join(datdir_primary, 'Seasonal')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    #File to save data to
    filename = '{}_{}-{}_{}.nc'.format(field_name, season_start_string, season_end_string, yearstr)
                    
    if not os.path.exists(join(outdir, filename)): #Only compute if file doesn't already exist
        seasonal_avg_field.to_netcdf(path=join(outdir, filename), engine="scipy")
    
    seasonal_avg_field.close()

##############################

def main(**kwargs):
    
    if kwargs:
        
        field_name = kwargs.get('field_name')
        initial_month, initial_year = kwargs.get('initial_month'), kwargs.get('initial_year')
        final_month, final_year = kwargs.get('final_month'), kwargs.get('final_year')
        time_ave_type = kwargs.get('time_ave_type')
        datdir_primary = kwargs.get('datdir_primary')
        time_kwargs = kwargs.get('time_kwargs')
        
        month, year = int(initial_month), int(initial_year)
    
        if time_ave_type == "monthly":
        
            field_shortname = get_monthly_shortname(get_field_variable(field_name))
            field_nc_string = get_monthly_nc_string(get_field_variable(field_name))
    
            while year < int(final_year):
                yearstr = str(year)
                while month <= 12:
                    monthstr = get_monthstr(month)
                    endmonth = get_month_end(monthstr, yearstr)
                    StartDate, EndDate = yearstr+"-"+monthstr+"-02", yearstr+"-"+monthstr+"-"+endmonth
                    ecco_podaac_download(ShortName=field_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir_primary, n_workers=6, force_redownload=False)
                    month += 1
                year += 1
                month = 1
        
            if year == int(final_year):
                while month <= int(final_month):
                    monthstr = get_monthstr(month)
                    endmonth = get_month_end(monthstr, final_year)
                    StartDate, EndDate = final_year+"-"+monthstr+"-02", final_year+"-"+monthstr+"-"+endmonth
                    ecco_podaac_download(ShortName=field_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir_primary, n_workers=6, force_redownload=False)
                    month += 1
                    
            print("Done downloading ECCO data.")
        
        elif time_ave_type == "seasonal":
        
            field_shortname = get_monthly_shortname(get_field_variable(field_name))
            field_nc_string = get_monthly_nc_string(get_field_variable(field_name))
            
            season_start_string, season_end_string = time_kwargs[0], time_kwargs[1]
            season_start, season_end = int(season_start_string), int(season_end_string)
            
            if season_start < season_end:
                
                while year <= int(final_year):
                    
                    yearstr = str(year)
                    monthly_fields = []
                    
                    while month <= season_end:
                        monthstr = get_monthstr(month)
                        endmonth = get_month_end(monthstr, yearstr)
                        date_string = yearstr + "-" + monthstr
                        StartDate, EndDate = date_string + "-02", date_string + "-" + endmonth
                        ecco_podaac_download(ShortName=field_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir_primary, n_workers=6, force_redownload=False)
                        ds_month = load_ECCO_data_file(field_name, date_string, datdir_primary, 'seasonal') #Load the DataSet
                        monthly_fields.append(ds_month)
                        month += 1
                        
                    #After loading data for all months in the season, compute seasonal average
                    compute_seasonal_average(monthly_fields, datdir_primary, field_name, season_start_string, season_end_string, yearstr)
                        
                    year += 1
                    month = season_start

            elif season_start > season_end:
                
                while year <= int(final_year):
                    
                    monthly_fields = []
                    
                    while (season_start <= month) or (month <= season_end):
                        
                        yearstr = str(year)
                        monthstr = get_monthstr(month)
                        endmonth = get_month_end(monthstr, yearstr)
                        date_string = yearstr + "-" + monthstr
                        StartDate, EndDate = date_string + "-02", date_string + "-" + endmonth
                        ecco_podaac_download(ShortName=field_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir_primary, n_workers=6, force_redownload=False)
                        ds_month = load_ECCO_data_file(field_name, date_string, datdir_primary, 'seasonal') #Load the DataSet
                        monthly_fields.append(ds_month)
                        if month == 12:
                            year += 1
                            month = 1
                        else:
                            month += 1
                    
                    #After loading data for all months in the season, compute seasonal average
                    compute_seasonal_average(monthly_fields, datdir_primary, field_name, season_start_string, season_end_string, '{}-{}'.format(str(year-1), yearstr))
                        
                    month = season_start
    
            print("Done downloading ECCO data and saving seasonal averages.")
            
##############################

if __name__ == "__main__":
    main()