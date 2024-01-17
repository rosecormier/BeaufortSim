"""
For a specified ECCO field (scalar or vector) and date range, iterates over time and downloads the field data.
If "seasonal" time-average type is specified, saves seasonal averages of primary data.

Rosalie Cormier, 2023
"""

import os
import xarray as xr

from os.path import join

import load_data_files

from functions_ecco_download import ecco_podaac_download
from functions_ecco_general import compute_temporal_mean, get_monthstr, get_month_end
from functions_field_variables import get_field_variable, get_vector_comps, get_monthly_shortname, get_monthly_nc_string, get_seasonal_shortname, get_seasonal_nc_string

##############################

def compute_seasonal_average(monthly_fields, datdir_primary, field_name, season_start_string, season_end_string, yearstr):
    
    """
    Compute seasonal average of field from monthly averages; save output to NetCDF file.
    """ 

    if monthly_fields[1] is None:
        seasonal_avg = compute_temporal_mean(monthly_fields[0])
    elif monthly_fields[1] is not None:
        seasonal_avg_0, seasonal_avg_1 = compute_temporal_mean(monthly_fields[0]).squeeze(), compute_temporal_mean(monthly_fields[1]).squeeze()
        seasonal_avg = xr.merge([seasonal_avg_0, seasonal_avg_1])
    
    outdir = join(datdir_primary, get_seasonal_shortname(get_field_variable(field_name)))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    #File to save data to
    filename = '{}{}-{}_{}.nc'.format(get_seasonal_nc_string(get_field_variable(field_name)), season_start_string, season_end_string, yearstr)
                    
    if not os.path.exists(join(outdir, filename)): #Only compute if file doesn't already exist
        seasonal_avg.to_netcdf(path=join(outdir, filename), engine="scipy")
    
    seasonal_avg.close()

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
        return None
        
    elif time_ave_type == "seasonal":
        
        field_shortname = get_monthly_shortname(get_field_variable(field_name))
        field_nc_string = get_monthly_nc_string(get_field_variable(field_name))
            
        season_start_string, season_end_string = time_kwargs[0], time_kwargs[1]
        season_start, season_end = int(season_start_string), int(season_end_string)
        
        #This list will be used to track the dates corresponding to the files we download
        all_seasons_all_date_strings = []
            
        if season_start < season_end:
                
            while year <= int(final_year):
                    
                yearstr = str(year)
                monthly_fields_0, monthly_fields_1 = None, None
                
                while month <= season_end:
                    monthstr = get_monthstr(month)
                    endmonth = get_month_end(monthstr, yearstr)
                    date_string = yearstr + "-" + monthstr
                    StartDate, EndDate = date_string + "-02", date_string + "-" + endmonth
                    ecco_podaac_download(ShortName=field_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir_primary, n_workers=6, force_redownload=False)
                    ds_month = load_data_files.main(field_name=field_name, date_string=date_string, datdir_primary=datdir_primary, time_ave_type='monthly') 
                        
                    try:
                        attributes = get_vector_comps(field_name) #For vector fields
                        if monthly_fields_0 is None:
                            monthly_fields_0, monthly_fields_1 = ds_month[attributes[0]], ds_month[attributes[1]]
                        elif monthly_fields_0 is not None:
                            monthly_fields_0 = xr.concat([monthly_fields_0, ds_month[attributes[0]]], 'time')
                            monthly_fields_1 = xr.concat([monthly_fields_1, ds_month[attributes[1]]], 'time')
                    except:
                        if monthly_fields_0 is None:
                            monthly_fields_0 = ds_month[get_field_variable(field_name)]
                        elif monthly_fields_0 is not None:
                            monthly_fields_0 = xr.concat([monthly_fields_0, ds_month[get_field_variable(field_name)]], 'time')
                    
                    all_seasons_all_date_strings.append(date_string)
                    
                    month += 1
                        
                #After loading data for all months in the season, compute seasonal average
                compute_seasonal_average([monthly_fields_0, monthly_fields_1], datdir_primary, field_name, season_start_string, season_end_string, yearstr)
                        
                year += 1
                month = season_start

        elif season_start > season_end:
            
            monthly_fields_0, monthly_fields_1 = None, None
            
            while year <= int(final_year):
                
                if (year != int(final_year) and season_start <= month) or (year != int(initial_year) and month <= season_end):

                    yearstr = str(year)
                    monthstr = get_monthstr(month)
                    endmonth = get_month_end(monthstr, yearstr)
                    date_string = yearstr + "-" + monthstr
                    StartDate, EndDate = date_string + "-02", date_string + "-" + endmonth
                    ecco_podaac_download(ShortName=field_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir_primary, n_workers=6, force_redownload=False)
                    ds_month = load_data_files.main(field_name=field_name, date_string=date_string, datdir_primary=datdir_primary, time_ave_type='monthly') 

                    try:
                        attributes = get_vector_comps(field_name) #For vector fields
                        if monthly_fields_0 is None:
                            monthly_fields_0, monthly_fields_1 = ds_month[attributes[0]], ds_month[attributes[1]]
                        elif monthly_fields_0 is not None:
                            monthly_fields_0 = xr.concat([monthly_fields_0, ds_month[attributes[0]]], 'time')
                            monthly_fields_1 = xr.concat([monthly_fields_1, ds_month[attributes[1]]], 'time')
                    except:
                        if monthly_fields_0 is None:
                            monthly_fields_0 = ds_month[get_field_variable(field_name)]
                        elif monthly_fields_0 is not None:
                            monthly_fields_0 = xr.concat([monthly_fields_0, ds_month[get_field_variable(field_name)]], 'time')
                        
                    all_seasons_all_date_strings.append(date_string)
                        
                if month == season_end and monthly_fields_0 is not None:
                        
                    #After loading data for all months in the season, compute seasonal average
                    compute_seasonal_average([monthly_fields_0, monthly_fields_1], datdir_primary, field_name, season_start_string, season_end_string, '{}-{}'.format(str(year-1), yearstr))
                        
                    month = season_start
                    monthly_fields_0, monthly_fields_1 = None, None
                    
                elif month == 12:
                    year += 1
                    month = 1   
                else:
                    month += 1
    
        print("Done downloading ECCO data and saving seasonal averages.")
        return all_seasons_all_date_strings
            
##############################

if __name__ == "__main__":
    main()