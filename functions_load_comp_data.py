"""
Loads a specified NetCDF file containing a DataSet computed from ECCO data.

Rosalie Cormier, 2023
"""

import os
import xarray as xr

from os.path import join

from functions_ecco_general import load_dataset

import compute_monthly_avgs
import save_annual_avgs
import save_seasonal_avgs

##############################

def load_comp_file(monthly_file, latmin, latmax, lonmin, lonmax, year, datdir, compdatdir):
    
    """
    Checks that a computed (monthly) file exists, and creates it if it doesn't, then loads DataSet.
    """
    
    if os.path.exists(monthly_file): #Look for the file
        ds_month = xr.open_mfdataset(monthly_file, engine="scipy") #Load DataSet

    else: #If it doesn't exist, compute it (for the year)
                                
        compute_monthly_avgs.main(latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax, startyr=year, years=1, datdir=datdir, outdir=compdatdir)
        
        ds_month = load_dataset(monthly_file) #Load DataSet
        
    return ds_month

##############################

def load_annual_scalar_ds(yearlydatdir, scalar_attr, year, datdir, scalarECCO, usecompdata=True):
    
    """
    Checks that an annual datafile exists, and creates it if it doesn't, then loads DataSet.
    """
    
    yearstr = str(year)
    
    scalar_annual_file = join(yearlydatdir, "avg_"+scalar_attr+"_"+yearstr+".nc")
    
    if not os.path.exists(scalar_annual_file): #If it doesn't exist, compute it
                        
        if scalarECCO: #If variable comes from ECCO directly
            usecompdata = False
                            
        save_annual_avgs.main(years=[year], field=scalar_attr, datdir=datdir, usecompdata=usecompdata, outdir=yearlydatdir)
                    
    ds_scalar_year = xr.open_mfdataset(scalar_annual_file, engine="scipy")
    ds_scalar_year.load() #Load DataSet
    
    return ds_scalar_year

##############################

def load_seasonal_scalar_ds(seasonaldatdir, scalar_attr, season_start, year, season_end, endyearstr, scalarECCO):
    
    """
    Checks that a seasonal datafile exists, and creates it if it doesn't, then loads DataSet.
    """
    
    yearstr = str(year)
    
    scalar_seas_file = join(seasonaldatdir, "avg_"+scalar_attr+"_"+season_start+yearstr+"-"+season_end+endyearstr+".nc")
                    
    if not os.path.exists(scalar_seas_file): #If it doesn't exist, compute it
                        
        if scalarECCO: #If variable comes from ECCO directly
            datdirshort, usecompdata = 'Downloads', False
                            
        elif not scalarECCO:
            datdirshort, usecompdata = 'computed_monthly', True
                        
        save_seasonal_avgs.main(field=scalar_attr, years=[year], start_month=season_start, end_month=season_end, usecompdata=usecompdata, datdir=datdirshort, outdir=seasonaldatdir)
                    
    ds_scalar_seas = xr.open_mfdataset(scalar_seas_file, engine="scipy")
    ds_scalar_seas.load() #Load DataSet
    
    return ds_scalar_seas