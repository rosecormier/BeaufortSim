"""
Loads a specified NetCDF file containing a DataSet computed from ECCO data
(or a DataSet containing an ECCO vector).

Rosalie Cormier, 2023
"""

import os
import xarray as xr
import ecco_v4_py as ecco

from os.path import join

from functions_ecco_general import load_dataset, rotate_vector
from functions_field_variables import get_field_vars
from functions_geostrophy import rotate_comp_vector

import compute_monthly_avgs
import download_new_data
import save_annual_avgs
import save_seasonal_avgs

##############################

def check_for_ecco_file(variable_dir, variable, monthstr, year, datdir):
    
    """
    Checks that an ECCO file exists, and downloads it if it doesn't.
    """
    
    yearstr = str(year)
    
    variable_monthly_nc_str = get_field_vars(variable)[1]
    monthly_file = join(variable_dir, variable_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
     
    if not os.path.exists(monthly_file): #If the file doesn't exist, download data for that year
        download_new_data.main(startmo="01", startyr=year, months=12, variables=[variable], datdir=datdir)
        
    return monthly_file

##############################

def load_comp_file(monthly_file, lats_lons, year, datdir, compdatdir):
    
    """
    Checks that a computed (monthly) file exists, and creates it if it doesn't, then loads DataSet.
    """
    
    latmin, latmax, lonmin, lonmax = lats_lons
    
    if os.path.exists(monthly_file): #Look for the file
        ds_month = xr.open_mfdataset(monthly_file, engine="scipy") #Load DataSet

    else: #If it doesn't exist, compute it (for the year)
                                
        compute_monthly_avgs.main(latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax, startyr=year, years=1, datdir=datdir, outdir=compdatdir)
        
        ds_month = load_dataset(monthly_file) #Load DataSet
        
    return ds_month

##############################

def load_annual_scalar_ds(yearlydatdir, scalar_attr, year, datdir, ds_grid, scalarECCO):
    
    """
    Checks that an annual datafile exists, and creates it if it doesn't, then loads DataSet.
    """
    
    yearstr = str(year)
    
    scalar_annual_file = join(yearlydatdir, "avg_"+scalar_attr+"_"+yearstr+".nc")
    
    if not os.path.exists(scalar_annual_file): #If it doesn't exist, compute it
                        
        if scalarECCO: #If variable comes from ECCO directly
            usecompdata = False
        elif not scalarECCO:
            usecompdata = True
                            
        save_annual_avgs.main(years=[year], field=scalar_attr, datdir=datdir, usecompdata=usecompdata, outdir=yearlydatdir)
                    
    ds_scalar_year = xr.open_mfdataset(scalar_annual_file, engine="scipy")
    ds_scalar_year.load() #Load DataSet
    
    if scalar_attr == "WVEL": #If w, interpolate vertically  
        XGCM_grid = ecco.get_llc_grid(ds_grid)
        ds_scalar_year[scalar_attr] = XGCM_grid.interp(ds_scalar_year.WVEL, axis="Z")
    
    return ds_scalar_year

##############################

def load_annual_vector_ds(yearlydatdir, xvec_attr, yvec_attr, year, datdir, ds_grid, k, vectorECCO, compdatdir=None):
    
    """
    Checks that an annual datafile exists, and creates it if it doesn't, then loads DataSet and gets vector components.
    Can be used for ECCO or computed data.
    """
    
    yearstr = str(year)
    
    vector_annual_file = join(yearlydatdir, "avg_"+xvec_attr+yvec_attr+"_"+yearstr+".nc")

    if not os.path.exists(vector_annual_file): #If it doesn't exist, compute it
                        
        if vectorECCO: #If variable comes from ECCO directly
            usecompdata = False
        elif not vectorECCO:
            usecompdata = True
            datdir = compdatdir
                                    
        save_annual_avgs.main(years=[year], field=xvec_attr+yvec_attr, datdir=datdir, usecompdata=usecompdata, outdir=yearlydatdir)
                        
    ds_vector_year = xr.open_mfdataset(vector_annual_file, engine="scipy")
    ds_vector_year.load()
                    
    if vectorECCO: 
        vecE, vecN = rotate_vector(ds_grid, ds_vector_year, xvec_attr, yvec_attr)
        vecE, vecN = vecE.isel(k=k).squeeze(), vecN.isel(k=k).squeeze()
                        
    elif not vectorECCO:
        vecE, vecN = rotate_comp_vector(ds_grid, ds_vector_year[xvec_attr], ds_vector_year[yvec_attr], k)
        
    return vecE, vecN

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