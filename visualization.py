"""
Master script for Beaufort Sim visualization.

Rosalie Cormier, 2023
"""

##############################

#IMPORTS AND SETTINGS

import sys
import os
import argparse
import ecco_v4_py as ecco
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

from os.path import expanduser, join

from functions_ecco_general import load_grid, get_monthstr, load_dataset, ds_to_field, rotate_vector, get_vector_partner, ecco_resample, get_season_months_and_years, get_scalar_in_xy
from functions_visualization import plot_pcolormesh_k_plane, plot_pcm_quiver_k_plane
from functions_field_variables import get_field_vars, get_variable_str
from functions_load_comp_data import check_for_ecco_file, load_comp_file, load_annual_scalar_ds, load_annual_vector_ds, load_seasonal_scalar_ds, load_seasonal_vector_ds
from functions_geostrophy import rotate_comp_vector, comp_geos_metric

#The following are scripts that are imported as modules but may be run within this script

import compute_monthly_avgs
import save_annual_avgs
import save_seasonal_avgs

##############################

def get_parser():

    parser = argparse.ArgumentParser(description="Plot various attributes in Beaufort Gyre", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    #Spatial bounds

    parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, default=[70.5, 80.0])
    parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, default=[-155.0, -120.0])
    parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, default=0.25)
    parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 1])

    #Temporal bounds

    parser.add_argument("start", type=int, help="Start year") #This argument is required
    parser.add_argument("--years", type=int, help="Total number of years to plot", default=1)
    parser.add_argument("--seasonal", dest='seasonal', help="Whether to plot specific seasons", default=False, action='store_true')
    parser.add_argument("--seasonmonths", type=str, help="Start and end months of season", nargs=2, default=["01", "01"])

    #Attributes

    parser.add_argument("--scalar", type=str, help="Name of scalar attribute", default="ZETA")
    parser.add_argument("--scalarECCO", dest='scalarECCO', help="Whether scalar field comes from ECCO files", default=False, action='store_true')
    parser.add_argument("--vminmax", type=float, help="Minimum/maximum scalar values", nargs=2, default=[1, 1])
    parser.add_argument("--xvec", type=str, help="Name of vector attribute (x-comp)", default=None)
    parser.add_argument("--vectorECCO", dest='vectorECCO', help="Whether vector field comes from ECCO files", default=False, action='store_true')

    #Directories

    parser.add_argument("--datdir", type=str, help="Directory (rel. to home) with raw ECCO data", default="Downloads")
    parser.add_argument("--compdatdir", type=str, help="Directory (rel. to here) with computed monthly data", default="computed_monthly")
    parser.add_argument("--seasonaldatdir", type=str, help="Directory (rel. to here) with seasonal avgs", default="seasonal_averages")
    parser.add_argument("--yearlydatdir", type=str, help="Directory (rel. to here) with annual avgs", default="yearly_averages")
    parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", default="visualization")

    #Visualization
    
    parser.add_argument("--cmap", type=str, help="MPL colormap", default="viridis_r")
    parser.add_argument("--plane", type=str, help="Plane on which to plot (k/lon/lat)", default="k")
    
    return parser

##############################

def load_ECCO_ds_scalar_mo(scalar_dir, scalar_attr, monthstr, year, datdir, ds_grid):

    """
    Used in this file to load a monthly scalar ECCO DataSet.
    If scalar is w, interpolates vertical velocity component along z-axis.
    """
                            
    #Ensure file is downloaded
    scalar_file = check_for_ecco_file(scalar_dir, scalar_attr, monthstr, year, datdir)
    ds_scalar_mo = load_dataset(scalar_file) #Load data
                            
    if scalar_attr == "WVEL": #If w, interpolate vertically
        XGCM_grid = ecco.get_llc_grid(ds_grid)
        ds_scalar_mo[scalar_attr] = XGCM_grid.interp(ds_scalar_mo.WVEL, axis="Z")
        
    return ds_scalar_mo

##############################

def load_ECCO_vector_mo_components(vector_dir, monthstr, year, xvec_attr, yvec_attr, datdir, ds_grid, k):
    
    """
    Used in this file to load a monthly vector ECCO DataSet and return east/northward components at depth k.
    """
    
    variable = xvec_attr + yvec_attr
    
    vector_file = check_for_ecco_file(vector_dir, variable, monthstr, year, datdir)
    ds_vector_mo = load_dataset(vector_file) #Load data
    
    vecE, vecN = rotate_vector(ds_grid, ds_vector_mo, xvec_attr, yvec_attr) #Rotate vector
    vecE, vecN = vecE.isel(k=k).squeeze(), vecN.isel(k=k).squeeze() #Squeeze along time axis

    return vecE, vecN
    
##############################

def load_comp_scalar_ds_and_grid(scalar_dir, scalar_monthly_nc_str, monthstr, year, lats_lons, datdir, compdatdir, ds_grid, scalar_attr):
    
    """
    Used in this file to load a monthly computed scalar DataSet and return it alongside grid DataSet.
    """
    
    yearstr = str(year)
    
    scalar_file = join(scalar_dir, scalar_monthly_nc_str+yearstr+"-"+monthstr+".nc")
                            
    #Ensure file exists and load data
    ds_scalar_mo = load_comp_file(scalar_file, lats_lons, year, datdir, compdatdir)
    ds_grid = get_scalar_in_xy(ds_grid, ds_scalar_mo, scalar_attr) 
        
    return ds_scalar_mo, ds_grid

##############################

def load_comp_vector_ds_and_grid(vector_dir, vector_monthly_nc_str, xvec_attr, yvec_attr, monthstr, year, lats_lons, datdir, compdatdir, ds_grid, k, years):
    
    """
    Used in this file to load a monthly computed vector DataSet and return east/northward components at depth k.
    """
    
    yearstr = str(year)
    latmin, latmax, lonmin, lonmax = lats_lons
    
    vector_file = join(vector_dir, vector_monthly_nc_str+yearstr+"-"+monthstr+".nc")
                            
    if not os.path.exists(vector_file): #If it doesn't exist, compute it
        compute_monthly_avgs.main(latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax, startyr=year, years=years, datdir=datdir, outdir=compdatdir)
                                
    ds_vector_mo = xr.open_mfdataset(vector_file, engine="scipy") #Load monthly vector file
                            
    vecE, vecN = rotate_comp_vector(ds_grid, ds_vector_mo[xvec_attr], ds_vector_mo[yvec_attr], k)
    
    return vecE, vecN

##############################
            
def main():
    
    parser = get_parser()
    args = parser.parse_args()
    config = vars(args)

    #Spatial bounds

    latmin, latmax, lonmin, lonmax = config['lats'][0], config['lats'][1], config['lons'][0], config['lons'][1]
    kmin, kmax = config['kvals'][0], config['kvals'][1]
    
    resolution = config['res']

    #Temporal bounds

    startyr, years, seasonal = config['start'], config['years'], config['seasonal']

    if seasonal:
        season_start, season_end = config['seasonmonths']

    #Attributes

    scalar_attr = config['scalar']
    vmin, vmax = config['vminmax'][0], config['vminmax'][1]

    include_vector_field = False
    xvec_attr = config['xvec']

    if xvec_attr is not None:

        include_vector_field = True
        yvec_attr = get_vector_partner(xvec_attr)

        variables_str = get_variable_str(xvec_attr+yvec_attr) + "_" + get_variable_str(scalar_attr)

    elif xvec_attr is None:
        variables_str = get_variable_str(scalar_attr)

    scalarECCO, vectorECCO = config['scalarECCO'], config['vectorECCO']

    #Directories

    homedir = expanduser('~')
    sys.path.append(join(homedir, 'ECCOv4-py'))
    datdir = join(homedir, config['datdir'], 'ECCO_V4r4_PODAAC')

    compdatdir = join(".", config['compdatdir'])

    outdir = join(".", config['outdir'])

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if seasonal:
        subdirs = ["seasonal", "interannual"]
        seasonaldatdir = join(".", config['seasonaldatdir'])

    elif not seasonal:
        subdirs = ["monthly", "yearly"]
        yearlydatdir = join(".", config['yearlydatdir'])

    for subdir in subdirs:
        if not os.path.exists(join(outdir, subdir)):
            os.makedirs(join(outdir, subdir))

    #Visualization
            
    cmap = config['cmap']
    plane = config['plane']
    
    if plane == 'k':
        lats_lons = [latmin, latmax, lonmin, lonmax]
        
    elif plane == 'lon':
        lats_depths = [latmin, latmax, kmin, kmax]
        
    elif plane == 'lat':
        lons_depths = [lonmin, lonmax, kmin, kmax]

    ##############################

    #GET FILE NAMES

    scalar_monthly_shortname, scalar_monthly_nc_str = get_field_vars(scalar_attr)

    if scalarECCO:
        scalar_dir = join(datdir, scalar_monthly_shortname)

    elif not scalarECCO:
        scalar_dir = join(compdatdir, scalar_monthly_shortname)

    if include_vector_field:

        vector_monthly_shortname, vector_monthly_nc_str = get_field_vars(xvec_attr+yvec_attr)

        if vectorECCO:
            vector_dir = join(datdir, vector_monthly_shortname)

        elif not vectorECCO:
            vector_dir = join(compdatdir, vector_monthly_shortname)
 
    ds_grid = load_grid(datdir) #Load grid 
    
    ##############################

    #CASE WHERE ONLY A SCALAR IS PROVIDED

    if not include_vector_field:

        for k in range(kmin, kmax + 1): #Iterate over specified depths

            if not seasonal: #Case where we plot every month

                for i in range(years): #Iterate over specified years

                    year = startyr + i
                    yearstr = str(year)

                    Ro_l_list, OW_list = [], [] #Only used if scalar attribute is ZETA

                    for m in range(12): #Iterate over months
                        
                        monthstr = get_monthstr(m)
                        
                        #Load monthly DataSet
                        
                        if scalarECCO: #If variable comes from ECCO directly
                            ds_scalar_mo = load_ECCO_ds_scalar_mo(scalar_dir, scalar_attr, monthstr, year, config['datdir'], ds_grid)

                        elif not scalarECCO:
                            ds_scalar_mo, ds_grid = load_comp_scalar_ds_and_grid(scalar_dir, scalar_monthly_nc_str, monthstr, year, lats_lons, config['datdir'], compdatdir, ds_grid, scalar_attr)
                            
                        #File to save monthly plot to
                        outfile = join(outdir, 'monthly', '{}_k{}_{}{}.pdf'.format(variables_str, str(k), monthstr, yearstr))
                            
                        #Plot monthly data    
                        Ro_l_list, OW_list = plot_pcolormesh_k_plane(ds_grid, [ds_scalar_mo], k, scalar_attr, resolution, cmap, monthstr+"-"+yearstr, vmin, vmax, outfile, lats_lons, datdir, Ro_l_list, OW_list, yearstr, monthstr=monthstr, datdirname=config['datdir'], outdir=outdir)
       
                    #Get annually-averaged data, with interpolation if necessary
                    ds_scalar_year = load_annual_scalar_ds(yearlydatdir, scalar_attr, year, config['datdir'], ds_grid, scalarECCO)
                    
                    #File to save annual plot to
                    outfile = join(outdir, 'yearly', '{}_k{}_{}.pdf'.format(variables_str, str(k), yearstr))
                    
                    #Plot annual data
                    plot_pcolormesh_k_plane(ds_grid, [ds_scalar_year], k, scalar_attr, resolution, cmap, yearstr, vmin, vmax, outfile, lats_lons, datdir, Ro_l_list, OW_list, yearstr, outdir=outdir, annual=True)

            elif seasonal: #Case where we plot one season per year

                season_months, season_years = get_season_months_and_years(season_start, season_end)
                seas_monthstr = season_months[0] + "-" + season_months[-1] #For titles

                data_seasons = []
                Ro_l_list, OW_list = [], []

                for i in range(years): #Iterate over specified years

                    year = startyr + i
                    yearstr = str(year)
                    endyearstr = str(year + season_years[-1])
                    
                    #Get seasonally-averaged data
                    ds_scalar_seas = load_seasonal_scalar_ds(seasonaldatdir, scalar_attr, season_start, year, season_end, endyearstr, scalarECCO)
                    
                    if scalar_attr == "WVEL": #If w, interpolate vertically
                        ds_scalar_seas[scalar_attr] = XGCM_grid.interp(ds_scalar_seas.WVEL, axis="Z")

                    seas_yearstr = yearstr + "-" + str(year + season_years[-1]) #For titles
                    
                    outfile = join(outdir, 'seasonal', '{}_k{}_{}_{}.pdf'.format(variables_str, str(k), seas_monthstr, seas_yearstr))
                    
                    Ro_l_list, OW_list, data_seasons, lon_centers, lat_centers = plot_pcolormesh_k_plane(ds_grid, [ds_scalar_seas], k, scalar_attr, resolution, cmap, '{}, {}'.format(seas_monthstr, seas_yearstr), vmin, vmax, outfile, lats_lons, datdir, Ro_l_list, OW_list, yearstr, year=year, outdir=outdir, season_start=season_start, season_end=season_end, endyearstr=endyearstr, datdirname=config['datdir'], seasonal=True, seasonaldatdir=seasonaldatdir, data_seasons=data_seasons)
                
                if years != 1: #If there is more than one season to average over
                
                    datestr = '{}, {}'.format(seas_monthstr, seas_yearstr)
                    outfile = join(outdir, 'interannual', '{}_k{}_{}_{}.pdf'.format(variables_str, str(k), seas_monthstr, seas_yearstr))
          
                    #Plot average over all seasons, reusing seasonal lat/lon_centers
                    plot_pcolormesh_k_plane(ds_grid, data_seasons, k, scalar_attr, resolution, cmap, datestr, vmin, vmax, outfile, lats_lons, datdir, Ro_l_list, OW_list, yearstr, outdir=outdir, seas_monthstr=seas_monthstr, seas_yearstr=seas_yearstr, season_years=season_years, years=years, startyr=startyr, multiple_seas=True, lon_centers=lon_centers, lat_centers=lat_centers) 
                
    ##############################

    #CASE WHERE VECTOR AND SCALAR ARE PROVIDED

    elif include_vector_field:

        for k in range(kmin, kmax + 1): #Iterate over specified depths

            if not seasonal: #Case where we plot every month

                for i in range(years): #Iterate over specified years

                    year = startyr + i
                    yearstr = str(year)
                    
                    Ro_l_list, OW_list = [], [] #Only used if scalar attribute is ZETA

                    for m in range(12): #Iterate over months
                        
                        monthstr = get_monthstr(m)

                        #Load monthly scalar DataSet and vector components
                        
                        if scalarECCO: #If variable comes from ECCO directly
                            ds_scalar_mo = load_ECCO_ds_scalar_mo(scalar_dir, scalar_attr, monthstr, year, config['datdir'], ds_grid)

                        elif not scalarECCO:
                            ds_scalar_mo, ds_grid = load_comp_scalar_ds_and_grid(scalar_dir, scalar_monthly_nc_str, monthstr, year, lats_lons, config['datdir'], compdatdir, ds_grid, scalar_attr)

                        if vectorECCO: #If variable comes from ECCO directly
                            vecE, vecN = load_ECCO_vector_mo_components(vector_dir, monthstr, year, xvec_attr, yvec_attr, config['datdir'], ds_grid, k)
                            
                        elif not vectorECCO:
                            vecE, vecN = load_comp_vector_ds_and_grid(vector_dir, vector_monthly_nc_str, xvec_attr, yvec_attr, monthstr, year, lats_lons, config['datdir'], compdatdir, ds_grid, k, years)
                            
                        #File to save monthly plot to
                        outfile = join(outdir, 'monthly', '{}_k{}_{}{}.pdf'.format(variables_str, str(k), monthstr, yearstr))
                        
                        #Plot monthly data
                        Ro_l_list, OW_list = plot_pcm_quiver_k_plane(ds_grid, [ds_scalar_mo], k, scalar_attr, xvec_attr, vecE, vecN, resolution, cmap, monthstr+"-"+yearstr, vmin, vmax, outfile, lats_lons, datdir, Ro_l_list, OW_list, yearstr, outdir=outdir, monthstr=monthstr, datdirname=config['datdir'])
                        
                    #Get annually-averaged scalar data, with interpolation if necessary
                    ds_scalar_year = load_annual_scalar_ds(yearlydatdir, scalar_attr, year, config['datdir'], ds_grid, scalarECCO)
                   
                    #Get annually-averaged vector components
                    vecE, vecN = load_annual_vector_ds(yearlydatdir, xvec_attr, yvec_attr, year, config['datdir'], ds_grid, k, vectorECCO, compdatdir=compdatdir)
                    
                    #File to save annual plot to
                    outfile = join(outdir, 'yearly', '{}_k{}_{}.pdf'.format(variables_str, str(k), yearstr))
                    
                    #Plot annual average
                    plot_pcm_quiver_k_plane(ds_grid, [ds_scalar_year], k, scalar_attr, xvec_attr, vecE, vecN, resolution, cmap, yearstr, vmin, vmax, outfile, lats_lons, datdir, Ro_l_list, OW_list, yearstr, outdir=outdir, annual=True)
                    
            elif seasonal: #Case where we plot one season per year
                    
                season_months, season_years = get_season_months_and_years(season_start, season_end)
                seas_monthstr = season_months[0] + "-" + season_months[-1] #For titles

                scalar_data_seasons = []
                vecE_data_seasons, vecN_data_seasons = [], []
                
                Ro_l_list, OW_list = [], [] #Only used if scalar is ZETA
                
                for i in range(years): #Iterate over specified years
                
                    year = startyr + i
                    yearstr, endyearstr = str(year), str(year + season_years[-1])

                    #Get seasonally-averaged scalar data
                    ds_scalar_seas = load_seasonal_scalar_ds(seasonaldatdir, scalar_attr, season_start, year, season_end, endyearstr, scalarECCO)
                    
                    if scalar_attr == "WVEL": #If w, interpolate vertically
                        ds_scalar_seas[scalar_attr] = XGCM_grid.interp(ds_scalar_seas.WVEL, axis="Z")
                    
                    #Get seasonally-averaged vector data
                    vecE, vecN = load_seasonal_vector_ds(seasonaldatdir, xvec_attr, yvec_attr, season_start, year, season_end, endyearstr, vectorECCO, ds_grid, k, compdatdir=compdatdir)
                  
                    seas_yearstr = yearstr + "-" + str(year + season_years[-1]) #For titles
                    
                    outfile = join(outdir, 'seasonal', '{}_k{}_{}_{}.pdf'.format(variables_str, str(k), seas_monthstr, seas_yearstr))
                    
                    #Plot seasonal average
                    Ro_l_list, OW_list, scalar_data_seasons, vecE_data_seasons, vecN_data_seasons, lon_centers, lat_centers = plot_pcm_quiver_k_plane(ds_grid, [ds_scalar_seas], k, scalar_attr, xvec_attr, vecE, vecN, resolution, cmap, '{}, {}'.format(seas_monthstr, seas_yearstr), vmin, vmax, outfile, lats_lons, datdir, Ro_l_list, OW_list, yearstr, outdir=outdir, seas_monthstr=seas_monthstr, seas_yearstr=seas_yearstr, season_start=season_start, season_end=season_end, endyearstr=endyearstr, seasonal=True, seasonaldatdir=seasonaldatdir, scalar_data_seasons=scalar_data_seasons, vecE_data_seasons=vecE_data_seasons, vecN_data_seasons=vecN_data_seasons, datdirname=config['datdir'])
                    
                if years != 1: #If there is more than one season to average over     
                    
                    seas_yearstr = str(startyr) + "-" + str(startyr + (years-1) + season_years[-1]) #For titles
                    outfile = join(outdir, 'interannual', '{}_k{}_{}_{}.pdf'.format(variables_str, str(k), seas_monthstr, seas_yearstr))
                    
                    #Plot interannual average, reusing seasonal lat/lon_centers
                    plot_pcm_quiver_k_plane(ds_grid, [ds_scalar_seas], k, scalar_attr, xvec_attr, vecE, vecN, resolution, cmap, '{}, {}'.format(seas_monthstr, seas_yearstr), vmin, vmax, outfile, lats_lons, datdir, Ro_l_list, OW_list, yearstr, outdir=outdir, seas_monthstr=None, seas_yearstr=None, multiple_seas=True, years=years, startyr=startyr, scalar_data_seasons=scalar_data_seasons, vecE_data_seasons=vecE_data_seasons, vecN_data_seasons=vecN_data_seasons, lon_centers=lon_centers, lat_centers=lat_centers)
            
##############################

if __name__ == "__main__":
    main()