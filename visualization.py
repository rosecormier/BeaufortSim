"""
Master script for Beaufort Sim visualization.

Rosalie Cormier, 2023

Note -- if vector = u_g, need to also plot our metric
    -- if scalar = vorticity, need to also plot W
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

from ecco_general import load_grid, get_monthstr, load_dataset, ds_to_field, comp_residuals, rotate_vector, get_vector_partner, ecco_resample, get_season_months_and_years, get_scalar_in_xy
from ecco_visualization import ArcCir_pcolormesh, ArcCir_contourf_quiver, ArcCir_contourf_quiver_grid
from ecco_field_variables import get_field_vars, get_variable_str

vir_nanmasked = plt.get_cmap('viridis_r').copy()
vir_nanmasked.set_bad('black')

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot various attributes in Beaufort Gyre",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#Spatial bounds

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, \
                    default=[70.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, \
                    default=[-180.0, -90.0])
parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, \
                    default=1.0)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 1])

#Temporal bounds

parser.add_argument("start", type=int, help="Start year") #This argument is required
parser.add_argument("--years", type=int, help="Total number of years to plot", default=1)
parser.add_argument("--seasonal", dest='seasonal', help="Whether to plot specific seasons", \
                    default=False, action='store_true')
parser.add_argument("--seasonmonths", type=str, help="Start and end months of season", nargs=2, \
                    default=["01", "01"])

#Attributes

parser.add_argument("--scalar", type=str, help="Name of scalar attribute", default="ZETA")
parser.add_argument("--scalarECCO", dest='scalarECCO', help="Whether scalar field comes from ECCO files", \
                    default=False, action='store_true')
parser.add_argument("--vminmax", type=float, help="Minimum/maximum scalar values", nargs=2, \
                    default=[1, 1])
parser.add_argument("--xvec", type=str, help="Name of vector attribute (x-comp)", nargs=1, default=None)
parser.add_argument("--vectorECCO", dest='vectorECCO', help="Whether vector field comes from ECCO files", \
                    default=False, action='store_true')

#Directories

parser.add_argument("--datdir", type=str, help="Directory (rel. to home) with raw ECCO data", \
                    default="Downloads")
parser.add_argument("--compdatdir", type=str, help="Directory (rel. to here) with computed monthly data", \
                    default="computed_monthly")
parser.add_argument("--seasonaldatdir", type=str, help="Directory (rel. to here) with seasonal avgs", \
                    default="seasonal_averages")
parser.add_argument("--yearlydatdir", type=str, help="Directory (rel. to here) with annual avgs", \
                    default="yearly_averages")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", \
                    default="visualization")

#Visualization
parser.add_argument("--cmap", type=str, help="MPL colormap", default="viridis")

args = parser.parse_args()
config = vars(args)

#Spatial bounds

latmin, latmax = config['lats'][0], config['lats'][1]
lonmin, lonmax = config['lons'][0], config['lons'][1]
lats_lons = [latmin, latmax, lonmin, lonmax]

resolution = config['res']
kmin, kmax = config['kvals'][0], config['kvals'][1]

#Temporal bounds

startyr = config['start']
years = config['years']
seasonal = config['seasonal']

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

user_home_dir = expanduser('~')
sys.path.append(join(user_home_dir, 'ECCOv4-py'))
datdir = join(user_home_dir, config['datdir'], 'ECCO_V4r4_PODAAC')
  
compdatdir = join(".", config['compdatdir'])

outdir = join(".", config['outdir'])

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
if seasonal:
    
    subdirs = ["seasonal"]
    
    seasonaldatdir = join(".", config['seasonaldatdir'])
    
elif not seasonal:
    
    subdirs = ["monthly", "yearly"]
    
    yearlydatdir = join(".", config['yearlydatdir'])

for subdir in subdirs:
    if not os.path.exists(join(outdir, subdir)):
        os.makedirs(join(outdir, subdir))
        
cmap = config['cmap']

f_mean = 1e-4 #"Typical" Coriolis parameter (1/s)
        
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
                
                for m in range(12): #Iterate over months
                    
                    monthstr = get_monthstr(m)

                    if scalarECCO:
                        
                        curr_scalar_file = join(scalar_dir, scalar_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
                        
                        #Load monthly scalar file into workspace
                        ds_scalar_mo = load_dataset(curr_scalar_file)
                    
                    elif not scalarECCO:
                        
                        curr_scalar_file = join(scalar_dir, scalar_monthly_nc_str+yearstr+"-"+monthstr+".nc")
                        
                        #Load monthly scalar file into workspace
                        ds_scalar_mo = xr.open_mfdataset(curr_scalar_file, engine="scipy")
                        
                        ds_grid = get_scalar_in_xy(ds_grid, ds_scalar_mo, scalar_attr)

                    #Convert scalar DataSet to useful field
                    lon_centers, lat_centers, lon_edges, lat_edges, scalar = ds_to_field(ds_grid, ds_scalar_mo.isel(k=k), scalar_attr, latmin, latmax, lonmin, lonmax, resolution)
                    
                    #Plot monthly scalar data
                    ArcCir_pcolormesh(ds_grid, k, [scalar], resolution, cmap, lon_centers, lat_centers, monthstr+"-"+yearstr, scalar_attr, scalar_bounds=[vmin, vmax], extend='both', outfile=join(outdir, 'monthly', '{}_k{}_{}{}.pdf'.format(variables_str, str(k), monthstr, yearstr)), lats_lons=lats_lons) 
                    
                    if scalar_attr == 'ZETA': #If vorticity, also normalize, and compute and plot Ro_l, W
                        
                        scalar *= 1 / f_mean #Normalize vorticity by typical f
                        
                        #Plot normalized vorticity
                        ArcCir_pcolormesh(ds_grid, k, [scalar], resolution, cmap, lon_centers, lat_centers, monthstr+"-"+yearstr, 'zetanorm', scalar_bounds=[vmin, vmax], extend='both', outfile=join(outdir, 'monthly', 'norm_{}_k{}_{}{}.pdf'.format(variables_str, str(k), monthstr, yearstr)), lats_lons=lats_lons)
                    
                #Get annually-averaged data
                
                scalar_annual_file = join(yearlydatdir, "avg_"+scalar_attr+"_"+yearstr+".nc")
                ds_scalar_year = xr.open_mfdataset(scalar_annual_file, engine="scipy")
                ds_scalar_year.load()

                #Convert scalar DataSet to useful field
                lon_centers, lat_centers, lon_edges, lat_edges, scalar_year = ds_to_field(ds_grid, ds_scalar_year.isel(k=k), scalar_attr, latmin, latmax, lonmin, lonmax, resolution)

                #Plot annual average
                ArcCir_pcolormesh(ds_grid, k, [scalar_year], resolution, cmap, lon_centers, lat_centers, yearstr, scalar_attr, scalar_bounds=[vmin, vmax], extend='both', outfile=join(outdir, 'yearly', '{}_k{}_{}.pdf'.format(variables_str, str(k), yearstr)), lats_lons=lats_lons) 
                
        elif seasonal: #Case where we plot one season per year
            
            season_months, season_years = get_season_months_and_years(season_start, season_end)
            seas_monthstr = season_months[0] + "-" + season_months[-1] #For titles
            
            data_seasons = []
            
            for i in range(years): #Iterate over specified years
                
                year = startyr + i
                yearstr = str(year)
                
                #Get seasonally-averaged data
                
                scalar_seas_file = join(seasonaldatdir, "avg_"+scalar_attr+"_"+season_start+yearstr+"-"+season_end+"*")
                ds_scalar_seas = xr.open_mfdataset(scalar_seas_file, engine="scipy")
                ds_scalar_seas.load()
                
                #Convert scalar DataSet to useful field
                lon_centers, lat_centers, lon_edges, lat_edges, scalar_seas = ds_to_field(ds_grid, ds_scalar_seas.isel(k=k), scalar_attr, latmin, latmax, lonmin, lonmax, resolution)
                
                seas_yearstr = yearstr + "-" + str(year + season_years[-1]) #For titles
                
                #Plot seasonal average
                ArcCir_pcolormesh(ds_grid, k, [scalar_seas], resolution, cmap, lon_centers, lat_centers, '{}, {}'.format(seas_monthstr, seas_yearstr), scalar_attr, scalar_bounds=[vmin, vmax], extend='both', outfile=join(outdir, 'seasonal', '{}_k{}_{}_{}.pdf'.format(variables_str, str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)
                
                data_seasons.append(scalar_seas)
                
            seas_yearstr = str(startyr) + "-" + str(startyr + (years-1) + season_years[-1]) #For titles
            
            #Plot average over all seasons
            ArcCir_pcolormesh(ds_grid, k, data_seasons, resolution, cmap, lon_centers, lat_centers, '{}, {}'.format(seas_monthstr, seas_yearstr), scalar_attr, scalar_bounds=[vmin, vmax], extend='both', outfile=join(outdir, 'seasonal', '{}_k{}_{}_{}.pdf'.format(variables_str, str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)
            
##############################

#CASE WHERE VECTOR AND SCALAR ARE PROVIDED

if include_vector_field:
    
    for k in range(kmin, kmax + 1): #Iterate over specified depths
        
        if not seasonal: #Case where we plot every month
            
            for i in range(years): #Iterate over specified years
                
                year = startyr + i
                yearstr = str(year)
                
                for m in range(12): #Iterate over months
                    
                    monthstr = get_monthstr(m)