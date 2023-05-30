"""
Rosalie Cormier, 2023
"""

##############################

#IMPORTS AND SETTINGS

import sys
import os
import argparse
import ecco_v4_py as ecco
import glob
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

from os.path import expanduser, join

from ecco_general import load_grid, get_monthstr, load_dataset, ds_to_field, comp_residuals, get_starting_i, rotate_vector
from ecco_visualization import cbar_label, contourf_quiver_title, ArcCir_contourf_quiver, ArcCir_contourf_quiver_grid
from ecco_field_variables import get_scalar_field_vars, get_vector_field_vars

vir_nanmasked = plt.get_cmap('viridis_r').copy()
vir_nanmasked.set_bad('black')

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot scalar and vector fields in Beaufort Gyre",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, \
                    default=[70.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, \
                    default=[-180.0, -90.0])
parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 4])
parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, \
                    default=1.0)
parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", \
                    default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", \
                    default="visualization")

parser.add_argument("start", type=int, help="Start year")

args = parser.parse_args()
config = vars(args)

latmin, latmax = config['lats'][0], config['lats'][1]
lonmin, lonmax = config['lons'][0], config['lons'][1]
startmo, startyr, mos = config['month'], config['start'], config['months']
kmin, kmax = config['kvals'][0], config['kvals'][1]
resolution = config['res']

user_home_dir = expanduser('~')
sys.path.append(join(user_home_dir, 'ECCOv4-py'))
datdir = join(user_home_dir, config['datdir'], 'ECCO_V4r4_PODAAC')
  
outdir = join(".", config['outdir'])

if not os.path.exists(outdir):
    os.makedirs(outdir)

#Set parameters and get associated variables

scalar_attr = 'PHIHYDcR'
xvec_attr, yvec_attr = 'UVEL', 'VVEL'

vector_monthly_shortname, vector_monthly_nc_str, vector_variable = \
    get_vector_field_vars(xvec_attr, yvec_attr)
scalar_monthly_shortname, scalar_monthly_nc_str, scalar_variable = \
    get_scalar_field_vars(scalar_attr)
variables_str = vector_variable + '_' + scalar_variable
    
##############################

#LOAD GRID AND GET LISTS OF MONTHS/YEARS

ds_grid = load_grid(datdir)

year = startyr
monthstrs, yearstrs = [], []
    
i = get_starting_i(startmo)
    
#Iterate over all specified months
while i < get_starting_i(startmo) + mos:

    monthstr, yearstr = get_monthstr(i), str(year)
    monthstrs.append(monthstr)
    yearstrs.append(yearstr)

    if (i + 1) % 12 == 0 and (i + 1) != get_starting_i(startmo) + mos:
        year += 1 #Go to next year
        
    i += 1 #Go to next month

##############################

#GET FILE LISTS

vector_dir = join(datdir, vector_monthly_shortname)
scalar_dir = join(datdir, scalar_monthly_shortname)

##############################

#CREATE MONTHLY PLOTS OF VELOCITY AND PRESSURE ANOMALY

#Iterate over all specified depths
for k in range(kmin, kmax + 1):

    scalars, vecEs, vecNs = [], [], []

    for m in range(mos):

        monthstr, yearstr = monthstrs[m], yearstrs[m]

        curr_vector_file = join(vector_dir, vector_monthly_nc_str+yearstr+ \
                                "-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        curr_scalar_file = join(scalar_dir, scalar_monthly_nc_str+yearstr+ \
                                "-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")

        ds_vector_mo = load_dataset(curr_vector_file) #Load monthly vector file into workspace

        #Interpolate vectors to centres of grid cells
        (ds_vector_mo[xvec_attr]).data, (ds_vector_mo[yvec_attr]).data = \
            (ds_vector_mo[xvec_attr]).values, (ds_vector_mo[yvec_attr]).values
        vecE, vecN = rotate_vector(ds_grid, ds_vector_mo, xvec_attr, yvec_attr)
        vecE, vecN = (vecE.isel(k=k)).squeeze(), (vecN.isel(k=k)).squeeze()
        
        ds_scalar_mo = load_dataset(curr_scalar_file) #Load monthly scalar file into workspace
        
        #Convert scalar DataSet to useful field
        lon_centers, lat_centers, lon_edges, lat_edges, scalar = ds_to_field(ds_grid, ds_scalar_mo.isel(k=k), scalar_attr, k, latmin, latmax, lonmin, lonmax, resolution)
        
        scalars.append(scalar)
        vecEs.append(vecE) 
        vecNs.append(vecN)
        
        #Plot scalar with vector
        ArcCir_contourf_quiver(ds_grid, k, [scalar], [vecE], [vecN], resolution, vir_nanmasked, [93, 97], yearstr+"-"+monthstr, lon_centers, lat_centers, lon_edges, lat_edges, outfile=join(outdir, '{}_k{}_{}-{}.pdf'.format(variables_str, \
                                                                      str(k), \
                                                                      monthstr, \
                                                                      yearstr)))

    #Plot all months
    ArcCir_contourf_quiver_grid(ds_grid, k, scalars, vecEs, vecNs,
                                [93, 97], resolution, vir_nanmasked,  \
                                monthstrs, yearstrs, lon_centers, lat_centers, lon_edges, lat_edges, \
                                outfile=join(outdir, '{}_k{}_all{}.png'.format(variables_str, \
                                                                               str(k), \
                                                                               yearstr)))

    #Plot annual averages
    scalar_mean, vecE_mean, vecN_mean = ArcCir_contourf_quiver(ds_grid, k, scalars, vecEs, vecNs, \
                                                      resolution, vir_nanmasked, [93, 97], \
                                                      yearstrs[0]+" average", \
                                                      lon_centers, lat_centers, lon_edges, lat_edges, outfile=join(outdir, \
                                                                   '{}_k{}_avg{}.pdf'.format(variables_str, \
                                                                                             str(k), \
                                                                                             yearstr)))
    
    #Compute annual average magnitude of vector field
    magnitude_mean = np.sqrt(vecE_mean**2 + vecN_mean**2)
    print("Vector maximum magnitude =", np.nanmax(magnitude_mean), "; vector minimum magnitude =", np.nanmin(magnitude_mean))
    
    #Compute residuals of monthly averages
    scalar_residuals = comp_residuals(scalars, scalar_mean)
    vecE_residuals, vecN_residuals = comp_residuals(vecEs, vecE_mean), comp_residuals(vecNs, vecN_mean)

    #Plot residuals for all months
    ArcCir_contourf_quiver_grid(ds_grid, k, scalar_residuals, vecE_residuals, vecN_residuals,
                                [-2, 2], resolution, vir_nanmasked,  \
                                monthstrs, yearstrs, lon_centers, lat_centers, lon_edges, lat_edges, \
                                outfile=join(outdir, '{}_k{}_all{}.png'.format(variables_str, \
                                                                               str(k), \
                                                                               yearstr)), resid=True)