"""
Compute and plot vorticity and Okubo-Weiss parameter.

Rosalie Cormier, 2023
"""

#IMPORTS AND SETTINGS

import os
import sys
import argparse
import ecco_v4_py as ecco
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from os.path import expanduser, join

from ecco_general import load_grid, get_starting_i, get_monthstr, load_dataset, get_vector_in_xy, rotate_vector, ds_to_field, comp_temp_mean, comp_residuals, ecco_resample
from ecco_field_variables import get_vector_field_vars
from ecco_visualization import ArcCir_pcolormesh, ArcCir_contourf_quiver
from vorticity_functions import comp_vorticity, comp_total_strain, comp_OkuboWeiss

seis_nanmasked = plt.get_cmap("seismic").copy()
seis_nanmasked.set_bad('grey')

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot vorticity and Okubo-Weiss in Beaufort Gyre", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, \
                    default=[70.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, \
                    default=[-180.0, -90.0])
parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=24)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 1])
parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, default=0.25)
parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", default="visualization")

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

outdir = join(".", config['outdir'], 'vorticity')

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
#Get variables associated with velocity
vel_monthly_shortname, vel_monthly_nc_str, vel_variable = get_vector_field_vars('UVEL', 'VVEL')
    
##############################

#LOAD GRID AND GET LISTS OF MONTHS/YEARS

ds_grid = load_grid(datdir)

year = startyr
monthstrs, yearstrs = [], []

i = get_starting_i(startmo)

#Iterate over all specified months
while i <= mos:
    
    monthstr, yearstr = get_monthstr(i), str(year)
    monthstrs.append(monthstr)
    yearstrs.append(yearstr)

    if (i + 1) % 12 == 0 and (i + 1) != get_starting_i(startmo) + mos:
        year += 1 #Go to next year
        
    i += 1 #Go to next month

##############################

#GET FILE LIST
vel_dir = join(datdir, vel_monthly_shortname)

##############################

#CREATE PLOTS

#Iterate over all specified depths
for k in range(kmin, kmax + 1):
    
    ds_vels, zetas = [], []
    velEs, velNs = [], []
    
    for m in range(mos):
        
        monthstr, yearstr = monthstrs[m], yearstrs[m]
        
        curr_vel_file = join(vel_dir, vel_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        
        ds_vel_mo = load_dataset(curr_vel_file) #Load monthly velocity file into workspace
        ds_vels.append(ds_vel_mo)
        
        #Interpolate velocities to centres of grid cells
        (ds_vel_mo['UVEL']).data, (ds_vel_mo['VVEL']).data = (ds_vel_mo['UVEL']).values, (ds_vel_mo['VVEL']).values
        velE, velN = rotate_vector(ds_grid, ds_vel_mo, 'UVEL', 'VVEL')
        velE, velN = (velE.isel(k=k)).squeeze(), (velN.isel(k=k)).squeeze()
        
        velEs.append(velE)
        velNs.append(velN)

        xgcm_grid = ecco.get_llc_grid(ds_grid)
        
        #Compute vorticity and resample to lat-lon grid
        
        zeta = comp_vorticity(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz)
        lon_centers, lat_centers, lon_edges, lat_edges, zeta_field = ecco_resample(ds_grid, zeta.isel(k=k), latmin, latmax, lonmin, lonmax, resolution)

        zetas.append(zeta_field)
    
    #Compute and plot annual average vorticity
    zeta_mean = ArcCir_pcolormesh(ds_grid, k, zetas, resolution, 'PuOr', lon_centers, lat_centers, yearstr, scalar_attr='zeta', scalar_bounds=[-1e-7, 1e-7], outfile=join(outdir, 'zeta_k{}_avg{}.png'.format(str(k), yearstr)), lats_lons=[70.0, 85.0, -180.0, -90.0])
    
    #Compute residuals of monthly averages
    zeta_residuals = comp_residuals(zetas, zeta_mean) #Plot this?
    
    #Concatenate and average velocities
    
    ds_vels_concat = xr.concat(ds_vels, dim='time')
    ds_vels_mean = ds_vels_concat.sum(dim='time') / len(ds_vels)
    (ds_vels_mean['UVEL']).data, (ds_vels_mean['VVEL']).data = (ds_vels_mean['UVEL']).values, (ds_vels_mean['VVEL']).values

    #Compute annual average W
    W = comp_OkuboWeiss(xgcm_grid, ds_grid, ds_vels_mean['UVEL'], ds_vels_mean['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.dxG, ds_grid.dyG, ds_grid.rA, ds_grid.rAz, k, latmin, latmax, lonmin, lonmax, resolution)
    
    W[np.isnan(W)] = 0
    W = W * 1e12 #Useful for computing standard deviation
    W[W > 0] = 1.0
    W[W < 0] = -1.0
    W[W == 0] = np.nan
    
    #Set Okubo-Weiss threshold
    sigma_W = 0.2 * np.std(W) #Not used at the moment
    
    #Plot annual average W
    ArcCir_pcolormesh(ds_grid, k, [W], resolution, seis_nanmasked, lon_centers, lat_centers, yearstr, scalar_attr='W', scalar_bounds=[-1.5, 1.5], outfile=join(outdir, 'W_k{}_avg{}.png'.format(str(k), yearstr)), lats_lons=[70.0, 85.0, -180.0, -90.0])
    
    ArcCir_contourf_quiver(ds_grid, k, [W], velEs, velNs, resolution, seis_nanmasked, yearstr, lon_centers, lat_centers, lon_edges, lat_edges, scalar_attr='W', scalar_bounds=[-1.5, 1.5], outfile=join(outdir, 'W_k{}_avg{}_contour.png'.format(str(k), yearstr)), no_levels=3)