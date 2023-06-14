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

from ecco_general import load_grid, get_starting_i, get_monthstr, load_dataset, rotate_vector, comp_temp_mean, comp_residuals, ecco_resample
from ecco_field_variables import get_scalar_field_vars, get_vector_field_vars
from ecco_visualization import ArcCir_pcolormesh, ArcCir_contourf_quiver
from geostrophic_functions import get_density_and_pressure, comp_geos_vel
from vorticity_functions import comp_vorticity, comp_normal_strain, comp_shear_strain, comp_total_strain, comp_OkuboWeiss, comp_local_Ro

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
parser.add_argument("--iters", type=int, help="Total number of months/seasons to iterate over", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 1])
parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, default=0.25)
parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", default="visualization")
parser.add_argument("--seasonal", type=bool, help="Whether to take seasonal averages rather than continuous averages", \
                    default=False)
parser.add_argument("--season", type=str, help="Months marking start and end of 'season'", default=["01", "01"])

parser.add_argument("start", type=int, help="Start year")

args = parser.parse_args()
config = vars(args)

latmin, latmax = config['lats'][0], config['lats'][1]
lonmin, lonmax = config['lons'][0], config['lons'][1]
resolution = config['res']

kmin, kmax = config['kvals'][0], config['kvals'][1]

startmo, startyr = config['month'], config['start']
seasonal = config['seasonal']

if not seasonal:
    mos = config['iters']

elif seasonal:
    seasons = config['iters']
    season_start, season_end = config['season'][0], config['season'][1]

user_home_dir = expanduser('~')
sys.path.append(join(user_home_dir, 'ECCOv4-py'))
datdir = join(user_home_dir, config['datdir'], 'ECCO_V4r4_PODAAC')

outdir = join(".", config['outdir'], 'vorticity')

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
for subfolder in ["yearly"]:
    if not os.path.exists(join(outdir, subfolder)):
        os.makedirs(join(outdir, subfolder))
    
#Get variables associated with velocity
vel_monthly_shortname, vel_monthly_nc_str, vel_variable = get_vector_field_vars('UVEL', 'VVEL')

#Get variables associated with density/pressure
denspress_monthly_shortname, denspress_monthly_nc_str, denspress_variable = get_scalar_field_vars('PHIHYDcR')

f_mean = 1e-4 #Set "typical" Coriolis parameter (1/s)
    
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
    
datestr_start, datestr_end = monthstrs[0] + yearstrs[0], monthstrs[-1] + yearstrs[-1]

##############################

#GET FILE LISTS

vel_dir = join(datdir, vel_monthly_shortname)
denspress_dir = join(datdir, denspress_monthly_shortname)

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

        zetas.append(zeta_field/f_mean) #Normalize by f and save

    #Compute and plot annual average vorticity, normalized by f
    zeta_mean = ArcCir_pcolormesh(ds_grid, k, zetas, resolution, 'PuOr', lon_centers, lat_centers, yearstr, scalar_attr='zeta', scalar_bounds=[-1e-3, 1e-3], outfile=join(outdir, 'yearly', 'zeta_k{}_avg{}-{}.png'.format(str(k), datestr_start, datestr_end)), lats_lons=[70.0, 85.0, -180.0, -90.0])
    
    #Compute residuals of monthly averages
    zeta_residuals = comp_residuals(zetas, zeta_mean) #Plot this?
    
    zeta_mean *= f_mean #Recover actual vorticity
    
    #Compute and plot local Rossby number
    
    Ro_l = comp_local_Ro(zeta_mean, lat_centers)
    ArcCir_pcolormesh(ds_grid, k, [Ro_l], resolution, 'seismic', lon_centers, lat_centers, yearstr, scalar_attr='Ro_l', scalar_bounds=[-3e-3, 3e-3], outfile=join(outdir, 'yearly', 'localRo_k{}_avg{}.png'.format(str(k), yearstr)), lats_lons=[70.0, 85.0, -180.0, -90.0])
    
    #Concatenate and average velocities
    
    ds_vels_concat = xr.concat(ds_vels, dim='time')
    ds_vels_mean = ds_vels_concat.sum(dim='time') / len(ds_vels)
    (ds_vels_mean['UVEL']).data, (ds_vels_mean['VVEL']).data = (ds_vels_mean['UVEL']).values, (ds_vels_mean['VVEL']).values

    #Compute average strain components
    
    normal_strain = comp_normal_strain(xgcm_grid, ds_vels_mean['UVEL'], ds_vels_mean['VVEL'], ds_grid.dxG, ds_grid.dyG, ds_grid.rA)
    shear_strain = comp_shear_strain(xgcm_grid, ds_vels_mean['UVEL'], ds_vels_mean['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz)
    
    #Slice strain components in k- and time-dimensions
    normal_strain, shear_strain = normal_strain.isel(k=k).squeeze(), shear_strain.isel(k=k).squeeze()

    normal_strain= ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    shear_strain = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    
    #Compute annual average W
    W = comp_OkuboWeiss(zeta_mean, normal_strain, shear_strain)
    
    W[np.isnan(W)] = 0

    #W[W > 0] = 1.0 #Maps all positive W to single value
    #W[W < 0] = -1.0 #Maps all negative W to single value
    
    sigma_W = 0.01 * np.std(W) #Set Okubo-Weiss threshold
    W[abs(W) < sigma_W] = 0 #Mask small W

    W[W == 0] = np.nan
    
    #Plot annual average W
    ArcCir_pcolormesh(ds_grid, k, [W], resolution, seis_nanmasked, lon_centers, lat_centers, yearstr, scalar_attr='W', scalar_bounds=[-1.5e-14, 1.5e-14], outfile=join(outdir, 'yearly', 'W_k{}_avg{}-{}.png'.format(str(k), datestr_start, datestr_end)), lats_lons=[70.0, 85.0, -180.0, -90.0])
    
    ArcCir_contourf_quiver(ds_grid, k, [W], velEs, velNs, resolution, seis_nanmasked, yearstr, lon_centers, lat_centers, scalar_attr='W', scalar_bounds=[-1.5e-14, 1.5e-14], outfile=join(outdir, 'yearly', 'W_k{}_avg{}-{}_contour.png'.format(str(k), datestr_start, datestr_end)), no_levels=10)
    
##############################

#REPEAT FOR GEOSTROPHIC VELOCITIES

#Iterate over all specified depths
for k in range(kmin, kmax + 1):

    densities, pressures = [], []

    for m in range(mos):
        
        monthstr, yearstr = monthstrs[m], yearstrs[m]
        
        curr_denspress_file = join(denspress_dir, denspress_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        
        #Load monthly density-/pressure-anomaly file into workspace
        ds_denspress_mo = load_dataset(curr_denspress_file) 
        
        dens, press = get_density_and_pressure(ds_denspress_mo)
        densities.append(dens.isel(k=k).squeeze())
        pressures.append(press.isel(k=k).squeeze())
    
    vel_file = join(vel_dir, vel_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
    ds_vel = load_dataset(vel_file) #Load monthly velocity file into workspace
 
    #Compute average pressure and density fields
    
    press_mean = comp_temp_mean(pressures)
    dens_mean = comp_temp_mean(densities)
    
    #Compute annual average geostrophic velocity components
    u_g_mean, v_g_mean = comp_geos_vel(ds_grid, press_mean, dens_mean)
    
    #Save to velocity dataset
    ds_vel = ds_vel.squeeze()
    ds_vel = ds_vel.transpose(..., 'k') #Reshape
    (ds_vel['UVEL']).data, (ds_vel['VVEL']).data = u_g_mean.values, v_g_mean.values
    
    u_g_E, u_g_N = rotate_vector(ds_grid, ds_vel, 'UVEL', 'VVEL')
    u_g_E, u_g_N = (u_g_E.isel(k=k)).squeeze(), (u_g_N.isel(k=k)).squeeze()
    
    #Compute and plot annual average vorticity from geostrophic data, normalized by f
    
    zeta_mean = comp_vorticity(xgcm_grid, ds_vel['UVEL'], ds_vel['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz)
    lon_centers, lat_centers, lon_edges, lat_edges, zeta_field = ecco_resample(ds_grid, zeta_mean.isel(k=k), latmin, latmax, lonmin, lonmax, resolution)
    ArcCir_pcolormesh(ds_grid, k, [zeta_field/f_mean], resolution, 'PuOr', lon_centers, lat_centers, yearstr, scalar_attr='zeta_geos', scalar_bounds=[-1e-3, 1e-3], outfile=join(outdir, 'yearly', 'zeta_k{}_geos_avg{}-{}.png'.format(str(k), datestr_start, datestr_end)), lats_lons=[70.0, 85.0, -180.0, -90.0])
    
    #Compute average strain components
    
    normal_strain = comp_normal_strain(xgcm_grid, ds_vel['UVEL'], ds_vel['VVEL'], ds_grid.dxG, ds_grid.dyG, ds_grid.rA)
    shear_strain = comp_shear_strain(xgcm_grid, ds_vel['UVEL'], ds_vel['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz)
    
    #Slice strain components in k- and time-dimensions
    normal_strain, shear_strain = normal_strain.isel(k=k).squeeze(), shear_strain.isel(k=k).squeeze()

    normal_strain= ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    shear_strain = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    
    #Compute annual average W
    W = comp_OkuboWeiss(zeta_field, normal_strain, shear_strain)
    
    W[np.isnan(W)] = 0

    #W[W > 0] = 1.0 #Maps all positive W to single value
    #W[W < 0] = -1.0 #Maps all negative W to single value
    
    sigma_W = 0.01 * np.std(W) #Set Okubo-Weiss threshold
    W[abs(W) < sigma_W] = 0 #Mask small W
    
    W[W == 0] = np.nan
    
    #Plot annual average W
    ArcCir_pcolormesh(ds_grid, k, [W], resolution, seis_nanmasked, lon_centers, lat_centers, yearstr, scalar_attr='W_geos', scalar_bounds=[-1.5e-14, 1.5e-14], outfile=join(outdir, 'yearly', 'W_k{}_geos_avg{}-{}.png'.format(str(k), datestr_start, datestr_end)), lats_lons=[70.0, 85.0, -180.0, -90.0])
    
    ArcCir_contourf_quiver(ds_grid, k, [W], [u_g_E], [u_g_N], resolution, seis_nanmasked, yearstr, lon_centers, lat_centers, scalar_attr='W_geos', scalar_bounds=[-1.5e-14, 1.5e-14], outfile=join(outdir, 'yearly', 'W_k{}_geos_avg{}-{}_contour.png'.format(str(k), datestr_start, datestr_end)), no_levels=10)