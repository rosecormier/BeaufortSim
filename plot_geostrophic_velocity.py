"""
Compute and plot yearly average geostrophic velocity and Delta-u.

Rosalie Cormier, 2023, based on code by Andrew Delman
"""

##############################

#IMPORTS AND SETTINGS
 
import xgcm
import os
import sys
import argparse
import xarray as xr
import matplotlib.pyplot as plt
import ecco_v4_py as ecco
import numpy as np 

from os.path import expanduser, join

from ecco_general import load_grid, get_monthstr, get_starting_i, load_dataset, ds_to_field, get_vector_in_xy, comp_temp_mean, ecco_resample
from ecco_field_variables import get_scalar_field_vars, get_vector_field_vars
from geostrophic_functions import get_density_and_pressure, comp_geos_vel, rotate_u_g, comp_delta_u_norm, mask_delta_u
from ecco_visualization import ArcCir_contourf_quiver, ArcCir_pcolormesh

vir_nanmasked = plt.get_cmap('viridis_r').copy()
vir_nanmasked.set_bad('black')

red_nanmasked = plt.get_cmap("Reds").copy()
red_nanmasked.set_bad('grey')

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot geostrophic velocity in Beaufort Gyre", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, \
                    default=[70.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, \
                    default=[-180.0, -90.0])
parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[12, 13])
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
    
outdir = join(".", config['outdir'], 'geostrophic')

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
#Get variables associated with velocity and density/pressure

vel_monthly_shortname, vel_monthly_nc_str, vel_variable = get_vector_field_vars('UVEL', 'VVEL', geostrophic=True)
denspress_monthly_shortname, denspress_monthly_nc_str, denspress_variable = get_scalar_field_vars('PHIHYDcR')
variables_str = vel_variable + '_' + denspress_variable

rho_ref = 1029.0 #Reference density (kg/m^3)

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

#GET FILE LISTS

vel_dir = join(datdir, vel_monthly_shortname)
denspress_dir = join(datdir, denspress_monthly_shortname)

##############################

#CREATE MONTHLY PLOTS OF GEOSTROPHIC VELOCITY

#Iterate over all specified depths
for k in range(kmin, kmax + 1):
    
    u_list, u_g_list, u_a_list = [], [], []
    
    for m in range(mos):
        
        monthstr, yearstr = monthstrs[m], yearstrs[m]
        
        curr_vel_file = join(vel_dir, vel_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        curr_denspress_file = join(denspress_dir, denspress_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        
        ds_vel_mo = load_dataset(curr_vel_file) #Load monthly velocity file into workspace
        
        #Interpolate velocities to centres of grid cells
        
        (ds_vel_mo['UVEL']).data, (ds_vel_mo['VVEL']).data = (ds_vel_mo['UVEL']).values, (ds_vel_mo['VVEL']).values
        velocity_interp = get_vector_in_xy(ds_grid, ds_vel_mo, 'UVEL', 'VVEL') 
        u, v = velocity_interp['X'], velocity_interp['Y']
        u, v = (u.isel(k=k)).squeeze(), (v.isel(k=k)).squeeze()

        #Load monthly density-/pressure-anomaly file into workspace
        ds_denspress_mo = load_dataset(curr_denspress_file) 
        
        dens, press = get_density_and_pressure(ds_denspress_mo)
        
        u_g, v_g = comp_geos_vel(ds_grid, press, dens) #Compute geostrophic velocity components
        u_g, v_g = rotate_u_g(ds_grid, u_g, v_g, k)
        
        #Convert pressure-anomaly DataSet to useful field
        lon_centers, lat_centers, lon_edges, lat_edges, pressure = ds_to_field(ds_grid, ds_denspress_mo.isel(k=k), 'PHIHYDcR', latmin, latmax, lonmin, lonmax, resolution)
        
        #Compute ageostrophic velocity 
        u_a, v_a = u - u_g, v - v_g

        u_list.append(u + 1j * v)
        u_g_list.append(u_g + 1j * v_g)
        u_a_list.append(u_a + 1j * v_a)
        
    #Temporally average geostrophic and ageostrophic velocities, and regular velocity
    
    u_g_mean = comp_temp_mean(u_g_list)
    u_a_mean = comp_temp_mean(u_a_list)
    u_mean = comp_temp_mean(u_list)
    
    #Plot geostrophic velocity with pressure
    ArcCir_contourf_quiver(ds_grid, k, [pressure], [np.real(u_g_mean)], [np.imag(u_g_mean)], resolution, vir_nanmasked, yearstrs[0]+" average (geostrophic velocity)", lon_centers, lat_centers, outfile=join(outdir, '{}_k{}_all{}.pdf'.format(variables_str, \
                                                                      str(k), \
                                                                      yearstr)))

    #Compute Delta-u metric
    Delta_u = comp_delta_u_norm(ds_grid, k, u_mean, u_g_mean)
    
    #Plot Delta-u
    
    ds_grid_copy = ds_grid.copy()
    lon_centers, lat_centers, lon_edges, lat_edges, Delta_u_plot = ecco_resample(ds_grid_copy, Delta_u, latmin, latmax, lonmin, lonmax, resolution)
    
    ArcCir_pcolormesh(ds_grid, k, [Delta_u_plot], resolution, 'Reds', lon_centers, lat_centers, yearstr, scalar_attr="Delta_u", outfile=join(outdir, 'Delta_u_k{}_all{}.pdf'.format(str(k), yearstr)))
    
    #Repeat with small velocities masked
    Delta_u = comp_delta_u_norm(ds_grid, k, u_mean, u_g_mean, mask=mask_delta_u(0.005, u_mean))
    
    ds_grid_copy = ds_grid.copy()
    lon_centers, lat_centers, lon_edges, lat_edges, Delta_u_plot = ecco_resample(ds_grid_copy, Delta_u, latmin, latmax, lonmin, lonmax, resolution)
    
    ArcCir_pcolormesh(ds_grid, k, [Delta_u_plot], resolution, red_nanmasked, lon_centers, lat_centers, yearstr, scalar_attr="Delta_u", outfile=join(outdir, 'Delta_u_mask_k{}_all{}.pdf'.format(str(k), yearstr)))