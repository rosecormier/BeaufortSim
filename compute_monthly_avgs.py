"""
Compute and save monthly average geostrophic velocity, vorticity, and associated values.

Rosalie Cormier, 2023
"""

##############################

#IMPORTS AND SETTINGS

import os
import sys
import argparse
import ecco_v4_py as ecco

from os.path import expanduser, join

from ecco_general import get_monthstr, load_dataset, load_grid, get_vector_in_xy
from ecco_field_variables_new import get_field_vars, get_variable_str #
from geostrophic_functions import get_density_and_pressure, comp_geos_vel
from vorticity_functions import comp_vorticity, comp_normal_strain, comp_shear_strain

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Compute monthly geostrophic velocity, vorticity, etc. in Beaufort Gyre", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#Spatial bounds

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, \
                    default=[70.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, \
                    default=[-180.0, -90.0])
#parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, default=0.25)

#Temporal bounds

parser.add_argument("start", type=int, help="Start year") #Required
parser.add_argument("--years", type=int, help="Total number of years to iterate over", default=1)
    
#Directories

parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", default="computed_monthly")

args = parser.parse_args()
config = vars(args)

#Spatial bounds

latmin, latmax = config['lats'][0], config['lats'][1]
lonmin, lonmax = config['lons'][0], config['lons'][1]
#resolution = config['res']

#Temporal bounds

startyr = config['start']
years = config['years']

#Directories

user_home_dir = expanduser('~')
sys.path.append(join(user_home_dir, 'ECCOv4-py'))
datdir = join(user_home_dir, config['datdir'], 'ECCO_V4r4_PODAAC')

outdir = join(".", config['outdir'])

ug_monthly_shortname, ug_monthly_nc_str = get_field_vars('UG')
vg_monthly_shortname, vg_monthly_nc_str = get_field_vars('VG')
zeta_monthly_shortname, zeta_monthly_nc_str = get_field_vars('ZETA')
normal_monthly_shortname, normal_monthly_nc_str = get_field_vars('NORMAL')
shear_monthly_shortname, shear_monthly_nc_str = get_field_vars('SHEAR')

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
for subdir in [ug_monthly_shortname, zeta_monthly_shortname, \
              normal_monthly_shortname]:
    if not os.path.exists(join(outdir, subdir)):
        os.makedirs(join(outdir, subdir))
    
#Get variables associated with fields

denspress_monthly_shortname, denspress_monthly_nc_str = get_field_vars('PHIHYDcR')
vel_monthly_shortname, vel_monthly_nc_str = get_field_vars('UVELVVEL')

rho_ref = 1029.0 #Reference density (kg/m^3)
f_mean = 1e-4 #"Typical" Coriolis parameter (1/s)

##############################

#GET FILE LISTS

denspress_dir = join(datdir, denspress_monthly_shortname)
vel_dir = join(datdir, vel_monthly_shortname)

##############################

#LOAD GRID AND ITERATE TO COMPUTE MONTHLY VALUES

ds_grid = load_grid(datdir)

for i in range(years):

    year = startyr + i
    yearstr = str(year)

    for m in range(12):

        monthstr = get_monthstr(m)
            
        ##############################
            
        #GEOSTROPHIC VELOCITY
            
        curr_denspress_file = join(denspress_dir, denspress_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
            
        #Load monthly density-/pressure-anomaly file into workspace
        ds_denspress_mo = load_dataset(curr_denspress_file) 
        
        dens, press = get_density_and_pressure(ds_denspress_mo)
        
        u_g, v_g = comp_geos_vel(ds_grid, press, dens) #Compute geostrophic velocity components
        
        #Save geostrophic velocity components
        
        u_g.to_netcdf(path=join(outdir, ug_monthly_shortname, ug_monthly_nc_str+yearstr+"-"+monthstr+".nc"), engine="scipy")
        v_g.to_netcdf(path=join(outdir, vg_monthly_shortname, vg_monthly_nc_str+yearstr+"-"+monthstr+".nc"), engine="scipy")
        
        ##############################
        
        #VORTICITY
        
        curr_vel_file = join(vel_dir, vel_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        
        ds_vel_mo = load_dataset(curr_vel_file) #Load monthly velocity file into workspace
        
        (ds_vel_mo['UVEL']).data, (ds_vel_mo['VVEL']).data = (ds_vel_mo['UVEL']).values, (ds_vel_mo['VVEL']).values
        
        xgcm_grid = ecco.get_llc_grid(ds_grid)
        
        zeta = comp_vorticity(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz)
        
        zeta.to_netcdf(path=join(outdir, zeta_monthly_shortname, zeta_monthly_nc_str+yearstr+"-"+monthstr+".nc"), engine="scipy")
        
        ##############################
        
        #NORMAL STRAIN
        
        normal_strain = comp_normal_strain(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxG, ds_grid.dyG, ds_grid.rA)
        
        normal_strain.to_netcdf(path=join(outdir, normal_monthly_shortname, normal_monthly_nc_str+yearstr+"-"+monthstr+".nc"), engine="scipy")
        
        ##############################
        
        #SHEAR STRAIN
        
        shear_strain = comp_shear_strain(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz)
        
        shear_strain.to_netcdf(path=join(outdir, shear_monthly_shortname, shear_monthly_nc_str+yearstr+"-"+monthstr+".nc"), engine="scipy")