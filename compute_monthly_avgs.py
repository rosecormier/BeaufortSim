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
import xarray as xr

from os.path import expanduser, join

from ecco_general import get_monthstr, load_dataset, load_grid, get_vector_in_xy
from ecco_field_variables import get_field_vars, get_variable_str 
from geostrophic_functions import get_density_and_pressure, comp_geos_vel
from vorticity_functions import comp_vorticity, comp_normal_strain, comp_shear_strain

#To be called within this script
import download_new_data

##############################

def main(**kwargs):
    
    if not kwargs:

        parser = argparse.ArgumentParser(description="Compute monthly geostrophic velocity, vorticity, etc. in Beaufort Gyre", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        #Spatial bounds

        parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, \
                            default=[70.0, 85.0])
        parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, \
                            default=[-180.0, -90.0])

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

        #Temporal bounds

        startyr = config['start']
        years = config['years']

        #Directories

        datdirshort = config['datdir']
        outdir = join(".", config['outdir'])
        
    elif kwargs:
        
        latmin, latmax = kwargs.get('latmin'), kwargs.get('latmax')
        lonmin, lonmax = kwargs.get('lonmin'), kwargs.get('lonmax')
        startyr = kwargs.get('startyr')
        years = kwargs.get('years')
        datdirshort = kwargs.get('datdir')
        outdir = kwargs.get('outdir')
        
    user_home_dir = expanduser('~')
    sys.path.append(join(user_home_dir, 'ECCOv4-py'))
    datdir = join(user_home_dir, datdirshort, 'ECCO_V4r4_PODAAC')

    ug_monthly_shortname, ug_monthly_nc_str = get_field_vars('UGVG')
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
            
            if not os.path.exists(curr_denspress_file): #If the file doesn't exist, download it
                download_new_data.main(startmo="01", startyr=year, months=12, scalars=["PHIHYDcR"], xvectors=None, datdir="Downloads")

            #Load monthly density-/pressure-anomaly file into workspace
            ds_denspress_mo = load_dataset(curr_denspress_file) 

            dens, press = get_density_and_pressure(ds_denspress_mo)

            u_g, v_g = comp_geos_vel(ds_grid, press, dens) #Compute geostrophic velocity components
            u_g.name = 'u_g'
            v_g.name = 'v_g'
            
            #Save geostrophic velocity components

            vel_g_ds = xr.merge([u_g, v_g])
            vel_g_ds.to_netcdf(path=join(outdir, ug_monthly_shortname, ug_monthly_nc_str+yearstr+"-"+monthstr+".nc"), engine="scipy")

            ##############################

            #VORTICITY

            curr_vel_file = join(vel_dir, vel_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
            
            if not os.path.exists(curr_vel_file): #If the file doesn't exist, download it
                download_new_data.main(startmo="01", startyr=year, months=12, scalars=None, xvectors=["UVEL"], datdir="Downloads")
            
            ds_vel_mo = load_dataset(curr_vel_file) #Load monthly velocity file into workspace

            (ds_vel_mo['UVEL']).data, (ds_vel_mo['VVEL']).data = (ds_vel_mo['UVEL']).values, (ds_vel_mo['VVEL']).values

            xgcm_grid = ecco.get_llc_grid(ds_grid)

            ZETA = comp_vorticity(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz)

            ZETA.name = 'ZETA'
            ZETA.to_netcdf(path=join(outdir, zeta_monthly_shortname, zeta_monthly_nc_str+yearstr+"-"+monthstr+".nc"), engine="scipy")

            ##############################

            #NORMAL STRAIN

            normal_strain = comp_normal_strain(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxG, ds_grid.dyG, ds_grid.rA)

            normal_strain.to_netcdf(path=join(outdir, normal_monthly_shortname, normal_monthly_nc_str+yearstr+"-"+monthstr+".nc"), engine="scipy")

            ##############################

            #SHEAR STRAIN

            shear_strain = comp_shear_strain(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz)

            shear_strain.to_netcdf(path=join(outdir, shear_monthly_shortname, shear_monthly_nc_str+yearstr+"-"+monthstr+".nc"), engine="scipy")
            
##############################

if __name__ == "__main__":
    main()