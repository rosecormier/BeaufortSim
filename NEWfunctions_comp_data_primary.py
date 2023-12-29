"""
Contains functions that directly produce computed data.

To be added
    -Geostrophic velocity
    -Vorticity
    -Normal strain
    -Shear strain
    -Horizontal vel div
    -Ekman velocity

Only does monthly averaging at the moment; will update.

Rosalie Cormier, 2023
"""

import xarray as xr

from os.path import join

import download_data

from functions_ecco_general import load_dataset, load_grid
from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string

from functions_geostrophy import * #Ideally, I'd like to move all these function definitions to be defined within the respective functions here
    #But for now, just import them

##############################

def comp_geostrophic_vel(monthstr, yearstr, datdir_primary, datdir_secondary, time_ave_type='monthly'):
    
    """
    Computes and saves geostrophic velocity components (u_g, v_g) to DataSet in NetCDF format.
    """
    
    if time_ave_type == 'monthly':
        
        date_string = yearstr + '-' + monthstr
        
        #Look for the density/pressure file in primary directory and download if it doesnt exist
        download_data.main(field_name='density', initial_month=monthstr, initial_year=yearstr, final_month=monthstr, final_year=yearstr, time_ave_type=time_ave_type, datdir_primary=datdir_primary)
        
        denspress_shortname, denspress_nc_string = get_monthly_shortname(get_field_variable('density')), get_monthly_nc_string(get_field_variable('density'))
        denspress_file = join(datdir_primary, denspress_shortname, denspress_nc_string+date_string+"_ECCO_V4r4_native_llc0090.nc")
        ds_denspress = load_dataset(denspress_file) #Load the density-/pressure-anomaly DataSet into workspace
    
        #Call this function to extract the density and pressure from the DataSet
        density, pressure = get_density_and_pressure(ds_denspress, rho_ref)
    
        u_g, v_g = comp_geos_vel(ds_grid, pressure, density) #Compute geostrophic velocity components
        
        #Save the data
        
        u_g.name = 'UG'
        v_g.name = 'VG'
        vel_g_ds = xr.merge([u_g, v_g])
        vel_g_shortname, vel_g_nc_string = get_monthly_shortname(get_field_variable('geostrophic_vel')), get_monthly_nc_string(get_field_variable('geostrophic_vel'))
        vel_g_ds.to_netcdf(path=join(datdir_secondary, vel_g_shortname, vel_g_nc_string+date_string+".nc"), engine="scipy")

##############################                

#each field gets a separate function
    #some of them will call other functions, as appropriate

#horizontal velocity div

#normal strain

#shear strain

#Ekman velocity
    #First load ocean surface stresses (primary data)
    
##############################

def main(**kwargs):
    
    if kwargs:
        
        #Directories
        
        datdir_primary = kwargs.get('datdir_primary')
        datdir_secondary = kwargs.get('datdir_secondary')
    
    ds_grid = load_grid(datdir_primary) #Load the grid DataSet

    #then identify the input variable and run the appropriate function

##############################

if __name__ == "__main__":
    main()