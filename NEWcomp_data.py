"""
Contains functions that directly produce computed data.
    -Geostrophic velocity
    -Ekman velocity

To be added
    -Vorticity
    -Normal strain
    -Shear strain
    -Horizontal vel div
    -Ekman velocity div
    
The main function picks out the appropriate function to compute the requested field.

Only does monthly averaging at the moment; will update.

Rosalie Cormier, 2023
"""

import os
import xarray as xr

from os.path import join
from math import pi

import download_data

from functions_ecco_general import load_dataset, load_grid
from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string

from functions_geostrophy import * #Ideally, I'd like to move all these function definitions to be defined within the respective functions here
    #But for now, just import them

##############################

#AUXILIARY VALUES/FUNCTIONS; USED DURING SUBSEQUENT COMPUTATIONS

Omega = (2 * pi) / 86164 #Earth angular velocity

def to_radians(angle):
    
    """
    Converts degrees to radians.
    """
    
    return angle * pi / 180

def comp_f(y):
    
    """
    Computes local Coriolis parameter.
    """
    
    lat = to_radians(y)
    return 2 * Omega * np.sin(lat)

def get_density_and_pressure(ds_denspress, rho_ref):
    
    """
    Gets density and pressure-like attributes from a DataSet, given a reference density.
    """
    
    density_anom = ds_denspress['RHOAnoma']
    density = density_anom + rho_ref #Compute absolute density
    
    #Update density DataArray
    
    density.name = 'RHO'
    density.attrs.update({'long_name': 'In-situ seawater density', 'units': 'kg m-3'})
    
    pressure_anom = ds_denspress.PHIHYDcR.copy() #Copy pressure data
    pressure_like = rho_ref * pressure_anom #Compute quantity to differentiate (units of density * pressure)
    
    return density, pressure

def get_vel_g_components(ds_grid, pressure_like, density):
    
    """
    Differentiates pressure-like quantity to compute geostrophic-velocity components.
    """
    
    xgcm_grid = ecco.get_llc_grid(ds_grid)
    
    #Compute derivatives of pressure in x and y
    
    d_press_dx = (xgcm_grid.diff(pressure_like, axis="X", boundary='extend')) / ds_grid.dxC
    d_press_dy = (xgcm_grid.diff(pressure_like, axis="Y", boundary='extend')) / ds_grid.dyC
    
    #Convert DataArray content from dask to np arrays
    d_press_dx.data, d_press_dy.data = d_press_dx.values, d_press_dy.values
    
    #Interpolate to centres of grid cells
    
    press_grads_interp = xgcm_grid.interp_2d_vector({'X': d_press_dx, 'Y': d_press_dy}, boundary='extend')
    dp_dx, dp_dy = press_grads_interp['X'], press_grads_interp['Y']
    dp_dx.name = 'dp_dx'
    dp_dy.name = 'dp_dy'
    
    GB_RHS_1, GB_RHS_2 = dp_dx / density, - dp_dy / density #Compute RHS of geostrophic-balance equations
    GB_RHS_1, GB_RHS_2 = GB_RHS_1.where(ds_grid.maskC), GB_RHS_2.where(ds_grid.maskC) #Mask land
    
    f = comp_f(ds_grid.YC) #Compute Coriolis param. from latitudes of grid cell centres
    u_g, v_g = GB_RHS_2 / f, GB_RHS_1 / f #Compute geostrophic velocity components
    u_g, v_g = u_g.where(ds_grid.maskC), v_g.where(ds_grid.maskC) #Mask land on result
    
    return u_g, v_g

##############################

#GEOSTROPHIC VELOCITY

def comp_geostrophic_vel(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, rho_ref, time_ave_type='monthly'):
    
    """
    Computes and saves geostrophic velocity components (u_g, v_g) to DataSet in NetCDF format.
    """
    
    if time_ave_type == 'monthly':
        
        #Look for the density/pressure file in primary directory and download if it doesn't exist
        download_data.main(field_name='density', initial_month=monthstr, initial_year=yearstr, final_month=monthstr, final_year=yearstr, time_ave_type=time_ave_type, datdir_primary=datdir_primary)
        
        date_string = yearstr + '-' + monthstr
        
        denspress_shortname, denspress_nc_string = get_monthly_shortname(get_field_variable('density')), get_monthly_nc_string(get_field_variable('density'))
        denspress_file = join(datdir_primary, denspress_shortname, denspress_nc_string+date_string+"_ECCO_V4r4_native_llc0090.nc")
        ds_denspress = load_dataset(denspress_file) #Load the density-/pressure-anomaly DataSet into workspace
    
        #Call this function to extract the density and pressure from the DataSet
        density, pressure_like = get_density_and_pressure(ds_denspress, rho_ref)
        
        #Compute geostrophic velocity components
        u_g, v_g = get_vel_g_components(ds_grid, pressure_like, density)
        
        #Save the data
        
        u_g.name = 'UG'
        v_g.name = 'VG'
        vel_g_ds = xr.merge([u_g, v_g])
        vel_g_shortname, vel_g_nc_string = get_monthly_shortname(get_field_variable('geostrophic_vel')), get_monthly_nc_string(get_field_variable('geostrophic_vel'))
        if not os.path.exists(join(datdir_secondary, vel_g_shortname)):
            os.makedirs(join(datdir_secondary, vel_g_shortname))
        vel_g_ds.to_netcdf(path=join(datdir_secondary, vel_g_shortname, vel_g_nc_string+date_string+".nc"), engine="scipy")

############################## 

#EKMAN VELOCITY

def comp_Ek_vel(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, rho_ref, nu_E, time_ave_type='monthly'):
    
    """
    Computes and saves Ekman velocity components (u_Ek, v_Ek) to DataSet in NetCDF format.
    """

    if time_ave_type == 'monthly':
        
        #Look for the density/pressure file in primary directory and download if it doesn't exist
        download_data.main(field_name='density', initial_month=monthstr, initial_year=yearstr, final_month=monthstr, final_year=yearstr, time_ave_type=time_ave_type, datdir_primary=datdir_primary)
        
        date_string = yearstr + '-' + monthstr
        
        denspress_shortname, denspress_nc_string = get_monthly_shortname(get_field_variable('density')), get_monthly_nc_string(get_field_variable('density'))
        denspress_file = join(datdir_primary, denspress_shortname, denspress_nc_string+date_string+"_ECCO_V4r4_native_llc0090.nc")
        ds_denspress = load_dataset(denspress_file) #Load the density-/pressure-anomaly DataSet into workspace
    
        #Call this function to extract the density and pressure from the DataSet
        #density, pressure_like = get_density_and_pressure(ds_denspress, rho_ref)
        
        #Look for the surface wind-on-ocean-stress file in primary directory and download if it doesn't exist
        download_data.main(field_name='wind_stress', initial_month=monthstr, initial_year=yearstr, final_month=monthstr, final_year=yearstr, time_ave_type=time_ave_type, datdir_primary=datdir_primary)
        
        stress_shortname, stress_nc_string = get_monthly_shortname(get_field_variable('wind_stress')), get_monthly_nc_string(get_field_variable('wind_stress'))
        stress_file = join(datdir_primary, stress_shortname, stress_nc_string+date_string+"_ECCO_V4r4_native_llc0090.nc")
        ds_stress = load_dataset(stress_file) #Load the stress DataSet into workspace
    
        #Compute Ekman velocity components
        u_Ek, v_Ek = comp_Ekman_vel(ds_grid, ds_denspress, ds_stress, nu_E, rho_ref)

        #Save the data
        
        u_Ek.name = 'UEk'
        v_Ek.name = 'VEk'
        vel_Ek_ds = xr.merge([u_Ek, v_Ek])
        vel_Ek_shortname, vel_Ek_nc_string = get_monthly_shortname(get_field_variable('Ek_vel')), get_monthly_nc_string(get_field_variable('Ek_vel'))
        if not os.path.exists(join(datdir_secondary, vel_Ek_shortname)):
            os.makedirs(join(datdir_secondary, vel_Ek_shortname))
        vel_Ek_ds.to_netcdf(path=join(datdir_secondary, vel_Ek_shortname, vel_Ek_nc_string+date_string+".nc"), engine="scipy")
            
##############################

#each field gets a separate function
    #some of them will call other functions, as appropriate

#horizontal velocity div

#normal strain

#shear strain

#vorticity

#Ekman velocity divergence
    
##############################

def main(**kwargs):
    
    if kwargs:
        
        datdir_primary = kwargs.get('datdir_primary')
        datdir_secondary = kwargs.get('datdir_secondary')
        
        monthstr = kwargs.get('monthstr')
        yearstr = kwargs.get('yearstr')
        
        field_name = kwargs.get('field_name')
        
        rho_ref = float(kwargs.get('rho_ref'))
        nu_E = float(kwargs.get('nu_E'))
    
    ds_grid = load_grid(datdir_primary) #Load the grid DataSet
    
    #Identify the input field and run the appropriate function to compute it

    if field_name == 'geostrophic_vel':
        comp_geostrophic_vel(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, rho_ref, time_ave_type='monthly')
    elif field_name == 'Ek_vel':
        comp_Ek_vel(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, rho_ref, nu_E, time_ave_type='monthly')

##############################

if __name__ == "__main__":
    main()