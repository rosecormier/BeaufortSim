"""
Contains necessary auxiliary functions and functions that directly produce computed data.
    -Divergence of horizontal velocity;
    -Geostrophic velocity;
    -Ekman velocity;
    -Normal strain;
    -Shear strain;
    -Vorticity.
    
The main function picks out the appropriate function to compute the requested field.

Notes:
    -Only does monthly averaging at the moment; will update this.
    -Does not currently include divergence of Ekman velocity; could add this, but is a nontrivial process.
        -There's an old script that computes this field in a format suitable for plotting, if needed.
        -Okubo-Weiss computations are also in an old script.

Rosalie Cormier, 2023
"""

import os
import xarray as xr
import numpy as np
import ecco_v4_py as ecco

from os.path import join
from math import e, pi

import download_data

from functions_ecco_general import load_ECCO_data_file, load_grid
from functions_field_variables import get_field_variable, get_monthly_shortname, get_monthly_nc_string

##############################

#AUXILIARY VALUES/FUNCTIONS; USED DURING SUBSEQUENT COMPUTATIONS

Omega = (2 * pi) / 86164 #Earth angular velocity

def to_radians(angle): #Converts degrees to radians
    return angle * pi / 180

def comp_f(y): #Computes local Coriolis parameter
    lat = to_radians(y)
    return 2 * Omega * np.sin(lat)

#Gets density and pressure-like attributes from a DataSet, given a reference density
def get_density_and_pressure(ds_denspress, rho_ref):

    density_anom = ds_denspress['RHOAnoma']
    density = density_anom + rho_ref #Compute absolute density
    
    #Update density DataArray
    
    density.name = 'RHO'
    density.attrs.update({'long_name': 'In-situ seawater density', 'units': 'kg m-3'})
    
    pressure_anom = ds_denspress.PHIHYDcR.copy() #Copy pressure data
    pressure_like = rho_ref * pressure_anom #Compute quantity to differentiate (units of density * pressure)
    
    return density, pressure_like

#Differentiates pressure-like quantity to compute geostrophic-velocity components
def get_vel_g_components(ds_grid, pressure_like, density):
    
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

#Computes Ekman-velocity components, using surface density, surface stress, Ekman-layer depth
def get_vel_Ek_components(ds_grid, ds_denspress, ds_stress, rho_ref, nu_E):
    
    z = ds_grid.Z #Get z-coordinate
    
    density = get_density_and_pressure(ds_denspress, rho_ref)[0] #Get density field
    surface_density = density.isel(k=0) #Compute surface density field
    
    f = comp_f(ds_grid.YC) #Compute Coriolis param. from latitudes of grid cell centres
    Ek_depth = (2 * nu_E / f)**0.5 #Compute Ekman-layer depth
    
    #Get wind-stress components (x-y coordinates) and rotate wind-stress vector (to E-N coordinates)
    
    tau_x, tau_y = ds_stress.EXFtaux, ds_stress.EXFtauy
    tau_E = tau_x * ds_grid['CS'] - tau_y * ds_grid['SN']
    tau_N = tau_x * ds_grid['SN'] + tau_y * ds_grid['CS']
    
    Ek_coeff = (np.sqrt(2) / (surface_density * f * Ek_depth)) #Coefficient that multiplies Ekman velocity
    
    #Compute Ekman velocity components
    
    u_Ek = Ek_coeff * e**(z/Ek_depth) * (tau_E * np.cos((z/Ek_depth) - (pi/4)) - tau_N * np.sin((z/Ek_depth) - (pi/4)))
    v_Ek = Ek_coeff * e**(z/Ek_depth) * (tau_E * np.sin((z/Ek_depth) - (pi/4)) + tau_N * np.cos((z/Ek_depth) - (pi/4)))
    u_Ek, v_Ek = u_Ek.where(ds_grid.maskC), v_Ek.where(ds_grid.maskC) #Mask land on result
    
    return u_Ek, v_Ek

##############################

#DIVERGENCE OF HORIZONTAL VELOCITY

def comp_2D_div_vel(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, time_ave_type='monthly'):
    
    """
    Computes and saves divergence of horizontal velocity to DataSet in NetCDF format.
    """
    
    if time_ave_type == 'monthly':
        
        #Look for the velocity file in primary directory and download it if it doesn't exist
        download_data.main(field_name='horizontal_vel', initial_month=monthstr, initial_year=yearstr, final_month=monthstr, final_year=yearstr, time_ave_type=time_ave_type, datdir_primary=datdir_primary)
        
        date_string = yearstr + '-' + monthstr
        
        ds_velocity = load_ECCO_data_file('horizontal_vel', date_string, datdir_primary, time_ave_type) #Load DataSet
        ds_velocity['UVEL'].data, ds_velocity['VVEL'].data = ds_velocity['UVEL'].values, ds_velocity['VVEL'].values
    
        xgcm_grid = ecco.get_llc_grid(ds_grid)

        u, v = ds_velocity['UVEL'], ds_velocity['VVEL']
        dx, dy = ds_grid.dxG, ds_grid.dyG
        cell_area = ds_grid.rA

        div_vel = (xgcm_grid.diff(u*dy, 'X') + xgcm_grid.diff(v*dx, 'Y')).squeeze() / cell_area #Compute divergence

        #Save the data

        div_vel.name = 'DIVU'
        div_vel_shortname, div_vel_nc_string = get_monthly_shortname(get_field_variable('2D_div_vel')), get_monthly_nc_string(get_field_variable('2D_div_vel'))
        if not os.path.exists(join(datdir_secondary, div_vel_shortname)):
            os.makedirs(join(datdir_secondary, div_vel_shortname))
        div_vel.to_netcdf(path=join(datdir_secondary, div_vel_shortname, div_vel_nc_string+date_string+".nc"), engine="scipy")

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
        
        ds_denspress = load_ECCO_data_file('density', date_string, datdir_primary, time_ave_type) #Load DataSet
        
        #Extract density and pressure-like from the DataSet
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
        
        ds_denspress = load_ECCO_data_file('density', date_string, datdir_primary, time_ave_type) #Load DataSet
       
        #Look for the surface wind-on-ocean-stress file in primary directory and download if it doesn't exist
        download_data.main(field_name='wind_stress', initial_month=monthstr, initial_year=yearstr, final_month=monthstr, final_year=yearstr, time_ave_type=time_ave_type, datdir_primary=datdir_primary)
        
        ds_stress = load_ECCO_data_file('wind_stress', date_string, datdir_primary, time_ave_type) #Load DataSet
        
        #Compute Ekman velocity components
        u_Ek, v_Ek = get_vel_Ek_components(ds_grid, ds_denspress, ds_stress, rho_ref, nu_E)

        #Save the data
        
        u_Ek.name = 'UEk'
        v_Ek.name = 'VEk'
        vel_Ek_ds = xr.merge([u_Ek, v_Ek])
        vel_Ek_shortname, vel_Ek_nc_string = get_monthly_shortname(get_field_variable('Ek_vel')), get_monthly_nc_string(get_field_variable('Ek_vel'))
        if not os.path.exists(join(datdir_secondary, vel_Ek_shortname)):
            os.makedirs(join(datdir_secondary, vel_Ek_shortname))
        vel_Ek_ds.to_netcdf(path=join(datdir_secondary, vel_Ek_shortname, vel_Ek_nc_string+date_string+".nc"), engine="scipy")
        
##############################

#NORMAL STRAIN

def comp_normal_strain(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, time_ave_type='monthly'):
    
    """
    Computes and saves normal strain to DataSet in NetCDF format.
    """
    
    if time_ave_type == 'monthly':
        
        #Look for the velocity file in primary directory and download it if it doesn't exist
        download_data.main(field_name='horizontal_vel', initial_month=monthstr, initial_year=yearstr, final_month=monthstr, final_year=yearstr, time_ave_type=time_ave_type, datdir_primary=datdir_primary)
        
        date_string = yearstr + '-' + monthstr
        
        ds_velocity = load_ECCO_data_file('horizontal_vel', date_string, datdir_primary, time_ave_type) #Load DataSet
        
        ds_velocity['UVEL'].data, ds_velocity['VVEL'].data = ds_velocity['UVEL'].values, ds_velocity['VVEL'].values

        xgcm_grid = ecco.get_llc_grid(ds_grid)
        
        u_mean, v_mean = ds_velocity['UVEL'], ds_velocity['VVEL']
        dx, dy = ds_grid.dxG, ds_grid.dyG
        cell_area = ds_grid.rA
 
        normal_strain = (xgcm_grid.diff(u_mean*dy, 'X') - xgcm_grid.diff(v_mean*dx, 'Y')) / cell_area #Compute strain
        
        #Save the data

        normal_strain.name = 'NORMAL'
        normal_shortname, normal_nc_string = get_monthly_shortname(get_field_variable('normal_strain')), get_monthly_nc_string(get_field_variable('normal_strain'))
        if not os.path.exists(join(datdir_secondary, normal_shortname)):
            os.makedirs(join(datdir_secondary, normal_shortname))
        normal_strain.to_netcdf(path=join(datdir_secondary, normal_shortname, normal_nc_string+date_string+".nc"), engine="scipy")

##############################

#SHEAR STRAIN

def comp_shear_strain(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, time_ave_type='monthly'):
    
    """
    Computes and saves shear strain to DataSet in NetCDF format.
    """
    
    if time_ave_type == 'monthly':
        
        #Look for the velocity file in primary directory and download it if it doesn't exist
        download_data.main(field_name='horizontal_vel', initial_month=monthstr, initial_year=yearstr, final_month=monthstr, final_year=yearstr, time_ave_type=time_ave_type, datdir_primary=datdir_primary)
        
        date_string = yearstr + '-' + monthstr
        
        ds_velocity = load_ECCO_data_file('horizontal_vel', date_string, datdir_primary, time_ave_type) #Load DataSet
        ds_velocity['UVEL'].data, ds_velocity['VVEL'].data = ds_velocity['UVEL'].values, ds_velocity['VVEL'].values

        xgcm_grid = ecco.get_llc_grid(ds_grid)
        
        u_mean, v_mean = ds_velocity['UVEL'], ds_velocity['VVEL']
        dx, dy = ds_grid.dxC, ds_grid.dyC
        cell_area = ds_grid.rAz
        
        shear_strain = (xgcm_grid.diff(v_mean*dy, 'X') + xgcm_grid.diff(u_mean*dx, 'Y')) / cell_area #Compute strain
        
        #Save the data

        shear_strain.name = 'SHEAR'
        shear_shortname, shear_nc_string = get_monthly_shortname(get_field_variable('shear_strain')), get_monthly_nc_string(get_field_variable('shear_strain'))
        if not os.path.exists(join(datdir_secondary, shear_shortname)):
            os.makedirs(join(datdir_secondary, shear_shortname))
        shear_strain.to_netcdf(path=join(datdir_secondary, shear_shortname, shear_nc_string+date_string+".nc"), engine="scipy")

##############################

#VORTICITY

def comp_vorticity(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, time_ave_type='monthly'):
    
    """
    Computes and saves vorticity to DataSet in NetCDF format.
    """
    
    if time_ave_type == 'monthly':
        
        #Look for the velocity file in primary directory and download it if it doesn't exist
        download_data.main(field_name='horizontal_vel', initial_month=monthstr, initial_year=yearstr, final_month=monthstr, final_year=yearstr, time_ave_type=time_ave_type, datdir_primary=datdir_primary)
        
        date_string = yearstr + '-' + monthstr
        
        ds_velocity = load_ECCO_data_file('horizontal_vel', date_string, datdir_primary, time_ave_type) #Load DataSet
        ds_velocity['UVEL'].data, ds_velocity['VVEL'].data = ds_velocity['UVEL'].values, ds_velocity['VVEL'].values

        xgcm_grid = ecco.get_llc_grid(ds_grid)
        
        u_mean, v_mean = ds_velocity['UVEL'], ds_velocity['VVEL']
        dx, dy = ds_grid.dxC, ds_grid.dyC
        cell_area = ds_grid.rAz

        vorticity = (-xgcm_grid.diff(u_mean*dx, 'Y') + xgcm_grid.diff(v_mean*dy, 'X')).squeeze() / cell_area #Compute vorticity

        #Save the data
    
        vorticity.name = 'ZETA'
        vorticity_shortname, vorticity_nc_string = get_monthly_shortname(get_field_variable('vorticity')), get_monthly_nc_string(get_field_variable('vorticity'))
        if not os.path.exists(join(datdir_secondary, vorticity_shortname)):
            os.makedirs(join(datdir_secondary, vorticity_shortname))
        vorticity.to_netcdf(path=join(datdir_secondary, vorticity_shortname, vorticity_nc_string+date_string+".nc"), engine="scipy") 
        
##############################

def main(**kwargs):
    
    if kwargs:
        
        datdir_primary = kwargs.get('datdir_primary')
        datdir_secondary = kwargs.get('datdir_secondary')
        
        time_ave_type = kwargs.get('time_ave_type')
        monthstr = kwargs.get('monthstr')
        yearstr = kwargs.get('yearstr')
        
        field_name = kwargs.get('field_name')
        
        rho_ref = float(kwargs.get('rho_ref'))
        nu_E = float(kwargs.get('nu_E'))
    
    ds_grid = load_grid(datdir_primary) #Load the grid DataSet
    
    #Identify the input field and run the appropriate function to compute it

    if field_name == '2D_div_vel':
        comp_2D_div_vel(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, time_ave_type=time_ave_type)
    elif field_name == 'geostrophic_vel':
        comp_geostrophic_vel(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, rho_ref, time_ave_type=time_ave_type)
    elif field_name == 'Ek_vel':
        comp_Ek_vel(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, rho_ref, nu_E, time_ave_type=time_ave_type)
    elif field_name == 'normal_strain':
        comp_normal_strain(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, time_ave_type=time_ave_type)
    elif field_name == 'shear_strain':
        comp_shear_strain(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, time_ave_type=time_ave_type)
    elif field_name == 'vorticity':
        comp_vorticity(ds_grid, monthstr, yearstr, datdir_primary, datdir_secondary, time_ave_type=time_ave_type)

##############################

if __name__ == "__main__":
    main()