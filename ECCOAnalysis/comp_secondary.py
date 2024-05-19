"""
Contains necessary auxiliary functions and functions that directly produce 
secondary data.
    -Divergence of horizontal velocity;
    -Geostrophic velocity;
    -Ekman velocity;
    -Normal strain;
    -Shear strain;
    -Vorticity;
    -Vertical velocity times potential temperature.
The main function checks whether computed data already exist and computes them 
if not.

Notes:
    -Does not currently include divergence of Ekman velocity; could add this, 
    but is a nontrivial process.
        -There's an old script that computes this field in a format suitable for
        plotting, if needed.
        -Okubo-Weiss computations are also in an old script.

R. Cormier, 2024
"""

import os
import xarray as xr
import numpy as np
import ecco_v4_py as ecco

from math import e, pi
from os.path import join

import download_data
import load_data_files

from functions_ecco_general import load_grid, get_args_from_date_string, \
get_monthstr
from functions_field_variables import get_variable_str, get_monthly_shortname, \ 
get_monthly_nc_string, get_seasonal_shortname, get_seasonal_nc_string

##############################

#AUXILIARY VALUES/FUNCTIONS; USED DURING SUBSEQUENT COMPUTATIONS

Omega = (2 * pi) / 86164 #Earth angular velocity

def to_radians(angle):
    return angle * pi / 180

def comp_f(y): #Computes local Coriolis parameter
    lat = to_radians(y)
    return 2 * Omega * np.sin(lat)

def get_density_and_pressure(ds_denspress, rho_ref):
    """
    Gets density and pressure-like attributes from a DataSet, given a reference 
    density.
    """

    density_anom = ds_denspress['RHOAnoma']
    density = density_anom + rho_ref #Compute absolute density
    
    #Update density DataArray
    density.name = 'RHO'
    density.attrs.update({'long_name': 'In-situ seawater density', 
                          'units': 'kg m-3'})
    
    pressure_anom = ds_denspress.PHIHYDcR.copy() #Copy pressure data
    
    #Compute quantity to differentiate (units of density * pressure)
    pressure_like = rho_ref * pressure_anom 
    
    return density, pressure_like

def get_vel_g_components(ds_grid, pressure_like, density):
    """
    Differentiates pressure-like quantity to compute geostrophic-velocity 
    components.
    """
    
    xgcm_grid = ecco.get_llc_grid(ds_grid)
    
    #Compute derivatives of pressure in x and y
    d_press_dx = ((xgcm_grid.diff(pressure_like, axis="X", boundary='extend')) / 
                  ds_grid.dxC)
    d_press_dy = ((xgcm_grid.diff(pressure_like, axis="Y", boundary='extend')) /
                  ds_grid.dyC)
    
    #Convert DataArray content from dask to np arrays
    d_press_dx.data, d_press_dy.data = d_press_dx.values, d_press_dy.values
    
    #Interpolate to centres of grid cells
    press_grads_interp = xgcm_grid.interp_2d_vector({'X': d_press_dx, 
                                                     'Y': d_press_dy}, 
                                                    boundary='extend')
    dp_dx, dp_dy = press_grads_interp['X'], press_grads_interp['Y']
    
    dp_dx.name = 'dp_dx'
    dp_dy.name = 'dp_dy'
    
    #Compute RHS of geostrophic-balance equations; mask land
    GB_RHS_1, GB_RHS_2 = dp_dx / density, - dp_dy / density 
    GB_RHS_1 = GB_RHS_1.where(ds_grid.maskC)
    GB_RHS_2 = GB_RHS_2.where(ds_grid.maskC)
    
    #Compute Coriolis parameter from latitudes of grid cell centres
    f = comp_f(ds_grid.YC) 
    #Compute geostrophic-velocity components
    u_g, v_g = GB_RHS_2 / f, GB_RHS_1 / f
    #Mask land on result
    u_g, v_g = u_g.where(ds_grid.maskC), v_g.where(ds_grid.maskC) 
    
    return u_g, v_g

def get_vel_Ek_components(ds_grid, ds_denspress, ds_stress, rho_ref, nu_E):
    """
    Computes Ekman-velocity components, using surface density, surface stress, 
    Ekman-layer depth.
    """
    
    z = ds_grid.Z
    
    density = get_density_and_pressure(ds_denspress, rho_ref)[0] 
    surface_density = density.isel(k=0) 
    
    #Compute Coriolis parameter from latitudes of grid cell centres
    f = comp_f(ds_grid.YC) 
    Ek_depth = (2 * nu_E / f)**0.5 #Ekman-layer depth
    
    #Get wind-stress components (x-y coordinates)
    tau_x, tau_y = ds_stress.EXFtaux, ds_stress.EXFtauy
    
    #Rotate wind-stress vector (to E-N coordinates)
    tau_E = tau_x * ds_grid['CS'] - tau_y * ds_grid['SN']
    tau_N = tau_x * ds_grid['SN'] + tau_y * ds_grid['CS']
    
    #Coefficient that multiplies Ekman velocity
    Ek_coeff = (np.sqrt(2) / (surface_density * f * Ek_depth)) 
    
    #Compute Ekman velocity component
    u_Ek = (Ek_coeff * e**(z/Ek_depth) * (tau_E * np.cos((z/Ek_depth) - (pi/4)) 
                                          - tau_N * np.sin((z/Ek_depth) 
                                                           - (pi/4))))
    v_Ek = (Ek_coeff * e**(z/Ek_depth) * (tau_E * np.sin((z/Ek_depth) - (pi/4)) 
                                          + tau_N * np.cos((z/Ek_depth) 
                                                           - (pi/4))))
    #Mask land on result
    u_Ek, v_Ek = u_Ek.where(ds_grid.maskC), v_Ek.where(ds_grid.maskC) 
    
    return u_Ek, v_Ek

##############################

def get_time_ave_type_attrs(date_string, time_ave_type, time_kwargs, 
                            field_name):
    """
    Gets variables associated with the time-average type specified.
    Note - only handles time_ave_type == 'monthly' or 'seasonal'.
    """
    
    args = get_args_from_date_string(date_string, time_ave_type, time_kwargs)
    initial_month, initial_year = args[0], args[1]
    final_month, final_year = args[2], args[3]
    
    if time_ave_type == 'monthly':
        month_strings = [date_string]
        
    elif time_ave_type == 'seasonal':
        
        month_strings = []
        month, year = int(initial_month), int(initial_year)
        
        while year < int(final_year):
            while month <= 12:
                
                month_string = str(year) + '-' + get_monthstr(month)
                month_strings.append(month_string)
                month += 1
                
                if month == int(final_month) + 1:
                    break
                if month == 13:
                    year += 1
                    month = 1
                    
        if year == int(final_year):
            while month <= int(final_month):
                month_string = str(year) + '-' + get_monthstr(month)
                month_strings.append(month_string)
                month += 1
                
    month_string_dict = {}
    for month_string in month_strings:
        month_string_dict[month_string] = [month_string[5:7], month_string[0:4]]
            #This creates a dictionary whose keys are of the form 'YYYY-MM' and 
            #whose values are of the form ['MM', 'YYYY'].
            #This is useful in the loops used in the following functions.
        
    if time_ave_type == 'monthly':
        var_shortname = get_monthly_shortname(get_variable_str(field_name))
        var_nc_string = get_monthly_nc_string(get_variable_str(field_name))
    elif time_ave_type == 'seasonal':
        var_shortname = get_seasonal_shortname(get_variable_str(field_name))
        var_nc_string = get_seasonal_nc_string(get_variable_str(field_name))
    
    return month_string_dict, var_shortname, var_nc_string

##############################

#GEOSTROPHIC VELOCITY

def comp_geostrophic_vel(ds_grid, date_string, datdir_primary, datdir_secondary,
                        time_ave_type, time_kwargs, rho_ref):
    """
    Compute and save geostrophic velocity components (u_g, v_g) to DataSet in 
    NetCDF format.
    """
    
    time_ave_type_attrs = get_time_ave_type_attrs(date_string, time_ave_type, 
                                                  time_kwargs, 
                                                  'geostrophic_vel')
    month_string_dict = time_ave_type_attrs[0]
    monthly_file_list = []
    
    for month_string in month_string_dict:
        
        month = month_string_dict[month_string][0]
        year = month_string_dict[month_string][1]
        
        vel_g_shortname = get_monthly_shortname('UGVG')
        vel_g_nc_string = get_monthly_nc_string('UGVG')
    
        #Look for dens/press file in primary directory and download if needed
        download_data.main(field_name='density_anom', initial_month=month, 
                           initial_year=year, final_month=month, 
                           final_year=year, time_ave_type='monthly', 
                           datdir_primary=datdir_primary, time_kwargs=None)
        #Load monthly DataSet
        ds_denspress = load_data_files.main(field_name='density_anom', 
                                            time_ave_type='monthly', 
                                            date_string=month_string, 
                                            datdir_primary=datdir_primary, 
                                            datdir_secondary=datdir_secondary) 

        #Extract density and pressure-like from the DataSet
        density, pressure_like = get_density_and_pressure(ds_denspress, rho_ref)

        #Compute geostrophic velocity components
        u_g, v_g = get_vel_g_components(ds_grid, pressure_like, density) 
        u_g, v_g = u_g.squeeze(), v_g.squeeze()
        
        #Save the monthly data
        u_g.name = 'UG'
        v_g.name = 'VG'
        vel_g_ds = xr.merge([u_g, v_g], combine_attrs='no_conflicts')
        if not os.path.exists(join(datdir_secondary, vel_g_shortname)):
            os.makedirs(join(datdir_secondary, vel_g_shortname))
        file_path = join(datdir_secondary, vel_g_shortname,
                        vel_g_nc_string+month_string+".nc")
        vel_g_ds.to_netcdf(path=file_path, engine="scipy")
        monthly_file_list.append(file_path)

    #If seasonal, compute the seasonally-averaged data if missing
    if time_ave_type == 'seasonal':
        
        vel_g_shortname = time_ave_type_attrs[1]
        vel_g_nc_string = time_ave_type_attrs[2]
        
        if not os.path.exists(join(datdir_secondary, vel_g_shortname)):
            os.makedirs(join(datdir_secondary, vel_g_shortname))
        
        vel_g_seasonal = xr.open_mfdataset(monthly_file_list, 
                                               combine='nested', 
                                               concat_dim='month')
        vel_g_seasonal_avg = (vel_g_seasonal.sum('month') 
                                  / len(monthly_file_list))
        vel_g_seasonal_avg.to_netcdf(path=join(datdir_secondary, 
                                        vel_g_shortname, 
                                        vel_g_nc_string+date_string+".nc"), 
                                         engine="scipy")

############################## 

#RATE OF VERTICAL TEMPERATURE TRANSPORT

def comp_T_zflux(ds_grid, date_string, datdir_primary, datdir_secondary,
                 time_ave_type, time_kwargs):
    
    xgcm_grid = ecco.get_llc_grid(ds_grid)
    
    time_ave_type_attrs = get_time_ave_type_attrs(date_string, time_ave_type, 
                                                  time_kwargs, 'T_zflux')
    month_string_dict = time_ave_type_attrs[0]
    monthly_file_list = []
    
    for month_string in month_string_dict:
        
        month = month_string_dict[month_string][0]
        year = month_string_dict[month_string][1]
        
        T_ZFLUX_shortname = get_monthly_shortname('T_ZFLUX')
        T_ZFLUX_nc_string = get_monthly_nc_string('T_ZFLUX')

        #Look for theta-file in primary directory and download if needed
        download_data.main(field_name='potential_temp', initial_month=month, 
                               initial_year=year, final_month=month, 
                               final_year=year, time_ave_type='monthly', 
                               datdir_primary=datdir_primary, time_kwargs=None)
        #Load monthly DataSet
        ds_temp = load_data_files.main(field_name='potential_temp', 
                                            time_ave_type='monthly', 
                                            date_string=month_string, 
                                            datdir_primary=datdir_primary, 
                                            datdir_secondary=datdir_secondary)
        
        #Look for velocity file in primary directory and download if needed
        download_data.main(field_name='vertical_vel', initial_month=month, 
                               initial_year=year, final_month=month, 
                               final_year=year, time_ave_type='monthly', 
                               datdir_primary=datdir_primary, time_kwargs=None)
        #Load monthly DataSet
        ds_vel = load_data_files.main(field_name='vertical_vel', 
                                            time_ave_type='monthly', 
                                            date_string=month_string, 
                                            datdir_primary=datdir_primary, 
                                            datdir_secondary=datdir_secondary)
        
        wvel = xgcm_grid.interp(ds_vel.WVEL.values, axis="Z").squeeze()
        wT = ds_temp.THETA.values.squeeze() * wvel
        wT_integ = wT.sum('k')
        
        wT_integ.name = 'T_ZFLUX'
        
        #Save the monthly data
        if not os.path.exists(join(datdir_secondary, T_ZFLUX_shortname)):
            os.makedirs(join(datdir_secondary, T_ZFLUX_shortname))
        file_path = join(datdir_secondary, T_ZFLUX_shortname,
                        T_ZFLUX_nc_string+month_string+".nc")
        wT_integ.to_netcdf(path=file_path, engine="scipy")
        monthly_file_list.append(file_path)
        
    #If seasonal, compute the seasonally-averaged data if missing
    if time_ave_type == 'seasonal':
        
        T_ZFLUX_shortname = time_ave_type_attrs[1]
        T_ZFLUX_nc_string = time_ave_type_attrs[2]
        
        if not os.path.exists(join(datdir_secondary, T_ZFLUX_shortname)):
            os.makedirs(join(datdir_secondary, T_ZFLUX_shortname))
        
        T_ZFLUX_seasonal = xr.open_mfdataset(monthly_file_list, 
                                               combine='nested', 
                                               concat_dim='month')
        T_ZFLUX_seasonal_avg = (T_ZFLUX_seasonal.sum('month') 
                                  / len(monthly_file_list))
        T_ZFLUX_seasonal_avg.to_netcdf(path=join(datdir_secondary, 
                                        T_ZFLUX_shortname, 
                                        T_ZFLUX_nc_string+date_string+".nc"), 
                                         engine="scipy")
        
############################## 

#EKMAN VELOCITY

def comp_Ek_vel(ds_grid, date_string, datdir_primary, datdir_secondary, 
                time_ave_type, time_kwargs, rho_ref, nu_E):
    """
    Compute and save Ekman velocity components (u_Ek, v_Ek) to DataSet in NetCDF
    format.
    """
    
    time_ave_type_attrs = get_time_ave_type_attrs(date_string, time_ave_type, 
                                                  time_kwargs, 'Ek_vel')
    month_string_dict = time_ave_type_attrs[0]
    monthly_file_list = []
    
    for month_string in month_string_dict:
        
        month = month_string_dict[month_string][0]
        year = month_string_dict[month_string][1]
        
        vel_Ek_shortname = get_monthly_shortname('UEkVEk')
        vel_Ek_nc_string = get_monthly_nc_string('UEkVEk')
    
        #Look for dens/press file in primary directory and download if needed
        download_data.main(field_name='density_anom', initial_month=month, 
                           initial_year=year, final_month=month, 
                           final_year=year, time_ave_type='monthly', 
                           datdir_primary=datdir_primary, time_kwargs=None)
        #Load monthly DataSet
        ds_denspress = load_data_files.main(field_name='density_anom', 
                                            time_ave_type='monthly', 
                                            date_string=month_string, 
                                            datdir_primary=datdir_primary, 
                                            datdir_secondary=datdir_secondary) 

        #Look for surface-stress file in primary directory; download if needed
        download_data.main(field_name='wind_stress', initial_month=month, 
                           initial_year=year, final_month=month, 
                           final_year=year, time_ave_type='monthly', 
                           datdir_primary=datdir_primary, time_kwargs=None)
        #Load monthly DataSet
        ds_stress = load_data_files.main(field_name='wind_stress', 
                                         time_ave_type='monthly', 
                                         date_string=month_string, 
                                         datdir_primary=datdir_primary, 
                                         datdir_secondary=datdir_secondary) 
        
        #Compute Ekman velocity components
        u_Ek, v_Ek = get_vel_Ek_components(ds_grid, ds_denspress, ds_stress, 
                                           rho_ref, nu_E) 
        u_Ek, v_Ek = u_Ek.squeeze(), v_Ek.squeeze()
        
        u_Ek.name = 'UEk'
        v_Ek.name = 'VEk'
        
        #Save the monthly data
        vel_Ek_ds = xr.merge([u_Ek, v_Ek], combine_attrs='no_conflicts')
        if not os.path.exists(join(datdir_secondary, vel_Ek_shortname)):
            os.makedirs(join(datdir_secondary, vel_Ek_shortname))
        file_path = join(datdir_secondary, vel_Ek_shortname,
                        vel_Ek_nc_string+month_string+".nc")
        vel_Ek_ds.to_netcdf(path=file_path, engine="scipy")
        monthly_file_list.append(file_path)

    #If seasonal, compute the seasonally-averaged data if missing
    if time_ave_type == 'seasonal':
        
        vel_Ek_shortname = time_ave_type_attrs[1]
        vel_Ek_nc_string = time_ave_type_attrs[2]
        
        if not os.path.exists(join(datdir_secondary, vel_Ek_shortname)):
            os.makedirs(join(datdir_secondary, vel_Ek_shortname))
        
        vel_Ek_seasonal = xr.open_mfdataset(monthly_file_list, 
                                               combine='nested', 
                                               concat_dim='month')
        vel_Ek_seasonal_avg = (vel_Ek_seasonal.sum('month') 
                                  / len(monthly_file_list))
        vel_Ek_seasonal_avg.to_netcdf(path=join(datdir_secondary, 
                                        vel_Ek_shortname, 
                                        vel_Ek_nc_string+date_string+".nc"), 
                                         engine="scipy")
        
##############################

#DIVERGENCE OF HORIZONTAL VELOCITY

def comp_2D_div_vel(ds_grid, date_string, datdir_primary, datdir_secondary, 
                    time_ave_type, time_kwargs):
    """
    Compute and save divergence of horizontal velocity to DataSet in NetCDF 
    format.
    """
    
    time_ave_type_attrs = get_time_ave_type_attrs(date_string, time_ave_type, 
                                                  time_kwargs, '2D_div_vel')
    month_string_dict = time_ave_type_attrs[0]
    monthly_file_list = []
    
    for month_string in month_string_dict:
        
        month = month_string_dict[month_string][0] 
        year = month_string_dict[month_string][1]
        
        div_vel_shortname = get_monthly_shortname('DIVU')
        div_vel_nc_string = get_monthly_nc_string('DIVU')
        
        #Look for velocity file in primary directory and download if needed
        download_data.main(field_name='horizontal_vel', initial_month=month, 
                           initial_year=year, final_month=month, 
                           final_year=year, time_ave_type='monthly', 
                           datdir_primary=datdir_primary, time_kwargs=None)
        #Load monthly DataSet
        ds_velocity = load_data_files.main(field_name='horizontal_vel', 
                                           time_ave_type='monthly', 
                                           date_string=month_string, 
                                           datdir_primary=datdir_primary, 
                                           datdir_secondary=datdir_secondary)
        ds_velocity['UVEL'].data = ds_velocity['UVEL'].values
        ds_velocity['VVEL'].data = ds_velocity['VVEL'].values

        xgcm_grid = ecco.get_llc_grid(ds_grid)

        u, v = ds_velocity['UVEL'], ds_velocity['VVEL']
        dx, dy = ds_grid.dxG, ds_grid.dyG
        cell_area = ds_grid.rA

        #Compute divergence
        div_vel = ((xgcm_grid.diff(u*dy, 'X') 
                    + xgcm_grid.diff(v*dx, 'Y')).squeeze() / cell_area) 
        
        div_vel.name = 'DIVU'
        
        #Save the monthly data
        if not os.path.exists(join(datdir_secondary, div_vel_shortname)):
            os.makedirs(join(datdir_secondary, div_vel_shortname))
        file_path = join(datdir_secondary, div_vel_shortname,
                        div_vel_nc_string+month_string+".nc")
        div_vel.to_netcdf(path=file_path, engine="scipy")
        monthly_file_list.append(file_path)

    #If seasonal, compute the seasonally-averaged data if missing
    if time_ave_type == 'seasonal':
        
        div_vel_shortname = time_ave_type_attrs[1]
        div_vel_nc_string = time_ave_type_attrs[2]
        
        if not os.path.exists(join(datdir_secondary, div_vel_shortname)):
            os.makedirs(join(datdir_secondary, div_vel_shortname))
        
        div_vel_seasonal = xr.open_mfdataset(monthly_file_list, 
                                               combine='nested', 
                                               concat_dim='month')
        div_vel_seasonal_avg = (div_vel_seasonal.sum('month') 
                                  / len(monthly_file_list))
        div_vel_seasonal_avg.to_netcdf(path=join(datdir_secondary, 
                                        div_vel_shortname, 
                                        div_vel_nc_string+date_string+".nc"), 
                                         engine="scipy")

##############################

#NORMAL STRAIN

def comp_normal_strain(ds_grid, date_string, datdir_primary, datdir_secondary, 
                       time_ave_type, time_kwargs):
    """
    Compute and save normal strain to DataSet in NetCDF format.
    """
    
    time_ave_type_attrs = get_time_ave_type_attrs(date_string, time_ave_type, 
                                                  time_kwargs, 'normal_strain')
    month_string_dict = time_ave_type_attrs[0]
    monthly_file_list = []
    
    for month_string in month_string_dict:
        
        month = month_string_dict[month_string][0]
        year = month_string_dict[month_string][1]
        
        normal_shortname = get_monthly_shortname('NORMAL')
        normal_nc_string = get_monthly_nc_string('NORMAL')
        
        #Look for velocity file in primary directory and download if needed
        download_data.main(field_name='horizontal_vel', initial_month=month, 
                           initial_year=year, final_month=month, 
                           final_year=year, time_ave_type='monthly', 
                           datdir_primary=datdir_primary, time_kwargs=None)
        #Load monthly DataSet
        ds_velocity = load_data_files.main(field_name='horizontal_vel', 
                                           time_ave_type='monthly', 
                                           date_string=month_string, 
                                           datdir_primary=datdir_primary,
                                           datdir_secondary=datdir_secondary)
        ds_velocity['UVEL'].data = ds_velocity['UVEL'].values
        ds_velocity['VVEL'].data = ds_velocity['VVEL'].values

        xgcm_grid = ecco.get_llc_grid(ds_grid)

        u_mean, v_mean = ds_velocity['UVEL'], ds_velocity['VVEL']
        dx, dy = ds_grid.dxG, ds_grid.dyG
        cell_area = ds_grid.rA

        #Compute strain
        normal_strain = ((xgcm_grid.diff(u_mean*dy, 'X') - 
                          xgcm_grid.diff(v_mean*dx, 'Y')).squeeze() / cell_area)
        
        normal_strain.name = 'NORMAL'
        
        #Save the monthly data
        if not os.path.exists(join(datdir_secondary, normal_shortname)):
            os.makedirs(join(datdir_secondary, normal_shortname))
        file_path = join(datdir_secondary, normal_shortname,
                        normal_nc_string+month_string+".nc")
        normal_strain.to_netcdf(path=file_path, engine="scipy")
        monthly_file_list.append(file_path)
        
    #If seasonal, compute the seasonally-averaged data if missing
    if time_ave_type == 'seasonal':
        
        normal_shortname = time_ave_type_attrs[1]
        normal_nc_string = time_ave_type_attrs[2]
        
        if not os.path.exists(join(datdir_secondary, normal_shortname)):
            os.makedirs(join(datdir_secondary, normal_shortname))
        
        normal_seasonal = xr.open_mfdataset(monthly_file_list, 
                                               combine='nested', 
                                               concat_dim='month')
        normal_seasonal_avg = (normal_seasonal.sum('month') 
                                  / len(monthly_file_list))
        normal_seasonal_avg.to_netcdf(path=join(datdir_secondary, 
                                        normal_shortname, 
                                        normal_nc_string+date_string+".nc"), 
                                         engine="scipy")

##############################

#SHEAR STRAIN

def comp_shear_strain(ds_grid, date_string, datdir_primary, datdir_secondary, 
                      time_ave_type, time_kwargs):
    """
    Compute and save shear strain to DataSet in NetCDF format.
    """
    
    time_ave_type_attrs = get_time_ave_type_attrs(date_string, time_ave_type, 
                                                  time_kwargs, 'shear_strain')
    month_string_dict = time_ave_type_attrs[0]
    monthly_file_list = []
    
    for month_string in month_string_dict:
        
        month = month_string_dict[month_string][0]
        year = month_string_dict[month_string][1]
        
        shear_shortname = get_monthly_shortname('SHEAR')
        shear_nc_string = get_monthly_nc_string('SHEAR')
    
        #Look for velocity file in primary directory and download if needed
        download_data.main(field_name='horizontal_vel', initial_month=month, 
                           initial_year=year, final_month=month, 
                           final_year=year, time_ave_type='monthly', 
                           datdir_primary=datdir_primary, time_kwargs=None)
        #Load monthly DataSet
        ds_velocity = load_data_files.main(field_name='horizontal_vel', 
                                           time_ave_type='monthly', 
                                           date_string=month_string, 
                                           datdir_primary=datdir_primary,
                                           datdir_secondary=datdir_secondary) 
        ds_velocity['UVEL'].data = ds_velocity['UVEL'].values
        ds_velocity['VVEL'].data = ds_velocity['VVEL'].values

        xgcm_grid = ecco.get_llc_grid(ds_grid)

        u_mean, v_mean = ds_velocity['UVEL'], ds_velocity['VVEL']
        dx, dy = ds_grid.dxC, ds_grid.dyC
        cell_area = ds_grid.rAz

        #Compute strain
        shear_strain = ((xgcm_grid.diff(v_mean*dy, 'X') + 
                         xgcm_grid.diff(u_mean*dx, 'Y')).squeeze() / cell_area)
        
        shear_strain.name = 'SHEAR'
        
        #Save the monthly data
        if not os.path.exists(join(datdir_secondary, shear_shortname)):
            os.makedirs(join(datdir_secondary, shear_shortname))
        file_path = join(datdir_secondary, shear_shortname,
                        shear_nc_string+month_string+".nc")
        shear_strain.to_netcdf(path=file_path, engine="scipy")
        monthly_file_list.append(file_path)
    
    #If seasonal, compute the seasonally-averaged data if missing
    if time_ave_type == 'seasonal':
        
        shear_shortname = time_ave_type_attrs[1]
        shear_nc_string = time_ave_type_attrs[2]
        
        if not os.path.exists(join(datdir_secondary, shear_shortname)):
            os.makedirs(join(datdir_secondary, shear_shortname))
        
        shear_seasonal = xr.open_mfdataset(monthly_file_list, 
                                               combine='nested', 
                                               concat_dim='month')
        shear_seasonal_avg = (shear_seasonal.sum('month') 
                                  / len(monthly_file_list))
        shear_seasonal_avg.to_netcdf(path=join(datdir_secondary, 
                                        shear_shortname, 
                                        shear_nc_string+date_string+".nc"), 
                                         engine="scipy")

##############################

#VORTICITY

def comp_vorticity(ds_grid, date_string, datdir_primary, datdir_secondary, 
                   time_ave_type, time_kwargs):
    """
    Compute and save vorticity to DataSet in NetCDF format.
    """
    
    time_ave_type_attrs = get_time_ave_type_attrs(date_string, time_ave_type, 
                                                  time_kwargs, 'vorticity')
    month_string_dict = time_ave_type_attrs[0]
    monthly_file_list = []
    
    for month_string in month_string_dict:
        
        month = month_string_dict[month_string][0] 
        year = month_string_dict[month_string][1]
        
        vorticity_shortname = get_monthly_shortname('ZETA')
        vorticity_nc_string = get_monthly_nc_string('ZETA')
    
        #Look for velocity file in primary directory and download if needed
        download_data.main(field_name='horizontal_vel', initial_month=month, 
                           initial_year=year, final_month=month, 
                           final_year=year, time_ave_type='monthly', 
                           datdir_primary=datdir_primary, time_kwargs=None)
        #Load monthly DataSet
        ds_velocity = load_data_files.main(field_name='horizontal_vel', 
                                           time_ave_type='monthly', 
                                           date_string=month_string, 
                                           datdir_primary=datdir_primary,
                                           datdir_secondary=datdir_secondary) 
        ds_velocity['UVEL'].data = ds_velocity['UVEL'].values
        ds_velocity['VVEL'].data = ds_velocity['VVEL'].values

        xgcm_grid = ecco.get_llc_grid(ds_grid)

        u_mean, v_mean = ds_velocity['UVEL'], ds_velocity['VVEL']
        dx, dy = ds_grid.dxC, ds_grid.dyC
        cell_area = ds_grid.rAz

        #Compute vorticity
        vorticity = ((-xgcm_grid.diff(u_mean*dx, 'Y') 
                      + xgcm_grid.diff(v_mean*dy, 'X')).squeeze() / cell_area)
        
        vorticity.name = 'ZETA'
        
        #Save the monthly data
        if not os.path.exists(join(datdir_secondary, vorticity_shortname)):
            os.makedirs(join(datdir_secondary, vorticity_shortname))
        file_path = join(datdir_secondary, vorticity_shortname,
                        vorticity_nc_string+month_string+".nc")
        vorticity.to_netcdf(path=file_path, engine="scipy")
        monthly_file_list.append(file_path)
    
    #If seasonal, compute the seasonally-averaged data if missing
    if time_ave_type == 'seasonal':
        
        vorticity_shortname = time_ave_type_attrs[1]
        vorticity_nc_string = time_ave_type_attrs[2]
        
        if not os.path.exists(join(datdir_secondary, vorticity_shortname)):
            os.makedirs(join(datdir_secondary, vorticity_shortname))
        
        vorticity_seasonal = xr.open_mfdataset(monthly_file_list, 
                                               combine='nested', 
                                               concat_dim='month')
        vorticity_seasonal_avg = (vorticity_seasonal.sum('month') 
                                  / len(monthly_file_list))
        vorticity_seasonal_avg.to_netcdf(path=join(datdir_secondary, 
                                        vorticity_shortname, 
                                        vorticity_nc_string+date_string+".nc"), 
                                         engine="scipy")
        
##############################

def main(**kwargs):
   
    if kwargs:
        datdir_primary = kwargs.get('datdir_primary')
        datdir_secondary = kwargs.get('datdir_secondary')
        date_string = kwargs.get('date_string')
        time_ave_type = kwargs.get('time_ave_type')
        time_kwargs = kwargs.get('time_kwargs')
        field_name = kwargs.get('field_name')
        rho_ref = float(kwargs.get('rho_ref'))
        nu_E = float(kwargs.get('nu_E'))
        
    if time_ave_type == 'monthly':
        field_shortname = get_monthly_shortname(get_variable_str(field_name))
        field_nc_string = get_monthly_nc_string(get_variable_str(field_name))
    elif time_ave_type == 'seasonal':
        field_shortname = get_seasonal_shortname(get_variable_str(field_name))
        field_nc_string = get_seasonal_nc_string(get_variable_str(field_name))
        
    filename = field_nc_string + date_string + '.nc'
    path_to_file = os.path.join(datdir_secondary, field_shortname, filename)

    if not os.path.exists(path_to_file): #Create file if it doesn't exist
    
        ds_grid = load_grid(datdir_primary) #Load the grid DataSet

        #Identify the input field and run the appropriate function to compute it
        if field_name == '2D_div_vel':
            comp_2D_div_vel(ds_grid, date_string, datdir_primary, 
                            datdir_secondary, time_ave_type, time_kwargs)
        elif field_name == 'geostrophic_vel':
            comp_geostrophic_vel(ds_grid, date_string, datdir_primary, 
                                 datdir_secondary, time_ave_type, time_kwargs, 
                                 rho_ref)
        elif field_name == 'Ek_vel':
            comp_Ek_vel(ds_grid, date_string, datdir_primary, datdir_secondary, 
                        time_ave_type, time_kwargs, rho_ref, nu_E)
        elif field_name == 'normal_strain':
            comp_normal_strain(ds_grid, date_string, datdir_primary, 
                               datdir_secondary, time_ave_type, time_kwargs)
        elif field_name == 'shear_strain':
            comp_shear_strain(ds_grid, date_string, datdir_primary, 
                              datdir_secondary, time_ave_type, time_kwargs)
        elif field_name == 'vorticity':
            comp_vorticity(ds_grid, date_string, datdir_primary, 
                           datdir_secondary, time_ave_type, time_kwargs)
        elif field_name == 'T_zflux':
            comp_T_zflux(ds_grid, date_string, datdir_primary, datdir_secondary,
                         time_ave_type, time_kwargs)

##############################

if __name__ == "__main__":
    main()