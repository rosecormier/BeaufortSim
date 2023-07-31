"""
Rosalie Cormier, 2023
"""

def get_field_vars(attribute):
    
    monthly_shortnames = {'PHIHYDcR': 'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4', \ #Density, strat, pressure
                         'UVELVVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4', \ #u, v
                         'WVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4', \ #w
                         'UGVG': 'GEOS_VEL_MONTHLY', \ #Geostrophic velocity (doesn't come from ECCO)
                         'ZETA': 'VORTICITY_MONTHLY', \ #Vorticity (doesn't come from ECCO)
                         'NORMAL': 'STRAIN_MONTHLY', \ #Normal strain (doesn't come from ECCO)
                         'SHEAR': 'STRAIN_MONTHLY', \ #Shear strain (doesn't come from ECCO)
                         'EXFtauxEXFtauy': 'ECCO_L4_STRESS_LLC0090GRID_MONTHLY_V4R4', \ #Tau_x and Tau_y
                         'UEkVEk': 'EK_VEL_MONTHLY'} #Ekman velocity (doesn't come from ECCO)
    
    monthly_nc_strings = {'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_mon_mean_', \
                         'UVELVVEL': 'OCEAN_VELOCITY_mon_mean_', \
                         'WVEL': 'OCEAN_VELOCITY_mon_mean_', \
                         'UGVG': 'OCEAN_GEOS_UVEL_mon_mean_', \
                         'ZETA': 'OCEAN_VORTICITY_mon_mean', \
                         'NORMAL': 'OCEAN_NORMAL_STRAIN_mon_mean', \
                         'SHEAR': 'OCEAN_SHEAR_STRAIN_mon_mean', \
                         'EXFtauxEXFtauy': 'OCEAN_AND_ICE_SURFACE_STRESS_mon_mean_', \
                         'UEkVEk': 'OCEAN_EK_VEL_mon_mean'}
    
    return monthly_shortnames[attribute], monthly_nc_strings[attribute]

def get_variable_str(attribute, geostrophic=False):
    
    variables = {'PHIHYDcR': 'p_anom', \
                'UVELVVEL': 'u', \
                'WVEL': 'w', \
                'ZETA': 'zeta', \
                'UGVG': 'u_g', \
                'EXFtauxEXFtauy': 'tau', \
                'UEkVEk': 'u_Ek'}

    variable_string = variables[attribute]
        
    return variable_string