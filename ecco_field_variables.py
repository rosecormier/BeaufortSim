"""
Rosalie Cormier, 2023
"""

def get_field_vars(attribute):
    
    monthly_shortnames = {'PHIHYDcR': 'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4', \
                         'UVELVVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4', \
                         'UGVG': 'GEOS_VEL_MONTHLY', \
                         'ZETA': 'VORTICITY_MONTHLY', \
                         'NORMAL': 'STRAIN_MONTHLY', \
                         'SHEAR': 'STRAIN_MONTHLY'}
    
    monthly_nc_strings = {'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_mon_mean_', \
                         'UVELVVEL': 'OCEAN_VELOCITY_mon_mean_', \
                         'UGVG': 'OCEAN_GEOS_UVEL_mon_mean_', \
                         'ZETA': 'OCEAN_VORTICITY_mon_mean', \
                         'NORMAL': 'OCEAN_NORMAL_STRAIN_mon_mean', \
                         'SHEAR': 'OCEAN_SHEAR_STRAIN_mon_mean'}
    
    return monthly_shortnames[attribute], monthly_nc_strings[attribute]

def get_variable_str(attribute, geostrophic=False):
    
    variables = {'PHIHYDcR': 'p_anom', \
                'UVELVVEL': 'u', \
                'ZETA': 'zeta', \
                'UGVG': 'u_g'}

    variable_string = variables[attribute]
        
    return variable_string