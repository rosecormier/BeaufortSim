"""
Rosalie Cormier, 2023
"""

##############################

def get_field_variable(field_name):
    
    field_variables = {'density': 'PHIHYDcR', \
                      'pressure': 'PHIHYDcR', \
                      'vertical_vel': 'WVEL', \
                      'horizontal_vel': 'UVELVVEL', \
                      'vorticity': 'ZETA', \
                      'normal_strain': 'NORMAL', \
                      'shear_strain': 'SHEAR', \
                      'div_vel': 'DIVU', \
                      'div_Ek_vel': 'DIVUEk', \
                      'geostrophic_vel': 'UGVG', \
                      'Ek_vel': 'UEkVEk', \
                      'wind_stress': 'EXFtauxEXFtauy'}
    
    return field_variables[field_name]

##############################

def field_is_primary(field_name):
    
    field_variable = get_field_variable(field_name)
    
    primary_variables = ['PHIHYDcR', 'WVEL', 'UVELVVEL', 'EXFtauxEXFtauy']
    secondary_variables = ['ZETA', 'NORMAL', 'SHEAR', 'DIVU', 'DIVUEk', 'UGVG', 'UEkVEk']
    
    if field_variable in primary_variables:
        return True
    elif field_variable in secondary_variables:
        return False

##############################

def get_monthly_shortname(field_variable):
    
    monthly_shortnames = {'PHIHYDcR': 'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4', \
                         'UVELVVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4', \
                         'WVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4', \
                         'UGVG': 'GEOS_VEL_MONTHLY', \
                         'ZETA': 'VORTICITY_MONTHLY', \
                         'NORMAL': 'STRAIN_MONTHLY', \
                         'SHEAR': 'STRAIN_MONTHLY', \
                         'EXFtauxEXFtauy': 'ECCO_L4_STRESS_LLC0090GRID_MONTHLY_V4R4', \
                         'UEkVEk': 'EK_VEL_MONTHLY', \
                         'DIVU': 'DIVU_MONTHLY', \
                         'DIVUEk': 'DIVUEk_MONTHLY'}
    
    return monthly_shortnames[variable]

##############################
    
def get_monthly_nc_string(field_variable):
    
    monthly_nc_strings = {'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_mon_mean_', \
                         'UVELVVEL': 'OCEAN_VELOCITY_mon_mean_', \
                         'WVEL': 'OCEAN_VELOCITY_mon_mean_', \
                         'UGVG': 'OCEAN_GEOS_UVEL_mon_mean_', \
                         'ZETA': 'OCEAN_VORTICITY_mon_mean_', \
                         'NORMAL': 'OCEAN_NORMAL_STRAIN_mon_mean_', \
                         'SHEAR': 'OCEAN_SHEAR_STRAIN_mon_mean_', \
                         'EXFtauxEXFtauy': 'OCEAN_AND_ICE_SURFACE_STRESS_mon_mean_', \
                         'UEkVEk': 'OCEAN_EK_VEL_mon_mean_', \
                         'DIVU': 'OCEAN_DIVU_mon_mean_', \
                         'DIVUEk': 'OCEAN_DIVUEk_mon_mean_'}
    
    return monthly_nc_strings[variable]

##############################

def get_field_string(field_variable):
    
    """
    Strings to be used in filenames.
    """
    
    field_strings = {'PHIHYDcR': 'p_anom', \
                    'UVELVVEL': 'u', \
                    'WVEL': 'w', \
                    'ZETA': 'zeta', \
                    'UGVG': 'u_g', \
                    'EXFtauxEXFtauy': 'tau', \
                    'UEkVEk': 'u_Ek', \
                    'DIVU': 'div_u',
                    'DIVUEk': 'div_u_Ek'}

    field_string = field_strings[field_variable]
        
    return field_string