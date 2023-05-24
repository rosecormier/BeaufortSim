"""
Rosalie Cormier, 2023
"""

def get_scalar_field_vars(scalar_attr):
    
    monthly_shortnames = {'PHIHYDcR': 'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4'}
    monthly_nc_strings = {'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_mon_mean_'}
    variables = {'PHIHYDcR': 'p_anom'}
    
    scalar_monthly_shortname = monthly_shortnames[scalar_attr]
    scalar_monthly_nc_string = monthly_nc_strings[scalar_attr]
    scalar_variable = variables[scalar_attr]
    
    return scalar_monthly_shortname, scalar_monthly_nc_string, scalar_variable
    
def get_vector_field_vars(xvec_attr, yvec_attr):
    
    if xvec_attr == 'UVEL' and yvec_attr == 'VVEL':
        vec_field_name = 'velocity'
        
    else:
        print('Invalid vector component combination.')
    
    monthly_shortnames = {'velocity': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4'}
    monthly_nc_strings = {'velocity': 'OCEAN_VELOCITY_mon_mean_'}
    variables = {'velocity': 'u'}
    
    vector_monthly_shortname = monthly_shortnames[vec_field_name]
    vector_monthly_nc_string = monthly_nc_strings[vec_field_name]
    vector_variable = variables[vec_field_name]
    
    return vector_monthly_shortname, vector_monthly_nc_string, vector_variable