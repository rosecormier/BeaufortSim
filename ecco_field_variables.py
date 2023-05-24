def get_month_end(monthstr, yearstr):
    
    """
    Return string representing last day of specified month.
    """
    
    month_end_dict = {"01": "31", "03": "31", "04": "30",
                 "05": "31", "06": "30", "07": "31", "08": "31",
                 "09": "30", "10": "31", "11": "30", "12": "31"}
    
    if monthstr != "02":
        endmonth = month_end_dict[monthstr]
        
    elif monthstr == "02":
        if int(yearstr) % 4 == 0:
            endmonth = "29"
        else: 
            endmonth = "28"
    
    return endmonth

def get_scalar_field_vars(scalar_attr):
    
    monthly_shortnames = {'PHIHYDcR': 'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4'}
    monthly_nc_strings = {'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_mon_mean_'}
    
    scalar_monthly_shortname = monthly_shortnames[scalar_attr]
    scalar_monthly_nc_string = monthly_nc_strings[scalar_attr]
    
    return scalar_monthly_shortname, scalar_monthly_nc_string
    
def get_vector_field_vars(xvec_attr, yvec_attr):
    
    if xvec_attr == 'UVEL' and yvec_attr == 'VVEL':
        vec_field_name = 'velocity'
        
    else:
        print('Invalid vector component combination.')
    
    monthly_shortnames = {'velocity': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4'}
    monthly_nc_strings = {'velocity': 'OCEAN_VELOCITY_mon_mean_'}
    
    vector_monthly_shortname = monthly_shortnames[vec_field_name]
    vector_monthly_nc_string = monthly_nc_strings[vec_field_name]
    
    return vector_monthly_shortname, vector_monthly_nc_string