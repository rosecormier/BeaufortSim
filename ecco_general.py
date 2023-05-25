"""
Rosalie Cormier, 2023
"""

import xarray as xr
import ecco_v4_py as ecco

from os.path import join
from ecco_download import ecco_podaac_download

def load_grid(datdir):
    
    """
    Loads ECCO grid.
    
    datdir = directory where data is stored
    """
    
    grid_params_shortname = "ECCO_L4_GEOMETRY_LLC0090GRID_V4R4"

    #Download ECCO grid parameters (date is arbitrary)
    ecco_podaac_download(ShortName=grid_params_shortname, StartDate="2000-01-01", \
                     EndDate="2000-01-02", download_root_dir=datdir, n_workers=6, \
                     force_redownload=False)

    grid_params_file = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
    grid_params_file_path = join(datdir, grid_params_shortname, grid_params_file)

    #Load grid parameters
    ds_grid = xr.open_dataset(grid_params_file_path)
    
    return ds_grid

def get_starting_i(startmo):

    """
    Get index representing start month.
    
    startmo = starting month
    """
    
    month_dict = {0: "01", 1: "02", 2: "03", 3: "04", 4: "05", 5: "06",
              6: "07", 7: "08", 8: "09", 9: "10", 10: "11", 11: "12"}
    month_key_list, month_val_list = list(month_dict.keys()), list(month_dict.values())
    
    month_index = month_val_list.index(startmo)
    i = month_key_list[month_index]
    
    return i

def get_monthstr(i):
    
    """
    Returns a string corresponding to month i.
    """
    
    month_dict = {1: "01", 2: "02", 3: "03", 4: "04", 5: "05", 6: "06",
              7: "07", 8: "08", 9: "09", 10: "10", 11: "11", 0: "12"}
    
    return month_dict[i % 12]

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

def load_dataset(curr_file):

    """
    Opens ECCO dataset given a file.
    """
    
    dataset = xr.open_mfdataset(curr_file, parallel=True, data_vars='minimal', coords='minimal', compat='override')
    
    return dataset

def get_scalar_in_xy(ecco_ds_grid, k_val, ecco_ds_scalar, scalar_attr, skip_k=False):
    
    """
    Loads scalar field in xy-grid.
    
    ecco_ds_grid = ECCO grid
    k_val = depth index of interest
    ecco_ds_scalar = DataSet containing field
    scalar_attr = name of field
    skip_k = whether to skip isolating a value of k (True if value has already been isolated)
    """
    
    ds_grid = ecco_ds_grid.copy()

    if skip_k:
        ds_grid[scalar_attr] = ecco_ds_scalar[scalar_attr]
    
    else:
        ecco_ds_scalar_k = ecco_ds_scalar.isel(k=k_val)
        ds_grid[scalar_attr] = ecco_ds_scalar_k[scalar_attr]
        
    ds_grid = ds_grid.load()
    
    return ds_grid
    
def get_vector_in_xy(ecco_ds_grid, k_val, ecco_ds_vector, xvec_attr, yvec_attr, skip_k=False):
    
    """
    Loads vector field in xy-grid.
    
    ecco_ds_grid = ECCO grid
    k_val = depth index of interest
    ecco_ds_vector = DataSet containing vector field
    xvec_attr = name of x-comp of vector field
    yvec_attr = name of y-comp of vector field
    skip_k = whether to skip isolating a value of k (True if already isolated)
    """

    ds_grid = ecco_ds_grid.copy()
    ds_grid = ds_grid.load()
    XGCM_grid = ecco.get_llc_grid(ds_grid)
    
    if skip_k:
        velc = XGCM_grid.interp_2d_vector({'X': ecco_ds_vector[xvec_attr], \
                                           'Y': ecco_ds_vector[yvec_attr]}, \
                                           boundary='fill')
                                           
    else:
        velc = XGCM_grid.interp_2d_vector({'X': (ecco_ds_vector[xvec_attr]).isel(k=k_val), \
                                                 'Y': (ecco_ds_vector[yvec_attr]).isel(k=k_val)}, \
                                                 boundary='fill')

    return velc

def comp_temp_mean_scalar(k_val, ecco_ds_scalars, scalar_attr):
    
    """
    Computes temporal mean of a scalar field.

    k_val = depth value of interest
    ecco_ds_scalars = scalar Datasets
    scalar_attr = string corresponding to scalar attribute of interest
    """ 
    
    mean_field = ((ecco_ds_scalars[0]).copy()) * 0 
    concat_field = ((ecco_ds_scalars[0]).copy()).isel(k=k_val) * 0
    
    for dataset in ecco_ds_scalars:
        curr_field = dataset.isel(k=k_val)
        concat_field = xr.concat((concat_field, curr_field), dim='time')

    mean_scalar_field = (concat_field[scalar_attr]).sum(dim=['time']) / len(ecco_ds_scalars)
    mean_scalar_field = mean_scalar_field.where(mean_scalar_field != 0)
    mean_field[scalar_attr] = mean_scalar_field
    (mean_field[scalar_attr]).data = (mean_field[scalar_attr]).values

    skip_k = True
    
    return mean_field, skip_k

def comp_temp_mean_vector(k_val, ecco_ds_vectors, xvec_attr, yvec_attr):
    
    """
    Computes temporal mean of a vector field.
    
    k_val = depth value of interest
    ecco_ds_vectors = vector Datasets
    xvec_attr = string corresponding to x-comp of attribute of interest
    yvec_attr = string corresponding to y-comp of attribute of interest
    """
    
    mean_fields = ((ecco_ds_vectors[0]).copy()) * 0
    concat_x_field = ((ecco_ds_vectors[0]).copy()).isel(k=k_val) * 0
    concat_y_field = ((ecco_ds_vectors[0]).copy()).isel(k=k_val) * 0
    
    for dataset in ecco_ds_vectors:
        
        curr_x_field = dataset.isel(k=k_val)
        curr_y_field = dataset.isel(k=k_val)
        
        concat_x_field = xr.concat((concat_x_field, curr_x_field), dim='time')
        concat_y_field = xr.concat((concat_y_field, curr_y_field), dim='time')
        
    mean_x_field = (concat_x_field[xvec_attr]).sum(dim=['time']) / len(ecco_ds_vectors)
    mean_x_field = mean_x_field.where(mean_x_field != 0)
    mean_y_field = (concat_y_field[yvec_attr]).sum(dim=['time']) / len(ecco_ds_vectors)
    mean_y_field = mean_y_field.where(mean_y_field != 0)
    
    mean_fields[xvec_attr] = mean_x_field
    (mean_fields[xvec_attr]).data = (mean_fields[xvec_attr]).values
    mean_fields[yvec_attr] = mean_y_field
    (mean_fields[yvec_attr]).data = (mean_fields[yvec_attr]).values
    
    skip_k = True
    
    return mean_fields, skip_k

def comp_residuals(attributes, datasets, mean_dataset):
    
    """
    Computes residuals relative to a mean.
    
    attributes = variables to compute residuals of
    datasets = datasets to compute residuals for
    mean_dataset = mean of all attributes
    """
        
    residual_list = []
        
    for dataset in datasets:
            
        residual = dataset.copy() * 0 
            
        concat_attr = xr.concat((dataset, -1*mean_dataset), dim='time')
        
        for attribute in attributes:
            residual_field = (concat_attr[attribute]).sum(dim=['time'])
            residual_field = residual_field.where(residual_field != 0)
            residual[attribute] = residual_field
            (residual[attribute]).data = (residual[attribute]).values
            
        residual_list.append(residual)
    
    return residual_list