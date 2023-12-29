"""
Contains functions (associated with computed data) to:
    -Create a file containing computed data (if it doesn't already exist);
    -Load a computed DataSet.

Rosalie Cormier, 2023
"""

import os
import xarray as xr

##############################

def create_comp_data_file():
    
    """
    Checks that a computed file exists, and creates it if it doesn't.
    Note - only does monthly avgs at the moment - will fix in future.
    """
    
    #define variable 'filename' here
    
    if not os.path.exists(filename): #Execute only if the file doesn't already exist
        #create the file - to be added
        
##############################

def load_comp_data_file():
    
    """
    Loads computed DataSet.
    """
    
    #define variable 'filename' here
    
    try:
        comp_ds = xr.open_mfdataset(filename, engine="scipy") #Load DataSet
        return comp_ds
    
    except:
        print("DataSet does not exist.")