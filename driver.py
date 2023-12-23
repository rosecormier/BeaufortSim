"""
This script uses the data provided in the input file (input.py) to do all of the following:
    -Download data corresponding to specified fields (skips if done);
    -Computes and saves secondary fields (skips if unnecessary or done);
    -Plots data for all fields specified and saves plots; and
    -Removes primary data files.

R. Cormier, F. Poulin, 2023
"""

import numpy as np 
import matplotlib.pyplot as plt 

### Read input.txt

# read_input.py

### read primary data files to be plotted

scalar_fields = list(plot_fields.keys())
vector_fields = list(plot_fields.values())

# To-do: what if either is empty?

for scalar_field in scalar_fields:
    # check if primary vs secondary
    # may have to download
    # may have to compute (secondary)

for vector_field in vector_fields:
    # check if primary vs secondary
    # may have to download
    # may have to compute (secondary)

# try and read in secondary fields to be plotted
#   if yes, good!
#   else compute secondary fields and save them
#       much work could be done here
#   return


### Visualize data

for scalar_field in scalar_fields:
    vector_field = plot_fields[scalar_field]

    # plot and save


### Remove data files

