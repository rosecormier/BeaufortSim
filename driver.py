import numpy as np 
import matplotlib.pyplot as plt 

### Read input.txt

# read_input.py

### read primary data files to be plotted

scalar_fields = list(plot_fields.keys())
vector_fields = list(plot_fields.values())

# To-do: what if either is empyt?

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
#   else compute seconary fields and save them
#       much work could be done here
#   return


### Visualize data

for scalar_field in scalar_fields:
    vector_field = plot_fields[scalar_field]

    # plot and save


### Remove data files

