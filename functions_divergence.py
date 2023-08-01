"""
Functions to compute velocity divergence.

Rosalie Cormier, 2023
"""

##############################

def comp_2d_divergence(grid_llc, u, v, dx, dy, cell_area):
    return (grid_llc.diff(u * dy, 'X') + grid_llc.diff(v * dx, 'Y')).squeeze() / cell_area