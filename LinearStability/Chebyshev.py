"""
Modification of Storer's code cheb.py.

I reference "Spectral Methods in MATLAB" by L. Trefethen.
"""

import numpy as np

from math import pi

def Chebyshev(N, xTransform = None):
    """
    Computes the Chebyshev differentiation matrix on N+1 points 
     (i.e., N intervals).
    Returns:
      D = (N+1) x (N+1) Chebyshev differentiation matrix
      x = Chebyshev grid
    """

    if N == 0:
        D, x = 0, 1

    else:

        #Create a vector of N+1 Chebyshev-spaced components from 1 to -1
        x = np.cos(pi * np.arange(0, (N+1)) / N)

        if xTransform is not None:
            x = xTransform(x)

        #Define the Chebyshev coeffs c_{ij}
        c = np.hstack([2, np.ones(N-1), 2]) * (-1)**np.arange(0, (N+1))
        
        X  = np.tile(x, (N+1, 1))
        dX = X.T - X #Element-wise difference between X^T and X

        #Initialize D with off-diag entries computed by eq. 6.5 (Trefethen)
        D = ((c[:, np.newaxis] * (1.0 / c)[np.newaxis, :])
             / (dX + (np.identity(N+1))))

        #Update diagonal entries (currently = 0) using identity 6.6 (Trefethen)
        D = D - np.diag(D.sum(axis = 1))

    return D, x
