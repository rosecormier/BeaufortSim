"""
Modification of Storer's code cheb.py.

I reference "Spectral Methods in MATLAB" by L. Trefethen.
"""

import numpy as np

from math import pi

def Chebyshev(N):
    """
    Computes the Chebyshev differentiation matrix on N+1 points 
     (i.e., N intervals).
    Returns:
      D = (N+1) x (N+1) Chebyshev differentiation matrix
      x = Chebyshev grid
    """

    if N == 0:
        D, x = 0, 1
  
    #else:
        
        #Create a vector of N+1 Chebyshev-spaced components from 1 to -1
        x = (np.cos(pi * np.array(range(N+1)) / N)).reshape([N+1, 1])

        X = np.tile(x, (1, N+1))
        #Note: Let x = [[x0], [x1], ..., [xN]].
        #X is a square (N+1 x N+1) matrix with entries:
        # [[x0, x0, ..., x0],
        #  [x1, x1, ..., x1],
        #  ...
        #  [xN, xN, ..., xN]]

        dX = X - X.transpose() 
        #Element-wise difference between X and X^T

        #Define the Chebyshev coeffs c_{ij}
        c = (np.ravel(np.vstack([2, np.ones([N-1, 1]), 2]))
             * (-1.0)**np.ravel(np.array(range(-N,N+1,2))))
        #Note: * above is an element-wise product, not matrix product
        c = c.reshape(c.shape[0], 1)
        #The resulting c is an (N+1)-component vector with entries:
        # [[2.0], [-1.0], ..., [(-1.0)^(N-1)], [2(-1.0)^N]]

        #Initialize D with off-diag entries computed by eq. 6.5 (Trefethen)
        D = (c * (1/c).transpose()) / (dX + (np.eye(N+1)))
        
        #Update diagonal entries (currently = 0) using identity 6.6 (Trefethen)
        D = D - np.diag(np.sum(D, 1)) 
   
    else:

        x = np.cos(pi*np.arange(0,N+1)/N)
        c = np.hstack([2, np.ones(N-1), 2])*(-1)**np.arange(0,N+1)
        X = np.tile(x,(N+1,1))
        dX = X.T - X
        D = (c[:,np.newaxis]*(1.0/c)[np.newaxis,:])/(dX+(np.identity(N+1)))       # off-diagonal entries
        D = D - np.diag(D.sum(axis=1))

    return D, x
