import numpy as np

from math import pi

def cheb(N):
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
        
        #Create a vector of N+1 Chebyshev-spaced components ranging [0, N]
        x = (np.cos(pi * np.array(range(N+1)) / N)).reshape([N+1, 1])
        
        X = np.tile(x, (1, N+1))
        #Note: Let x = [[x0], [x1], ..., [xN]].
        #X is a square (N+1 x N+1) matrix with entries:
        # [[x0, x0, ..., x0],
        #  [x1, x1, ..., x1],
        #  ...
        #  [xN, xN, ..., xN]]

        dX = X - X.conj().transpose()
        #This is the element-wise difference between X and (X*)^T.
        #All entries of x are real, so X = X* and thus dX is antisymmetric.

        #Define the Chebyshev coeffs c_{ij}
        c = (np.ravel(np.vstack([2, np.ones([N-1, 1]), 2]))
             * (-1.0)**np.ravel(np.array(range(N+1))))
        #Note: * above is an element-wise product, not matrix product
        c = c.reshape(c.shape[0], 1)
        #The resulting c is an (N+1)-component vector with entries:
        # [[2.0], [-1.0], ..., [(-1.0)^(N-1)], [2(-1.0)^N]]

        print(c)
        print(1/c.conj())
        print(c * (1/c).conj().transpose())
        D = (c * (1/c).conj().transpose()) / (dX + (np.eye(N+1))) #Off-diags
        D = D - np.diag(np.sum(D, 1)) #Diagonal entries
   
    return D, x
