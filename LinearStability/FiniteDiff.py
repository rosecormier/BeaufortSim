import numpy as np
import numpy.linalg as nlg
import scipy
import scipy.sparse as sp
from scipy.special import factorial
import scipy.linalg as spalg

def FiniteDiff(x, n, spb, uniform):
    """
    Create a finite difference matrix of arbitrary order for an arbitrary grid.
    x       = discretized domain (grid)
    n       = order of stencil?
    spb     = whether matrix is sparse? 
    uniform = whether points in x are uniformly spaced

    Note that x can be provided in the form [a, b, c].
    This is interpreted as x = linspace(a, b, c).
    The advantage of this shorthand is that we don't actually need to
     generate x, so we can save on memory when using large grids.
    """

    n = int(n)
    
    if len(x) == 3: #If using shorthand for linspace
        Nx = x[2] #Number of points in discrete domain
    else:
        Nx = len(x)

    if spb: #If differentiation matrix will be sparse
        Dx = sp.lil_matrix((Nx, Nx)) #Construct an empty Nx x Nx matrix
    else:
        Dx = np.zeros([Nx, Nx]) #Initialize an Nx x Nx matrix with zeros

    if uniform: #If grid is uniformly spaced

        if len(x) == 3:
            dx = (x[1] - x[0]) / float((x[2] - 1)) #Grid spacing
        else:
            dx = x[1] - x[0]

        #Define entries in Dx for gridpoints that are closer to 
        # the boundary than half the stencil width
        for i in range(0, int(np.ceil(n/2))):
            
            A = np.zeros([n+1, n+1])
            
            for j in range(0, n+1):
                A[:, j] = (np.power((j-i) * dx * np.ones([1, n+1]), 
                                    range(0, n+1)) 
                           / factorial(range(0, n+1)))
            
            b    = np.zeros(n+1)
            b[1] = 1
            
            coeff = nlg.solve(A, b)
            coeff = coeff.conj().transpose()
            
            Dx[i, 0:(n+1)] = coeff

        #Define entries in Dx for gridpoints that are closer to
        # the boundary than half the stencil width
        for i in range(Nx - int(np.ceil(n/2)), Nx):

            A = np.zeros([n+1, n+1])
            
            for j in range(Nx-n-1,Nx):
                A[:, (j-Nx+n+1)] = (np.power((j-i) * dx * np.ones([1, n+1]),
                                             range(0, n+1))
                                    / factorial(range(0, n+1)))

            b    = np.zeros(n+1)
            b[1] = 1

            coeff = nlg.solve(A, b)
            coeff = coeff.conj().transpose()

            Dx[i, (Nx-n-1):Nx] = coeff

        #Now define entries in Dx for interior gridpoints

        A = np.zeros([n+1, n+1])
        
        if n % 2 == 0: #If finite-differencing order is even

            for j in range(-n//2 - 1, n//2):
                A[:, (j + (n//2) + 1)] = (np.power((j+1) * dx 
                                                   * np.ones([1, n+1]), 
                                                   range(0, n+1)) 
                                                 / factorial(range(0, n+1)))
      
            b    = np.zeros(n+1)
            b[1] = 1
            
            coeff = nlg.solve(A, b)
            coeff = coeff.conj().transpose()
            coeff = np.tile(coeff, [Nx, 1])

            Dx[(n//2):(Nx - n//2), :] = sp.spdiags(coeff.T, 
                                                   range(0, n+1), 
                                                   Nx-n, Nx).todense()

        elif n % 2 == 1: #If finite-differencing order is odd

            for j in range((int(-n//2) - 1), int(np.ceil(n/2))):
                A[:, (j + int(n//2) + 1)] = (((j*dx)**range(0, n+1))
                                             / factorial(range(0, n+1)))

            b    = np.zeros(n+1)
            b[1] = 1

            coeff = nlg.solve(A, b)
            coeff = coeff.conj().transpose()
            coeff = np.tile(coeff, [Nx, 1])

            Dx[int(np.ceil(n/2)):(Nx - int(np.ceil(n/2))), :] = sp.spdiags(
                             coeff.T, range(0, n+1), (Nx - (n+1)), Nx).todense()
    
    else: #If grid is not uniformly spaced

        for i in range(0, Nx):

            #Define entries in Dx for gridpoints that are closer to
            # the boundary than half the stencil width
            if i <= np.ceil(n/2):
                
                A = np.zeros([n+1, n+1])
                
                for j in range(0, n+1):
                    dx     = x[j] - x[i]
                    A[:,j] = (dx**range(0, n+1)) / factorial(range(0, n+1))
                
                b    = np.zeros(n+1)
                b[1] = 1
                
                coeff = nlg.solve(A, b)
                coeff = coeff.conj().transpose()

                Dx[i, 0:(n+1)] = coeff
            
            #Define entries in Dx for gridpoints that are closer to
            # the boundary than half the stencil width
            elif i > (Nx - ceil(n/2)):
                
                A = np.zeros([n+1, n+1])

                for j in range((Nx - (n+1)), Nx):
                    dx                     = x[j] - x[i]
                    A[:, (j - (Nx - (n+1)))] = (dx**range(0, n+1) 
                                              / factorial(range(0, n+1)))

                b    = np.zeros(n+1)
                b[1] = 1
                
                coeff = nlg.solve(A, b)
                coeff = coeff.conj().transpose()

                Dx[i, (Nx - (n+1)):Nx] = coeff

            #Now define entries in Dx for interior gridpoints
            elif i <= (Nx - np.ceil(n/2)):
                
                #If n is even, use a centred scheme
                if n % 2 == 0:
                    
                    A = np.zeros([n+1, n+1])
                    
                    for j in range((-(n/2) - 1), n/2):
                        dx                  = x[i+j] - x[i]
                        A[:, (j + n/2 + 1)] = (dx**range(0, n+1)
                                               / factorial(range(0, n+1)))

                    b    = np.zeros(n+1)
                    b[1] = 1

                    coeff = nlg.solve(A, b)
                    coeff = coeff.conj().transpose()

                    Dx[i, (i - n/2):(i + n/2)] = coeff

                #If n is odd, then bias to the side with the closest gridpoint
                elif n % 2 == 1:

                    if (abs(x[(i + int(np.ceil(n/2)))] - x[i]) 
                            <= abs(x[(i - int(ceil(n/2)))] - x[i])):
                        
                        A = np.zeros(n+1, n+1)

                        for j in range((-int(n//2) - 1), int(np.ceil(n/2))):
                            dx                                 = x[i+j] - x[i]
                            A[:, (j + int(n//2) + 1)] = (dx**range(0, n+1)
                                                         / factorial(
                                                                 range(0, n+1)))
                      
                        b    = np.zeros(n+1)
                        b[1] = 1

                        coeff = nlg.solve(A, b)
                        coeff = coeff.conj().transpose()

                        Dx[i, (i - int(n//2)):(i + int(np.ceil(n/2)))] = coeff

                    else:

                        A = np.zeros([n+1, n+1])
                        
                        for j in range((-int(np.ceil(n/2)) - 1), int(n//2)):
                            dx                                = x[i+j] - x[i]
                            A[:, (j + int(np.ceil(n/2)) + 1)] = (dx**range(0, n)
                                                                 / factorial(
                                                                  range(0, n)))

                        b    = np.zeros(n+1)
                        b[1] = 1

                        coeff = nlg.solve(A, b)
                        coeff = coeff.conj().tranpose()

                        Dx[i, (i - int(np.ceil(n/2))):(i + int(n//2))] = coeff

    if spb: #If making a sparse matrix
        Dx = Dx.tocsr() #Convert to Compressed Sparse Row format
            
    return Dx

if __name__ == '__main__': #For testing
    
    Nx = 10
    Lx = 2
    
    diff_ordx = 2
    
    x = np.linspace(0, Lx, Nx+1)
    
    print('Sparse Test')
    Dx = FiniteDiff(x, diff_ordx, True, True)
    print(Dx)

    print('Non-Sparse Test')
    Dx = FiniteDiff(x, diff_ordx, False, True)
    print(Dx)
