"""
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

Some of the notation follows "Spectral Methods in MATLAB" by L.N. Trefethen.
All variables are assumed to be dimensionless unless otherwise specified.
"""

import argparse
import numpy as np
import scipy.linalg as spalg
#Note: I think we should aim to switch from matrix objects (e.g. built by sp.spdiags) to array objects, as the matrix functionality is now deprecated
import sys
import time
import timeit

from math import pi

from BuildLaplacian import BuildLaplacian
from BuildBkgdOperators import BuildBkgdOperators, GridInterior
from Chebyshev import Chebyshev
from SaveToNetCDF import SaveToNetCDF

#Parse command-line inputs

parser = argparse.ArgumentParser()
parser.add_argument('--NrEig',
                    help = 'Number (must be ODD) of r-gridpoints for eig computations',
                    type = int, default = 201)
parser.add_argument('--NzEig',
                    help = 'Number (must be EVEN) of z-gridpoints for eig computations',
                    type = int, default = 100)
parser.add_argument('-Lr', 
                    help = 'DIMENSIONLESS radius of physical domain',
                    type = float, default = 10.0)
parser.add_argument('-Lz',
                    help = 'DIMENSIONLESS depth (> 0) of physical domain',
                    type = float, default = 3.3)
parser.add_argument('-Ro', '--Rossby',
                    help = 'Rossby number of background flow', 
                    type = float, default = 4e-3)
parser.add_argument('-Bu', '--Burger',
                    help = 'Burger number of background flow',
                    type = float, default = 2.5e-3)
parser.add_argument('-f0', '--Coriolis',
                    help = 'DIMENSIONAL Coriolis frequency f0 (Hz)',
                    type = float, default = 1.4e-4)
parser.add_argument('--strat_shape',
                    help = 'Shape of ambient squared buoyancy frequency profile',
                    type = str, default = "constant")
parser.add_argument('-p', '--PrintOutputs',
                    help = 'Flag to turn on display for each computation',
                    action = 'store_true')
parser.add_argument('-kp', '--k_phi', 
                    help = 'Azimuthal wavenumbers; enter as -kp start stop step',
                    type = float, default = [1, 10, 1], nargs = 3)
parser.add_argument('--modes', 
                    help = 'Number of modes of instability to be considered',
                    type = int, default = 1)
parser.add_argument('-Np',
                    help = 'Number of points for discretization of phi',
                    type = int, default = 50)
parser.add_argument('--z_idx',
                    help = 'Constant z-index to plot 2D slices at',
                    type = int, default = 0)
#parser.add_argument('--r_idx',
#                    help = 'Constant r-index to plot 2D slices at',
#                    type = int, default = args.Nr // 2 - args.Nr // (2 * args.Lr))
args = parser.parse_args()
                    
class Parameters:
    
    Ro   = args.Rossby
    Bu   = args.Burger
    f0   = args.Coriolis
    bkgd = "BG"

    stratification_kw = args.strat_shape

    Lr     = args.Lr         #Max. r in physical space; half of computational domain
    Nr     = args.NrEig      #Number of computational gridpoints in r
    halfNr = args.NrEig // 2 #Number of physical r-gridpoints

    Lz = args.Lz    #Max. depth (i.e., -min(z)) in physical domain
    Nz = args.NzEig #Number of computational gridpoints in z

    dφ = 2 * pi / args.Np #φ-increment; for visualization only
    φ  = dφ * np.arange(0, args.Np)

    z_idx = args.z_idx
    r_idx = int(halfNr - halfNr // Lr)

    kφs    = np.arange(args.k_phi[0], args.k_phi[1], args.k_phi[2])        
    nmodes = args.modes

    printout = args.PrintOutputs

    def display(self):
        print(f"Ro = {self.Ro}")
        print(f"Bu = {self.Bu}")
        print(f"f0 = {self.f0}")
        print(f"Lr = {self.Lr}")
        print(f"Nr = {self.Nr}")
        print(f"halfNr = {self.halfNr}")
        print(f"Np = {self.Np}")
        print(f"Nz = {self.Nz}")
        print(f"Lz = {self.Lz}")
        print(f"kps = {self.kps}")
        print(f"nmodes = {self.nmodes}")

def zTransform():
    return lambda z : (z - 1) / 2

class Geometry:

    def __init__(self, method, params):
        
        self.method = method
        
        if method == "Chebyshev":
           
            #Compute r- and z-diff. matrices and Chebyshev-spaced grids
            Dr, r = Chebyshev(params.Nr)
            Dz, z = Chebyshev(params.Nz, xTransform = zTransform())
            
            #Scale gridpoints and variables of differentiation to fit domain
            self.r, self.Dr = r * params.Lr, Dr / params.Lr
            self.z, self.Dz = z * params.Lz, Dz / (params.Lz)
 
            self.Dr2 = np.matmul(self.Dr, self.Dr) #Second-order diff. matrix
 
def QG_Vortex_Stability():
 
    #Initialize parameters and set up geometry for Chebyshev solver
    paramsCh   = Parameters()
    geomCh     = Geometry("Chebyshev", paramsCh)
    geomCh.Lap = BuildLaplacian(paramsCh, geomCh, discretizeVertical = True)

    #Perform 2D discretization of grid
    GridInterior(paramsCh, geomCh, discretizeVertical = True)

    #Discretize background-state-flow operators
    geomCh.Ψ_op, geomCh.Q_op = BuildBkgdOperators(paramsCh, geomCh, 
                                                  discretizeVertical = True)

    #Information about wavenumbers and modes
    kφs, nmodes = paramsCh.kφs, paramsCh.nmodes

    #Initialize arrays to store results of eigen-computations
    growthNondimCh      = np.zeros([kφs.shape[0], nmodes])
    growthDimensionalCh = np.copy(growthNondimCh)
    propNondimCh        = np.zeros([kφs.shape[0], nmodes])
    propDimensionalCh   = np.copy(propNondimCh)
    eigModesCh          = np.zeros([kφs.shape[0], (paramsCh.halfNr * paramsCh.Nz), 
                                    nmodes], dtype = complex)

    ########################################
    # SOLVE GENERALIZED EIGENVALUE PROBLEM #
    ########################################
    
    for kφ_idx in range(0, kφs.shape[0]):

        kφ  = kφs[kφ_idx]
        kφ2 = kφ**2
    
        ##############################
        # BUILD MATRICES 'A' and 'B' #
        ##############################

        B_Ch = (geomCh.Lap - (kφ2 * geomCh.r2Recip)
                + ((1 / paramsCh.Bu) * (geomCh.Dz 
                    * (geomCh.N2Recip * geomCh.Dz)))
               )
        A_Ch = np.matmul(geomCh.Ψ_op, B_Ch) - geomCh.Q_op

        ############################
        # FIND EIGENSPACE DIRECTLY #
        ############################
            
        t0Ch = timeit.timeit()
        
        #Compute eigvals c and eigvecs psi with direct solver
        eigValsCh, eigVecsCh = spalg.eig(A_Ch, B_Ch)

        solveTimeCh = timeit.timeit() - t0Ch #Time for direct solver

        #Indexing that sorts eigvals by ASCENDING Im(c)
        indSortCh = np.argsort(eigValsCh.imag)
            
        eigValsCh = eigValsCh[indSortCh]    #Sort eigvals
        eigVecsCh = eigVecsCh[:, indSortCh] #Sort eigvecs in the same order
        ωsCh      = eigValsCh * kφ          #Corresponding omegas for this kφ
        
        ωsDimensionalCh = ωsCh * paramsCh.Ro * paramsCh.f0
        #print(eigValsCh[0:nmodes])

        #Store results of direct solve
        growthNondimCh[kφ_idx, :]           = -ωsCh[0:nmodes].imag
        growthDimensionalCh[kφ_idx, :]      = -ωsDimensionalCh[0:nmodes].imag
        propNondimCh[kφ_idx, :]             = ωsCh[0:nmodes].real
        propDimensionalCh[kφ_idx, :]        = -ωsDimensionalCh[0:nmodes].real
        eigModesCh[kφ_idx, :, :]            = eigVecsCh[:, 0:nmodes]

    SaveToNetCDF(paramsCh, np.arange(1, nmodes + 1, 1), growthDimensionalCh,
                 propDimensionalCh, eigModesCh, geomCh.r, paramsCh.φ, geomCh.z)

if __name__ == '__main__': #For testing
   QG_Vortex_Stability()
