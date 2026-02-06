"""
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

Some of the notation follows "Spectral Methods in MATLAB" by L.N. Trefethen.
All variables are assumed to be dimensionless.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as nlg
import scipy.linalg as spalg

import sys
import time
import timeit

from math import pi
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from BuildLaplacian import BuildLaplacian
from BuildBkgdOperators import BuildBkgdOperators, GridInterior
from Chebyshev import Chebyshev
from Streamfunctions import Streamfunction, EigenvelocityFrom1DEigvec

parser = argparse.ArgumentParser()
parser.add_argument('--Neig', 
                    help = 'Number (must be ODD) of grid points for eig computations',
                    type = int, default = 801)
parser.add_argument('-Lr', 
                    help = 'DIMENSIONLESS radius of the physical domain',
                    type = float, default = 8.0)
parser.add_argument('-Ro', '--Rossby',
                    help = 'Rossby number of background flow', 
                    type = float, default = 4e-3)
parser.add_argument('-Bu', '--Burger',
                    help = 'Burger number of background flow',
                    type = float, default = 2.5e-3)
parser.add_argument('-f0', '--Coriolis',
                    help = 'Coriolis frequency f0 (Hz)',
                    type = float, default = 1.4e-4)
parser.add_argument('--bkgd',
                    help = 'Background flow to use ("GM" or "BG")',
                    type = str, default = "BG")
parser.add_argument('-p', '--PrintOutputs',
                    help = 'Flag to turn on display for each computation',
                    action = 'store_true')
parser.add_argument('-Np', 
                    help = 'Number of points for discretization of phi', 
                    type = int, default = 50)
parser.add_argument('-kp', '--k_phi', 
                    help = 'Azimuthal wavenumbers; enter as -kp start stop step',
                    type = float, default = [1, 3, 1], nargs = 3)
parser.add_argument('-kz', '--k_z', 
                    help = 'DIMENSIONLESS vertical wavenumbers; enter as -kz start stop step',
                    type = float, default = [0, 2e-1, 3e-3], nargs = 3)
parser.add_argument('--modes', 
                    help = 'Number of modes of instability to be considered',
                    type = int, default = 1)
args = parser.parse_args()
                    
class Parameters:
    
    Ro   = args.Rossby
    Bu   = args.Burger
    f0   = args.Coriolis
    bkgd = args.bkgd

    Lr     = args.Lr #Max. r in physical space; half of computational domain
    Nr     = args.Neig #Number (odd) of gridpoints in computational domain
    halfNr = args.Neig // 2   
    Nφ     = args.Np #Number of azimuthal gridpoints; for visualization only
    kφs    = np.arange(args.k_phi[0], args.k_phi[1], args.k_phi[2])
    kzs    = np.arange(args.k_z[0], args.k_z[1], args.k_z[2])
        
    nmodes   = args.modes
    printout = args.PrintOutputs

    def display(self):
        print(f"Ro = {self.Ro}")
        print(f"Bu = {self.Bu}")
        print(f"f0 = {self.f0}")
        print(f"Lr = {self.Lr}")
        print(f"Nr = {self.Nr}")
        print(f"halfNr = {self.halfNr}")
        print(f"Np = {self.Np}")
        print(f"kps = {self.kφs}")
        print(f"kzs = {self.kzs}")
        print(f"nmodes = {self.nmodes}")

class Geometry:

    def __init__(self, params):
        
        self.method = "cheb"
           
        #Compute differentiation matrix and Chebyshev-spaced grid
        Dr, r = Chebyshev(params.Nr)
                        
        #Scale gridpoints and variable of differentiation to fit domain
        self.r, self.Dr = r * params.Lr, Dr / params.Lr
            
        self.Dr2 = np.matmul(self.Dr, self.Dr) #Second-order diff. matrix

def QG_Vortex_Stability():

    #Initialize parameters and set up geometry for Chebyshev solver
    params         = Parameters()
    geom           = Geometry(params)
    geom.Lap       = BuildLaplacian(params, geom)
    geom.rInterior = GridInterior(params, geom)
    
    #Discretize background-state-flow operators on Chebyshev grid
    geom.Ψ_op, geom.Q_op = BuildBkgdOperators(params, geom)
    
    #Information about wavenumbers and modes
    kφs, kzs, nmodes = params.kφs, params.kzs, params.nmodes

    #Initialize arrays to store results of eigen-computation
    growth = np.zeros([kzs.shape[0], kφs.shape[0], nmodes])
    prop   = np.zeros([kzs.shape[0], kφs.shape[0], nmodes])
    modes  = np.zeros([kzs.shape[0], kφs.shape[0], params.halfNr, nmodes],
                      dtype = complex)

    ########################################
    # SOLVE GENERALIZED EIGENVALUE PROBLEM #
    ########################################

    for kz_idx in range(0, kzs.shape[0]):

        kz  = kzs[kz_idx]
        kz2 = kz**2
  
        for kφ_idx in range(0, kφs.shape[0]):
            
            kφ  = kφs[kφ_idx]
            kφ2 = kφ**2
    
            ##############################
            # BUILD MATRICES 'A' and 'B' #
            ##############################
            
            B = (geom.Lap - np.diag(kφ2 / geom.rInterior**2)
                      - (kz2 * (1 / params.Bu) * np.eye(params.halfNr)))
            A = np.matmul(geom.Ψ_op, B) - geom.Q_op

            ##############################
            # FIND EIGENSPACE (DIRECTLY) #
            ##############################
            
            t0 = timeit.timeit()
        
            #Compute eigvals c and eigvecs psi with direct solver
            eigVals, eigVecs = spalg.eig(A, B)

            solveTime = timeit.timeit() - t0 #Time for direct solver
            
            #Indexing that sorts eigvals by ASCENDING Im(c)
            indSort = np.argsort(eigVals.imag)
            
            eigVals = eigVals[indSort] #Sort eigvals
            eigVecs = eigVecs[:, indSort] #Sort eigvecs in the same order
            ωs      = eigVals * kφ #Corresponding ω values for this kφ
            
            growth[kz_idx, kφ_idx, :]   = -ωs[0:nmodes].imag
            prop[kz_idx, kφ_idx, :]     = ωs[0:nmodes].real
            modes[kz_idx, kφ_idx, :, :] = eigVecs[:, 0:nmodes]

    #################
    # VISUALIZATION #
    #################

    plt.rcParams.update({"text.usetex": True, "font.size": 17})

    #Dimensionalize eigenvalues for visualization
    growthDim = growth * params.f0 * params.Ro
    propDim   = prop * params.f0 * params.Ro
    
    nkφ, nkz = (np.ravel(kφs)).shape[0], (np.ravel(kzs)).shape[0]

    for jj in range(0, nmodes):
        
        if nkφ < 4:

            ############################################################
            # VISUALIZE GROWTH RATES AND PROP. SPEEDS FOR DIFFERENT kφ #
            ############################################################
        
            fig, axes = plt.subplots(nkφ, 2, figsize = (13, 7), sharex = "col")

            for ii in range(0, nkφ):
                
                ax_growth = axes[ii, 0]
                ax_growth.plot(kzs, np.ravel(growthDim[:, ii, jj]), 
                               ".-", color = "mediumpurple")
                ax_growth.set(title = f"Growth rate; $k_{{\phi}}$ = {kφs[ii]}",
                              ylabel = "Growth rate ($s^{-1}$)")
                ax_growth.grid(True)

                ax_prop = axes[ii, 1]
                ax_prop.plot(kzs, np.ravel(propDim[:, ii, jj]), 
                             ".-", color = "mediumpurple")
                ax_prop.set(title = 
                        f"Propagation speed; $k_{{\phi}}$ = {kφs[ii]}",
                            ylabel = "Angular velocity ($s^{-1}$)")
                ax_prop.grid(True)

            ax_growth.set(xlabel = "Vertical wavenumber (per domain depth $H$)")
            ax_prop.set(xlabel = "Vertical wavenumber (per domain depth $H$)")
            
            if nmodes == 1:
                fig.savefig(f"omega_vs_m_fastestgrowing_nondimBTgyre.png")
            elif nmodes > 1:
                fig.savefig(f"omega_vs_m_mode{jj}_nondimBTgyre.png")

            plt.close(fig)

        ###########################################
        # PLOT EIGENFUNCTION STRUCTURES AGAINST r #
        ###########################################

        kz_idx, kφ_idx = 0, 1 #8, 0
        kz, kφ         = kzs[kz_idx], kφs[kφ_idx] #Wavenumbers to plot for

        eigVec     = modes[kz_idx, kφ_idx, :, jj]
        eigVecAmp  = np.sqrt(eigVec.real**2 + eigVec.imag**2)
        eigVecNorm = eigVec / max(eigVecAmp) #Normalize eigenvector

        fig, ax = plt.subplots(figsize = (10, 8))

        ax.plot(geom.r[1:(params.halfNr + 1)], eigVecNorm.real,
                "-", color = "mediumpurple",
                label = "Re[$\hat{\psi}$]")
        ax.plot(geom.r[1:(params.halfNr + 1)], eigVecNorm.imag,
                "--", color = "mediumpurple", 
                label = "Im[$\hat{\psi}$]")
        
        ax.set(xlabel = "$r/\sigma_r$", 
               ylabel = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$",
               title = fr"Components of fastest-growing eigenvector for wavenumbers $k_{{\phi}} =$ {kφ}, $\tilde{{m}} =$ {kz}")
        ax.legend()
        #plt.show()
        fig.savefig(f"eigvec_1Dstructure_k{kφ}_m{kz}_nondimBTgyre.png")
        plt.close(fig)
        
        ######################################
        # PLOT EIGEN-STRUCTURES IN r-φ PLANE #
        ######################################

        #Discretize phi-domain
        dφ      = 2 * pi / params.Nφ
        φCoords = dφ * np.arange(1, (params.Nφ + 1))

        #Meshgrid of polar coordinates to plot
        φVis, rVis = np.meshgrid(φCoords, geom.r[1:(params.halfNr + 1)])
        
        #Array to hold streamfunction values
        ψ = np.zeros([params.halfNr, params.Nφ], dtype = complex)

        #Evaluate streamfunction at (r, φ)-coordinate pairs
        for φ_idx in range(params.Nφ):
            for r_idx in range(params.halfNr - 1):
                ψ[r_idx, φ_idx] = Streamfunction(eigVecNorm[r_idx + 1],
                                                 k = kφ, φ = φCoords[φ_idx])

        #Evaluate components of eigen-velocity
        eigVecMesh, φMesh = np.meshgrid(eigVecNorm, φCoords)
        ur, uφ            = EigenvelocityFrom1DEigvec(params, geom, eigVecNorm,
                                                      kφ, φ = φMesh)

        #Absolute maximum amplitudes of velocity components
        urMax = np.max(np.abs(np.sqrt(ur.real**2 + ur.imag**2)))
        uφMax = np.max(np.abs(np.sqrt(uφ.real**2 + uφ.imag**2)))

        #Plot streamfunction in r-φ plane
        
        fig, axs = plt.subplots(1, 2, figsize = (11, 7),
                                subplot_kw = {"projection": "polar"})

        for i in range(2):
            axs[i].grid(False) #Required for pcolormesh

        axs[0].pcolormesh(φVis, rVis, ψ.real, 
                          cmap = "RdBu_r", vmin = -1, vmax = 1)
        axs[0].set_title(f"Re[$\hat{{\psi}}(r)$ exp($ik\phi$)]")
        axs[1].pcolormesh(φVis, rVis, ψ.imag, 
                          cmap = "RdBu_r", vmin = -1, vmax = 1)
        axs[1].set_title(f"Im[$\hat{{\psi}}(r)$ exp($ik\phi$)]")

        for i in range(2):
            axs[i].grid(True) #Restore grid for final version

        fig.subplots_adjust(hspace = 0.5, wspace = 0.75)
        fig.suptitle(f"Components of fastest-growing eigen-streamfunction in $r\phi$-plane\n" 
                     + fr"for wavenumbers $k_{{\phi}} =$ {kφ}, $\tilde{{m}} =$ {kz}")
        fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1), 
                                    cmap = "RdBu_r"), 
                     ax = axs.ravel().tolist(), orientation = "horizontal",
                     shrink = 0.8)
        #plt.show()
        fig.savefig(f"streamfunc2D_k{kφ}_m{kz}_nondimBTgyre.png")
        plt.close(fig)

        #Plot velocities in r-phi plane

        fig, axs = plt.subplots(2, 2, figsize = (8, 9),
                                subplot_kw = {"projection": "polar"})

        for i in range(2):
            for j in range(2):
                axs[i, j].grid(False) #Required for pcolormesh

        pcm_ur = axs[0, 0].pcolormesh(φVis, rVis, np.transpose(ur.real),
                                      cmap = "RdBu_r", vmin = -urMax, 
                                      vmax = urMax)
        axs[0, 0].set_title(f"Re[$u_r'(r, \phi)$]")
        axs[0, 1].pcolormesh(φVis, rVis, np.transpose(ur.imag),
                             cmap = "RdBu_r", vmin = -urMax, vmax = urMax)
        axs[0, 1].set_title(f"Im[$u_r'(r, \phi)$]")
        
        pcm_uφ = axs[1, 0].pcolormesh(φVis, rVis, np.transpose(uφ.real),
                                      cmap = "RdBu_r", vmin = -uφMax, 
                                      vmax = uφMax)
        axs[1, 0].set_title(f"Re[$u_{{\phi}}'(r,\phi)$]")
        axs[1, 1].pcolormesh(φVis, rVis, np.transpose(uφ.imag),
                             cmap = "RdBu_r", vmin = -uφMax, vmax = uφMax)
        axs[1, 1].set_title(f"Im[$u_{{\phi}}'(r,\phi)$]")

        for i in range(2):
            for j in range(2):
                axs[i, j].grid(True) #Restore grid for final version

        fig.subplots_adjust(hspace = 0.4, wspace = 0.8)
        fig.suptitle(f"Velocities derived from fastest-growing "
                     + "eigen-streamfunction \n in $r\phi$-plane "
                     + fr"for wavenumbers $k_{{\phi}} =$ {kφ}, $\tilde{{m}} =$ {kz}")
        fig.colorbar(pcm_ur, ax = [axs[0, 0], axs[0, 1]], 
                     location = "right", shrink = 0.6)
        fig.colorbar(pcm_uφ, ax = [axs[1, 0], axs[1, 1]], 
                     location = "right", shrink = 0.6)
        plt.show()
        fig.savefig(f"eigvelocities_k{kφ}_m{kz}_nondimBTgyre.png")
        plt.close(fig)

if __name__ == '__main__': #For testing
   QG_Vortex_Stability()
