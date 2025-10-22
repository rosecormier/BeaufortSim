"""
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

Some of the notation follows "Spectral Methods in MATLAB" by L.N. Trefethen.
All variables should be specified in SI base units.
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
from BuildBkgdOperators import BuildBkgdOperators, rInterior
from Chebyshev import Chebyshev
from Streamfunctions import Streamfunction, EigenvelocityFromEigvec

parser = argparse.ArgumentParser()
parser.add_argument('--Neig', 
                    help = 'Number (must be ODD) of grid points for eig computations',
                    type = int, default = 801)
parser.add_argument('-Lr', 
                    help = 'DIMENSIONAL radius of the physical domain (m)',
                    type = float, default = 2e6)
parser.add_argument('-N', '--buoyancyfreq',
                    help = 'Buoyancy frequency (Hz)',
                    type = float, default = 1e-2)
parser.add_argument('-f0', '--Coriolis',
                    help = 'Coriolis frequency f0 (Hz)',
                    type = float, default = 1.4e-4)
parser.add_argument('-U', '--bkgdU',
                    help = 'Characteristic scale for background velocity (m/s)',
                    type = float, default = 1.5e-1)
parser.add_argument('-σr', '--sigmar',
                    help = 'Radial length scale of gyre (m)',
                    type = float, default = 2.5e5)
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
                    help = 'DIMENSIONAL vertical wavenumbers (m^{-1}); enter as -kz start stop step',
                    type = float, default = [0, 4e-4, 5e-6], nargs = 3)
parser.add_argument('--modes', 
                    help = 'Number of modes of instability to be considered',
                    type = int, default = 1)
args = parser.parse_args()
                    
class Parameters:
    
    N    = args.buoyancyfreq
    f0   = args.Coriolis
    U    = args.bkgdU
    σr   = args.sigmar
    bkgd = args.bkgd

    Ro = U / (σr * f0) #Rossby number

    Lr     = args.Lr #Max. r in physical space; half of computational domain
    Nr     = args.Neig #Number (odd) of gridpoints in computational domain
    halfNr = args.Neig // 2   
    Nφ     = args.Np #Number of azimuthal gridpoints; for visualization only
    kφs    = np.arange(args.k_phi[0], args.k_phi[1], args.k_phi[2])
    kzs    = np.arange(args.k_z[0], args.k_z[1], args.k_z[2])
        
    nmodes   = args.modes
    printout = args.PrintOutputs

    def display(self):
        print(f"N = {self.N}")
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
    geom.rInterior = rInterior(params, geom)
    
    #Discretize background-state-flow operators on Chebyshev grid
    geom.Ψ_op, geom.Q_op = BuildBkgdOperators(params, geom, 
                                              dimensional_U = params.U,
                                              dimensional_σr = params.σr)
     
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
                 - (kz2 * (params.f0**2 / params.N**2) * np.eye(params.halfNr)))
            A = (np.matmul(geom.Ψ_op, B) - geom.Q_op)

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
    
    nkφ, nkz = (np.ravel(kφs)).shape[0], (np.ravel(kzs)).shape[0]

    for jj in range(0, nmodes):
        
        if nkφ < 4:

            #Visualize growth rates and propagation speeds for different kphi
        
            fig, axes = plt.subplots(nkφ, 2, figsize = (15, 7), sharex = "col")

            for ii in range(0, nkφ):
                
                ax_growth = axes[ii, 0]
                ax_growth.plot(kzs, np.ravel(growth[:, ii, jj]), 
                               ".-", color = "mediumpurple")
                ax_growth.set(title = 
                              f"Growth rate; $k_{{\phi}}$ = {kφs[ii]}",
                              ylabel = "Growth rate (s$^{{-1}}$)")
                ax_growth.grid(True)

                ax_prop = axes[ii, 1]
                ax_prop.plot(kzs, np.ravel(prop[:, ii, jj]), 
                             ".-", color = "mediumpurple")
                ax_prop.set(title = 
                        f"Propagation speed; $k_{{\phi}}$ = {kφs[ii]}",
                            ylabel = "Azimuthal speed (s$^{{-1}}$)")
                ax_prop.grid(True)

            ax_growth.set(xlabel = r'Vertical wavenumber (m$^{{-1}}$)')
            ax_prop.set(xlabel = r'Vertical wavenumber (m$^{{-1}}$)')
            #plt.show()
            fig.savefig(f"omega_vs_m_mode{jj}_dimensionalBTgyre.png")
            plt.close(fig)
        
        """
        elif nkz < 4:

            for ii in range(0, nkz):
                
                plt.subplot(nkz, 2, (1 + 2*ii))
                plt.plot(np.ravel(kps), 4 * np.ravel(growthFD[ii, :, jj]), '-o',
                         np.ravel(kps), 4 * np.ravel(growthCheb[ii, :, jj]), '-*')
                plt.title('Growth rate for vertical wavenumber = {} (units?)'.format(ii))
                plt.xlabel('Azimuthal wavenumber')
                plt.ylabel('Growth rate (units?)')

                plt.subplot(nkz, 2, (2 + 2*ii))
                plt.plot(np.ravel(kps), 4 * np.ravel(freqFD[ii, :, jj]), '-o',
                         np.ravel(kps), 4 * np.ravel(freqCheb[ii, :, jj]), '-*')
                plt.title('Propagation speed for vertical wavenumber = {} (units?)'.format(ii))
                plt.xlabel('Azimuthal wavenumber')
                plt.ylabel('Propagation speed (units?)')

        else:

            plt.subplot(2, 2, 1)
            plt.contour(np.ravel(kps), np.ravel(kzs), 4 * growthFD[:, :, jj])
            plt.title('Growth rate (eigs)')

            plt.subplot(2, 2, 2)
            plt.contour(np.ravel(kps), np.ravel(kzs), 4 * freqFD[:, :, jj])
            plt.title('Propagation speed (eigs)')

            plt.subplot(2, 2, 3)
            plt.contour(np.ravel(kps), np.ravel(kzs), 4 * growthFD[:, :, jj])
            plt.title('Growth rate (eig)')

            plt.subplot(2, 2, 4)
            plt.contour(np.ravel(kps), np.ravel(kzs), 4 * freqFD[:, :, jj])
            plt.title('Propagation speed (eig)')
        """

        #Plot eigenfunction structures against r

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
        
        ax.set(xlabel = "$r$ (m)", 
               ylabel = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$",
               title = f"Components of fastest-growing eigenvector for wavenumbers $k_{{\phi}}$ = {kφ}, $m =$ {kz}")
        ax.legend()
        #plt.show()
        fig.savefig(f"eigvec_1Dstructure_k{kφ}_m{kz}_dimensionalBTgyre.png")
        plt.close(fig)
        
        #Plot streamfunction structures in r-φ plane
        
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

        fig, axs = plt.subplots(1, 2, figsize = (11, 7),
                                subplot_kw = {"projection": "polar"})

        for i in range(2):
            axs[i].grid(False) #Required for pcolormesh

        axs[0].pcolormesh(φVis, rVis, ψ.real, 
                             cmap = "RdBu_r", vmin = -1, vmax = 1)
        axs[0].set(title = f"Re[$\hat{{\psi}}(r)$ exp($ik\phi$)]")

        axs[1].pcolormesh(φVis, rVis, ψ.imag,
                             cmap = "RdBu_r", vmin = -1, vmax = 1)
        axs[1].set(title = f"Im[$\hat{{\psi}}(r)$ exp($ik\phi$)]")

        for i in range(2):
            axs[i].grid(True) #Restore grid for final version

        fig.subplots_adjust(hspace = 0.5, wspace = 0.75)
        fig.suptitle(f"Components of fastest-growing eigen-streamfunction in $r\phi$-plane\n for wavenumbers $k_{{\phi}}$ = {kφ}, $m =$ {kz} m$^{{-1}}$\n\n")
        fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1),
                                    cmap = "RdBu_r"), 
                     ax = axs.ravel().tolist(), orientation = "horizontal",
                     shrink = 0.8)
        #plt.show()
        fig.savefig(f"streamfunc2D_k{kφ}_m{kz}_dimensionalBTgyre.png")
        plt.close(fig)

if __name__ == '__main__': #For testing
   QG_Vortex_Stability()
