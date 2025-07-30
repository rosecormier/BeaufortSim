"""`
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

Some of the notation follows "Spectral Methods in MATLAB" by L.N. Trefethen.
All variables are assumed to have been non-dimensionalized.
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
from Streamfunctions import GetStreamfunc

parser = argparse.ArgumentParser()
parser.add_argument('--Neig', 
                    help = 'Number (must be ODD) of grid points for eig computations',
                    type = int, default = 181)
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
                    type = float, default = [0, 2, 0.1], nargs = 3)
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
    Np     = args.Np #Number of azimuthal gridpoints; for visualization only
    kps    = np.arange(args.k_phi[0], args.k_phi[1], args.k_phi[2])
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
        print(f"kps = {self.kps}")
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

def Print_npArray(fp, arr):
    for ii in xrange(0, arr.shape[0]):
        for jj in xrange(0, arr.shape[1]):
            if jj == (arr.shape[1] - 1):
                fp.write('{0:+2.2e}'.format(arr[ii, jj]))
            else:
                fp.write('{0:+2.2e}, '.format(arr[ii, jj]))
        fp.write('\n')
            
def QG_Vortex_Stability():

    ##########
    # SET-UP #
    ##########

    #Initialize parameters and set up geometry for Chebyshev solver

    paramsCheb         = Parameters()
    GeomCheb           = Geometry(paramsCheb)
    GeomCheb.Lap       = BuildLaplacian(paramsCheb, GeomCheb)
    GeomCheb.rInterior = rInterior(paramsCheb, GeomCheb)
    
    #Discretize background-state-flow operators on Chebyshev grid
    GeomCheb.Ψ_op, GeomCheb.Q_op = BuildBkgdOperators(paramsCheb, GeomCheb)
    
    #Information about wavenumbers and modes
    kps, kzs, nmodes = paramsCheb.kps, paramsCheb.kzs, paramsCheb.nmodes

    #Initialize arrays to store results of eigen-computation

    growthCheb = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    propCheb   = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    modesCheb  = np.zeros([kzs.shape[0], kps.shape[0], paramsCheb.halfNr, 
                           nmodes],
                          dtype = complex)

    ########################################
    # SOLVE GENERALIZED EIGENVALUE PROBLEM #
    ########################################

    for kz_idx in range(0, kzs.shape[0]):

        kz  = kzs[kz_idx]
        kz2 = kz**2
  
        for kp_idx in range(0, kps.shape[0]):

            kp  = kps[kp_idx]
            kp2 = kp**2
    
            ##############################
            # BUILD MATRICES 'A' and 'B' #
            ##############################
            
            B_Cheb = (GeomCheb.Lap - np.diag(kp2 / GeomCheb.rInterior**2)
                      - (kz2 * (1 / paramsCheb.Bu) 
                         * np.eye(paramsCheb.halfNr)
                        )
                     )
            A_Cheb = (np.matmul(GeomCheb.Ψ_op, B_Cheb) - GeomCheb.Q_op)

            ##############################
            # FIND EIGENSPACE (DIRECTLY) #
            ##############################
            
            t0Cheb = timeit.timeit()
        
            #Compute eigvals c and eigvecs psi with direct solver
            eigValCheb, eigVecCheb = spalg.eig(A_Cheb, B_Cheb)

            timeCheb = timeit.timeit() - t0Cheb #Time for direct Cheb solver
            
            #Indexing that sorts eigvals by ASCENDING Im(c)
            indCheb = np.argsort(eigValCheb.imag)
            
            eigValCheb = eigValCheb[indCheb] #Sort eigvals
            eigVecCheb = eigVecCheb[:, indCheb] #Sort eigvecs in the same order
            omegaCheb  = eigValCheb * kp #Corresponding omegas for this k_phi
            
            #Store results

            growthCheb[kz_idx, kp_idx, :]   = -omegaCheb[0:nmodes].imag
            propCheb[kz_idx, kp_idx, :]     = omegaCheb[0:nmodes].real
            modesCheb[kz_idx, kp_idx, :, :] = eigVecCheb[:, 0:nmodes]

    #################
    # VISUALIZATION #
    #################

    plt.rcParams.update({"text.usetex": True,
                         "font.size": 17})

    #Dimensionalize eigenvalues
    
    growthDimCheb = growthCheb * paramsCheb.Ro * paramsCheb.f0
    propDimCheb   = propCheb * paramsCheb.Ro * paramsCheb.f0
    
    nkp, nkz = (np.ravel(kps)).shape[0], (np.ravel(kzs)).shape[0]

    for jj in range(0, nmodes):
        
        if nkp < 4:

            #Visualize growth rates and propagation speeds for different kphi
        
            fig, axes = plt.subplots(nkp, 2, figsize = (10, 7), sharex = "col")
            fig_poster, ax_poster = plt.subplots(1, 1, figsize = (3, 4))

            for ii in range(0, nkp):
                
                ax_growth = axes[ii, 0]
                ax_growth.plot(kzs, 4 * np.ravel(growthDimCheb[:, ii, jj]), 
                               ".-", color = "mediumpurple", 
                               label = "Cheb solver")
                ax_growth.set(title = 
                              f"Growth rate; $k_{{\phi}}$ = {kps[ii]}",
                              ylabel = "Growth rate ($s^{-1}$)")
                
                ax_prop = axes[ii, 1]
                ax_prop.plot(kzs, 4 * np.ravel(propDimCheb[:, ii, jj]), 
                             ".-", color = "mediumpurple", 
                             label = "Cheb solver")
                ax_prop.set(title = 
                        f"Propagation speed; $k_{{\phi}}$ = {kps[ii]}",
                            ylabel = "Azimuthal speed ($s^{-1}$)")

            ax_poster.plot(kzs, np.ravel(growthDimCheb[:, 0, 0]),
                           "-", color = "#f49100",
                           label = "k = 1")
            ax_poster.plot(kzs, np.ravel(growthDimCheb[:, 1, 0]),
                           "-", color = "#0d82a8",
                           label = "k = 2")
            ax_poster.set(xlabel = r"$m$ (vertical wavenumber) $\times H$",
                          ylabel = "Growth rate ($s^{-1}$)",
                          title = "Growth rates of fastest-growing\n modes")
            plt.grid(True)
            ax_poster.legend()
            fig_poster.savefig("omega_vs_m_poster.pdf")
            plt.close(fig_poster)

            ax_growth.set(xlabel = r'Vertical wavenumber ($\times \sigma_z$)')
            ax_prop.set(xlabel = r'Vertical wavenumber ($\times \sigma_z$)')
            axes[0, 0].legend()
            #plt.show()
            fig.savefig(f"omega_vs_m_mode{jj}.png")
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

        kz_idx, kp_idx = 8, 0 #0, 1
        kz, kphi       = kzs[kz_idx], kps[kp_idx] #Wavenumbers to plot for

        eigvecCheb     = modesCheb[kz_idx, kp_idx, :, jj]
        eigvecChebAmp  = np.sqrt(eigvecCheb.real**2 + eigvecCheb.imag**2)
        eigvecChebNorm = eigvecCheb / max(eigvecChebAmp) #Normalize eigenvector
     
        fig, ax = plt.subplots(figsize = (10, 8))

        ax.plot(GeomCheb.r[1:(paramsCheb.halfNr + 1)], eigvecChebNorm.real,
                "-", color = "mediumpurple",
                label = "Re[$\hat{\psi}$]; Cheb solver")
        ax.plot(GeomCheb.r[1:(paramsCheb.halfNr + 1)], eigvecChebNorm.imag,
                "--", color = "mediumpurple", 
                label = "Im[$\hat{\psi}$]; Cheb solver")
        
        ax.set(xlabel = "$r/\sigma_r$", 
               ylabel = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$",
               title = f"Components of mode-{jj} eigenvector for wavenumbers $k_{{\phi}}$ = {kphi}, $m =$ {kz}")
        ax.legend()
        #plt.show()
        fig.savefig(f"eigvec_structure_k{kphi}_m{kz}_mode{jj}.png")
        plt.close(fig)
        
        #Plot streamfunction structures in r-phi plane
        
        dphi      = 2 * pi / paramsCheb.Np
        phiCoords = dphi * np.arange(1, (paramsCheb.Np + 1))

        #Meshgrid of polar coordinates to plot
        phiVisCheb, rVisCheb = np.meshgrid(phiCoords, 
                                    GeomCheb.r[1:(paramsCheb.halfNr + 1)])
        
        #Array to hold streamfunction values
        psiCheb = np.zeros([paramsCheb.halfNr, paramsCheb.Np], 
                           dtype = complex)
        
        #Evaluate streamfunction at (r, phi)-coordinate pairs
        
        for phi_idx in range(paramsCheb.Np):
            for r_idx in range(paramsCheb.halfNr - 1):
                psiCheb[r_idx, phi_idx] = GetStreamfunc(
                                            eigvecChebNorm[r_idx + 1],
                                            k = kphi, phi = phiCoords[phi_idx])

        fig, axs = plt.subplots(1, 2, figsize = (4, 6),
                                subplot_kw = dict(projection = "polar"))

        for i in range(2):
            axs[i].grid(False) #Required for pcolormesh

        axs[0].pcolormesh(phiVisCheb, rVisCheb, psiCheb.real, 
                             cmap = "bwr", vmin = -1, vmax = 1)
        axs[0].set(title = f"Re[$\hat{{\psi}}(r)$ exp($ik\phi$)]; Cheb solver")

        axs[1].pcolormesh(phiVisCheb, rVisCheb, psiCheb.imag, 
                             cmap = "bwr", vmin = -1, vmax = 1)
        axs[1].set(title = f"Im[$\hat{{\psi}}(r)$ exp($ik\phi$)]; Cheb solver")

        for i in range(2):
            axs[i].grid(True) #Restore grid for final version

        fig.subplots_adjust(hspace = 0.5, wspace = 0.1)
        fig.suptitle(f"Components of mode-{jj} eigen-streamfunction in $r\phi$-plane for wavenumbers $k_{{\phi}}$ = {kphi}, $m =$ {kz}")
        fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1), 
                                    cmap = "bwr"), 
                     ax = axs.ravel().tolist(), orientation = "horizontal",
                     shrink = 0.75)
        #plt.show()
        fig.savefig(f"streamfunc_2d_k{kphi}_m{kz}_mode{jj}.png")
        plt.close(fig)
        
        fig_poster, ax_poster = plt.subplots(1, 1, figsize = (4, 7),
                                subplot_kw = dict(projection = "polar"))
        ax_poster.grid(False)
        ax_poster.pcolormesh(phiVisCheb, rVisCheb, psiCheb.real, cmap = "RdBu_r", vmin = -1, vmax = 1)
        #ax_poster.set_title(f"$k=$ {kphi}; $m=$ {kz}")
        ax_poster.grid(True)

        fig_poster.suptitle(f"Re [$\hat{{\psi}}(r)$ exp$(ik\phi)$];\n fastest-growing mode with\n $k=$ {kphi} and $m=$ {kz}\n\n")
        fig_poster.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1), cmap = "RdBu_r"),
                     ax = ax_poster, orientation = "horizontal") #, shrink = 0.75)
        plt.show()
        fig_poster.savefig(f"poster_streamfunc_k{kphi}_m{kz}.pdf")
        plt.close(fig_poster)

        """
        xx, yy = rVisCheb * np.cos(phiVisCheb), rVisCheb * np.sin(phiVisCheb)
        
        fig, axs = plt.subplots(1, 2, figsize = (9, 5), 
                                subplot_kw = dict(projection = "3d"))

        axs[0].plot_surface(xx, yy, psiCheb.real, 
                            rstride = 2, cstride = 2, cmap = "RdYlBu_r", 
                            vmin = -1, vmax = 1, alpha = 0.8)
        axs[1].plot_surface(xx, yy, psiCheb.imag,
                            rstride = 2, cstride = 2, cmap = "RdYlBu_r",
                            vmin = -1, vmax = 1, alpha = 0.8)
        
        for i in range(2):
            axs[i].set_xlabel("x")
            axs[i].set_ylabel("y")
            axs[i].set_xlim(-paramsCheb.Lr, paramsCheb.Lr)
            axs[i].set_ylim(-paramsCheb.Lr, paramsCheb.Lr)
        
        axs[0].set_zlabel(f"Re[$\hat{{\psi}}(r)$]; Cheb solver")
        axs[1].set_zlabel(f"Im[$\hat{{\psi}}(r)$]; Cheb solver")

        fig.subplots_adjust(hspace = 0.5, wspace = 0.1)
        fig.suptitle(f"Components of mode-{jj} eigen-streamfunction in $r\phi$-plane for wavenumbers $k_{{\phi}}$ = {kphi}, $m =$ {kz}")
        fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1),
                                    cmap = "RdYlBu_r"),
                     ax = axs.ravel().tolist(), orientation = "horizontal",
                     shrink = 0.75)
        plt.show()
        fig.savefig(f"streamfunc_surface_k{kphi}_m{kz}_mode{jj}.png")
        plt.close(fig)
        """
if __name__ == '__main__': #For testing
   QG_Vortex_Stability()
