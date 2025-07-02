"""`
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

Some of the notation follows "Spectral Methods in MATLAB" by L.N. Trefethen.
All variables are assumed to have been non-dimensionalized.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as nlg
import scipy
import scipy.linalg as spalg

import scipy.sparse as sp
#Note: I think we should aim to switch from matrix objects (e.g. built by sp.spdiags) to array objects, as the matrix functionality is now deprecated

import sys
import time
import timeit

from math import e, pi
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.sparse.linalg import eigs
from scipy.interpolate import interp1d
from scipy.special import factorial

from BuildLaplacian import BuildLaplacian
from BuildBkgdOperators import BuildBkgdOperators
from Chebyshev import Chebyshev
from FiniteDiff import FiniteDiff
from GetStreamfunc import GetStreamfunc

#Parse command-line inputs

parser = argparse.ArgumentParser()
parser.add_argument('--Neig', 
                    help = 'Number (must be ODD) of grid points for eig computations',
                    type = int, default = 181)
parser.add_argument('--Neigs', 
                    help = 'Number of grid points for eigs computations',
                    type = int, default = 1001)
parser.add_argument('-Lr', 
                    help = 'DIMENSIONLESS radius of the physical domain',
                    type = float, default = 5.0) #8.0)
parser.add_argument('-Ro', '--Rossby',
                    help = 'Rossby number of background flow', 
                    type = float, default = 4e-3)
parser.add_argument('-Bu', '--Burger',
                    help = 'Burger number of background flow',
                    type = float, default = 1.0) #2.5e-3)
parser.add_argument('-f0', '--Coriolis',
                    help = 'Coriolis frequency f0 (Hz)',
                    type = float, default = 1.4e-4)
parser.add_argument('--bkgd',
                    help = 'Background flow to use ("GM" or "BG")',
                    type = str, default = "GM") #"BG")
parser.add_argument('-p', '--PrintOutputs',
                    help = 'Flag to turn on display for each computation',
                    action = 'store_true')
parser.add_argument('-Np', 
                    help = 'Number (must be EVEN) of points for discretization of phi', 
                    type = int, default = 4)
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
    
    Np  = args.Np
    kps = np.arange(args.k_phi[0], args.k_phi[1], args.k_phi[2])
    kzs = np.arange(args.k_z[0], args.k_z[1], args.k_z[2])
        
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

    def __init__(self, method, params):
        
        self.method = method
        
        if method == "cheb":
           
            #Compute differentiation matrix and Chebyshev-spaced grid
            Dr, r = Chebyshev(params.Nr)
                        
            #Scale gridpoints and variable of differentiation to fit domain
            self.r, self.Dr = r * params.Lr, Dr / params.Lr
            
            self.Dr2 = np.matmul(self.Dr, self.Dr) #Second-order diff. matrix
        
        elif method == "FD":

            dr = 2 * params.Lr / params.Nr #Uniform r-interval size

            #Discretized domain with gridpoints listed in descending order 
            self.r = np.arange(params.Lr, (-params.Lr - dr), -dr)
                 
            #Compute (sparse) diff. matrix using 8th-order stencil
            self.Dr = FiniteDiff(self.r, 8, True, True)

            self.Dr2 = np.dot(self.Dr, self.Dr) #Second-order diff. matrix

def Print_npArray(fp, arr):
    for ii in xrange(0, arr.shape[0]):
        for jj in xrange(0, arr.shape[1]):
            if jj == (arr.shape[1] - 1):
                fp.write('{0:+2.2e}'.format(arr[ii, jj]))
            else:
                fp.write('{0:+2.2e}, '.format(arr[ii, jj]))
        fp.write('\n')
            
def QG_Vortex_Stability():
 
    #Initialize parameters and set up geometries for Chebyshev and FD solvers

    paramsCheb   = Parameters()
    GeomCheb     = Geometry("cheb", paramsCheb)
    GeomCheb.Lap = BuildLaplacian(paramsCheb, GeomCheb)

    paramsFD   = Parameters()
    GeomFD     = Geometry("FD", paramsFD)
    GeomFD.Lap = BuildLaplacian(paramsFD, GeomFD)

    #Discretize background-state-flow operators on Chebyshev and FD grids

    GeomCheb.rInterior, GeomCheb.Ψ_op, GeomCheb.Q_op = BuildBkgdOperators(
                                                        paramsCheb, GeomCheb)

    GeomFD.rInterior, GeomFD.Ψ_op, GeomFD.Q_op = BuildBkgdOperators(paramsFD,
                                                                    GeomFD)
    
    #Information about wavenumbers and modes is the same for both solvers
    kps, kzs, nmodes = paramsFD.kps, paramsFD.kzs, paramsFD.nmodes

    #Initialize arrays to store results of eigen-computations

    growthCheb = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    propCheb   = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    modesCheb  = np.zeros([kzs.shape[0], kps.shape[0], 
                           (paramsCheb.halfNr * paramsCheb.Np), nmodes],
                          dtype = complex)

    growthFD = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    propFD   = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    modesFD  = np.zeros([kzs.shape[0], kps.shape[0],
                         (paramsFD.halfNr * paramsFD.Np), nmodes],
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

            #For Chebyshev solver
            
            B_Cheb = (GeomCheb.Lap - np.diag(kp2 / GeomCheb.rInterior**2)
                      - (kz2 * (1 / paramsCheb.Bu) 
                         * np.eye(paramsCheb.halfNr * paramsCheb.Np)
                        )
                     )
            A_Cheb = (np.matmul(GeomCheb.Ψ_op, B_Cheb) - GeomCheb.Q_op)

            
            #For finite-difference solver
            B_FD = (GeomFD.Lap - sp.diags(kp2 / GeomFD.rInterior**2).tocsr()
                    - sp.diags([kz2 * (1/paramsFD.Bu)], 
                               shape = ((paramsFD.halfNr * paramsFD.Np), 
                                        (paramsFD.halfNr * paramsFD.Np))
                              ).tocsr()
                    )
            A_FD = GeomFD.Ψ_op.dot(B_FD) - GeomFD.Q_op
            
            ############################
            # FIND EIGENSPACE DIRECTLY #
            ############################
            
            t0Cheb = timeit.timeit()
        
            #Compute eigvals c and eigvecs psi with direct solver
            eigValCheb, eigVecCheb = spalg.eig(A_Cheb, B_Cheb)

            timeCheb = timeit.timeit() - t0Cheb #Time for direct Cheb solver
            
            #Indexing that sorts eigvals by ASCENDING Im(c)
            indCheb = np.argsort(eigValCheb.imag)
            
            eigValCheb = eigValCheb[indCheb] #Sort eigvals
            eigVecCheb = eigVecCheb[:, indCheb] #Sort eigvecs in the same order
            omegaCheb  = eigValCheb * kp #Corresponding omegas for this k_phi
            
            #Store results of direct Cheb solve

            growthCheb[kz_idx, kp_idx, :]   = -omegaCheb[0:nmodes].imag
            propCheb[kz_idx, kp_idx, :]     = omegaCheb[0:nmodes].real
            modesCheb[kz_idx, kp_idx, :, :] = eigVecCheb[:, 0:nmodes]
            
            t0FD = timeit.timeit()

            #Compute eigvals c and eigvecs psi with direct solver ('eig')
            eigValFD, eigVecFD = sp.linalg.eigs(A_FD, nmodes, B_FD)

            timeFD = timeit.timeit() - t0FD #Time for direct FD solver

            #Indexing that sorts eigvals by ASCENDING Im(c)
            indFD = np.argsort(eigValFD.imag)

            eigValFD = eigValFD[indFD] #Sort eigvals
            eigVecFD = eigVecFD[indFD, :] #Sort eigvecs in the same order
            omegaFD  = eigValFD * kp #Corresponding omegas for this k_phi

            #Store results of direct FD solve

            growthFD[kz_idx, kp_idx, :]   = -omegaFD[0:nmodes].imag
            propFD[kz_idx, kp_idx, :]     = omegaFD[0:nmodes].real
            modesFD[kz_idx, kp_idx, :, :] = eigVecFD[0:nmodes]
            
            ##############################
            # FIND EIGENSPACE INDIRECTLY #
            ##############################
            """
            for ii in range(0, nmodes): #Loop over modes
            
                growth = omegaCheb[ii].imag
                freq   = omegaCheb[ii].real

                sig0 = eigValCheb[ii]

                X = np.hstack([np.array([paramsCheb.Lr]),
                               np.ravel(GeomCheb.r[1:paramsCheb.halfNr+1]),
                               np.array([0])])[::-1]
                Y = np.hstack([np.array([0]), 
                               eigVecCheb[:, ii], 
                               np.array([0])])[::-1]

                Y_normalization = (-np.abs(Y)).argsort()
                Y               = Y / Y[Y_normalization[0]]
                
                Xnew = np.ravel(GeomFD.r[1:(paramsFD.halfNr + 1)])[1:-1]

                interp_fcn = interp1d(X, Y, kind = "cubic")
                chebvec    = interp_fcn(Xnew)
                chebvec    = chebvec[::-1]
                
                Xnew = Xnew[::-1]

                tmp           = chebvec
                tmp[tmp == 0] = 1
                tmp           = tmp.conj()
                
                T    = np.diag(np.ravel(tmp))
                Tinv = nlg.inv(T)

                t0 = timeit.timeit()
                
                try:

                    sig1, vec1 = eigs(np.dot(A_FD, Tinv), 1, np.dot(B_FD, Tinv),
                                      sigma = sig0, v0 = np.dot(T, chebvec))

                    vec1_normalization = (-np.abs(vec1)).argsort(axis = None)
                    vec1               = vec1 / vec1[vec1_normalization[0]]

                    #plt.subplot(3,2,1)
                    #plt.plot(Xnew,chebvec.real,'-b', Xnew,chebvec.imag,'-r')
                    #plt.title('Original Eig Vector')

                    #plt.subplot(3,2,2)
                    #plt.plot(Xnew, np.dot(T,chebvec).real, '-b', Xnew, np.dot(T,chebvec).imag, '-r')
                    #plt.title('Transformed Eig Vector')
                    
                    #plt.subplot(3,2,3)
                    #plt.plot(Xnew, vec1.real, '-b', Xnew, vec1.imag, '-r')
                    #plt.title('Original Eigs Vector')
                    
                    #vec1 = np.dot(Tinv, vec1)
                    #plt.subplot(3,2,4)
                    #plt.plot(Xnew, vec1.real, '-b', Xnew, vec1.imag, '-r')
                    #plt.title('Inverse Transformed Eigs Vector')

                    #plt.subplot(3,2,5)
                    #plt.plot(Xnew, np.abs(np.ravel(vec1)-np.ravel(chebvec)))
                    #plt.title('Absolute difference')
                    
                    #plt.show()

                except:
                    
                    sig1 = [np.nan + 1j * np.nan]

                    print('Eigs failed for mode {0:.2f}, k_phi = {1:.2f}, kz = {2:.4f}.\n'.format(ii, kp, kz))
                    sys.stdout.flush()
                
                timeFD = timeit.timeit() - t0
                
                omegaFD  = kp * sig1[0]

                growthFD_temp = omegaFD.imag
                freqFD_temp   = omegaFD.real

                growthFD[kz_idx, kp_idx, ii] = growthFD_temp
                freqFD[kz_idx, kp_idx, ii]   = freqFD_temp
        
                if paramsCheb.printout: #Display results
                    print('----------')
                    print('kz =i {0:4f}, kp = {1:2f}'.format(kz, kp))
                    print('eig : growth rate = {0:+4e}, frequency = {1:+4e}, cputime = {2:+4e}'.format(growth, freq, timeCheb))
                    print('eigs: growth rate = {0:+4e}, frequency = {1:+4e}, cputime = {2:+4e}'.format(growthFD_temp, freqFD_temp, timeFD))
                    sys.stdout.flush()
            """

    #################
    # VISUALIZATION #
    #################

    plt.rcParams.update({"text.usetex": True})

    #Dimensionalize eigenvalues
    
    growthDimCheb = growthCheb * paramsCheb.Ro * paramsCheb.f0
    propDimCheb   = propCheb * paramsCheb.Ro * paramsCheb.f0

    growthDimFD = growthFD * paramsFD.Ro * paramsFD.f0
    propDimFD   = propFD * paramsFD.Ro * paramsFD.f0
    
    nkp, nkz = (np.ravel(kps)).shape[0], (np.ravel(kzs)).shape[0]

    for jj in range(0, nmodes):
        
        if nkp < 4:

            #Visualize growth rates and propagation speeds for different kphi
        
            fig, axes = plt.subplots(nkp, 2, figsize = (10, 7), sharex = "col")
            
            for ii in range(0, nkp):
                
                ax_growth = axes[ii, 0]
                ax_growth.plot(kzs, 4 * np.ravel(growthDimCheb[:, ii, jj]), 
                               ".-", color = "mediumpurple", 
                               label = "Cheb solver")
                #ax_growth.plot(kzs, 4 * np.ravel(growthDimFD[:, ii, jj]), 
                #               ".-", color = "orange", label = "FD solver")
                ax_growth.set(title = 
                              f"Growth rate; $k_{{\phi}}$ = {kps[ii]}",
                              ylabel = "Growth rate ($s^{-1}$)")
                
                ax_prop = axes[ii, 1]
                ax_prop.plot(kzs, 4 * np.ravel(propDimCheb[:, ii, jj]), 
                             ".-", color = "mediumpurple", 
                             label = "Cheb solver")
                #ax_prop.plot(kzs, 4 * np.ravel(propDimFD[:, ii, jj]), 
                #             ".-", color = "orange", label = "FD solver")
                ax_prop.set(title = 
                        f"Propagation speed; $k_{{\phi}}$ = {kps[ii]}",
                            ylabel = "Azimuthal speed ($s^{-1}$)")

            ax_growth.set(xlabel = r'Vertical wavenumber ($\times \sigma_z$)')
            ax_prop.set(xlabel = r'Vertical wavenumber ($\times \sigma_z$)')
            axes[0, 0].legend()
            plt.show()
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

        kz_idx, kp_idx = 0, 1
        kz, kphi       = kzs[kz_idx], kps[kp_idx] #Wavenumbers to plot for

        dphi      = 2 * pi / paramsCheb.Np
        phiCoords = dphi * np.arange(1, (paramsCheb.Np + 1))

        [rVisCheb, phiVisCheb] = np.meshgrid(GeomCheb.r[0:(paramsCheb.halfNr 
                                                           + 1)], 
                                             phiCoords)
        [rVisFD, phiVisFD]     = np.meshgrid(GeomFD.r[0:(paramsFD.halfNr + 1)],
                                             phiCoords)

        eigvecCheb    = modesCheb[kz_idx, kp_idx, :, jj]
        eigvecChebAmp = np.sqrt(eigvecCheb.real**2 + eigvecCheb.imag**2)

        #Reshape output for visualization

        tmp1       = np.reshape(eigvecCheb, 
                                (paramsCheb.halfNr, paramsCheb.Np)
                               ).T
        tmp2       = np.vstack((tmp1[(paramsCheb.Np - 1), :], 
                                tmp1[0:(paramsCheb.Np - 1), :]))
        eigvecCheb = np.hstack([np.zeros((paramsCheb.Np, 1)), tmp2])

        #Conjugate (because we used '.T') and normalize eigenvector
        eigvecChebNorm = eigvecCheb.conj() / max(eigvecChebAmp)

        #Is there a formula I can use to get this a priori?
        indEigvecCheb = np.argsort(
                           np.sum(abs(eigvecChebNorm), axis = 1))[-1]

        eigvecFD    = modesFD[kz_idx, kp_idx, :, jj]
        eigvecFDAmp = np.sqrt(eigvecFD.real**2 + eigvecFD.imag**2)
        
        tmp1     = np.reshape(eigvecFD, (paramsFD.halfNr, paramsFD.Np)).T
        tmp2     = np.vstack((tmp1[(paramsFD.Np - 1), :], 
                              tmp1[0:(paramsFD.Np - 1), :]))
        eigvecFD = np.hstack([np.zeros((paramsFD.Np, 1)), tmp2])

        #Conjugate (because we used '.T') and normalize eigenvector
        eigvecFDNorm = eigvecFD.conj() / max(eigvecFDAmp)
         
        fig, ax = plt.subplots(figsize = (10, 8))

        ax.plot(rVisCheb[1, :], eigvecChebNorm[indEigvecCheb, :].real,
                "-", color = "mediumpurple",
                label = "Re[$\hat{\psi}$]; Cheb solver")
        ax.plot(rVisCheb[1, :], eigvecChebNorm[indEigvecCheb, :].imag,
                "--", color = "mediumpurple", 
                label = "Im[$\hat{\psi}$]; Cheb solver")
        #ax.plot(rVisFD[1, :], eigvecFDNorm[-1, :].real,
        #        "-", color = "orange", label = "Re[$\hat{\psi}$]; FD solver")
        #ax.plot(rVisFD[1, :], eigvecFDNorm[-1, :].imag,
        #        "--", color = "orange", label = "Im[$\hat{\psi}$]; FD solver")
        
        ax.set(xlabel = "$r/\sigma_r$", 
               ylabel = "Component of $\hat{\psi}$, normalized by max. amplitude of $\hat{\psi}$",
               title = f"Components of mode-{jj} eigenvector for wavenumbers $k_{{\phi}}$ = {kphi}, $m =$ {kz}")
        ax.legend()
        plt.show()
        fig.savefig(f"eigvec_structure_k{kphi}_m{kz}_mode{jj}.png")
        plt.close(fig)
        
        #Plot streamfunction structures in r-phi plane
        
        #Meshgrids of polar coordinates to plot
        
        phiVisCheb, rVisCheb = np.meshgrid(phiCoords, 
                                    GeomCheb.r[1:(paramsCheb.halfNr + 1)])
        phiVisFD, rVisFD     = np.meshgrid(phiCoords, 
                                    GeomFD.r[0:(paramsFD.halfNr + 1)])
        
        #Arrays to hold streamfunction values
        
        psiCheb = np.zeros([paramsCheb.halfNr, paramsCheb.Np], 
                           dtype = complex)
        psiFD   = np.zeros([(paramsFD.halfNr + 1), paramsCheb.Np], 
                           dtype = complex)
        
        #Evaluate streamfunction at (r, phi)-coordinate pairs
        
        for phi_idx in range(paramsCheb.Np):
            for r_idx in range(paramsCheb.halfNr):
                psiCheb[r_idx, phi_idx] = GetStreamfunc(eigvecChebNorm[indEigvecCheb, 
                                                                       r_idx+1],
                                            k = kphi, phi = phiCoords[phi_idx])
            for r_idx in range(paramsFD.halfNr):
                psiFD[r_idx, phi_idx] = GetStreamfunc(eigvecFDNorm[-1, r_idx], 
                                            k = kphi, phi = phiCoords[phi_idx])

        #psiChebAmp  = np.sqrt(psiCheb.real**2 + psiCheb.imag**2)
        #psiChebNorm = psiCheb / max(psiChebAmp)

        fig, axs = plt.subplots(2, 2, figsize = (8, 10),
                                subplot_kw = dict(projection = "polar"))

        for i in range(2):
            for j in range(2):
                axs[i, j].grid(False) #Required for pcolormesh

        axs[0, 0].pcolormesh(phiVisCheb, rVisCheb, psiCheb.real, 
                             cmap = "bwr", vmin = -1, vmax = 1)
        axs[0, 0].set(title = f"Re[$\hat{{\psi}}(r)$ exp($ik\phi$)]; Cheb solver")

        axs[0, 1].pcolormesh(phiVisCheb, rVisCheb, psiCheb.imag, 
                             cmap = "bwr", vmin = -1, vmax = 1)
        axs[0, 1].set(title = f"Im[$\hat{{\psi}}(r)$ exp($ik\phi$)]; Cheb solver")

        axs[1, 0].pcolormesh(phiVisFD, rVisFD, psiFD.real,
                             cmap = "bwr", vmin = -1, vmax = 1)
        axs[1, 0].set(title = f"Re[$\hat{{\psi}}(r)$ exp($ik\phi$)]; FD solver")

        axs[1, 1].pcolormesh(phiVisFD, rVisFD, psiFD.imag, 
                             cmap = "bwr", vmin = -1, vmax = 1)
        axs[1, 1].set(title = f"Im[$\hat{{\psi}}(r)$ exp($ik\phi$)]; FD solver")

        for i in range(2):
            for j in range(2):
                axs[i, j].grid(True) #Restore grid for final version

        fig.subplots_adjust(hspace = 0.5, wspace = 0.1)
        fig.suptitle(f"Components of mode-{jj} eigen-streamfunction in $r\phi$-plane for wavenumbers $k_{{\phi}}$ = {kphi}, $m =$ {kz}")
        fig.colorbar(ScalarMappable(norm = Normalize(vmin = -1, vmax = 1), 
                                    cmap = "bwr"), 
                     ax = axs.ravel().tolist(), orientation = "horizontal",
                     shrink = 0.75)
        plt.show()
        fig.savefig(f"streamfunc_2d_k{kphi}_m{kz}_mode{jj}.png")
        plt.close(fig)
   
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
