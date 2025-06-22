"""
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
import sys
import time
import timeit

from math import e
from scipy.sparse.linalg import eigs
from scipy.interpolate import interp1d
from scipy.special import factorial

from Chebyshev import Chebyshev
from FiniteDiff import FiniteDiff
from GetStreamfunc import GetStreamfunc

#Parse command-line inputs

parser = argparse.ArgumentParser()
parser.add_argument('--Neig', 
                    help = 'Number of grid points for eig computations',
                    type = int, default = 181)
parser.add_argument('--Neigs', 
                    help = 'Number of grid points for eigs computations',
                    type = int, default = 1001)
parser.add_argument('-Lr', 
                    help = 'DIMENSIONLESS radius of the physical domain',
                    type = float, default = 8.0) #6.25)
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
parser.add_argument('-kp', '--k_phi', 
                    help = 'Azimuthal wavenumbers; Enter as -kp min max step',
                    type = float, default = [1, 4, 1], nargs = 3)
parser.add_argument('-kz', '--k_z', 
                    help = 'DIMENSIONLESS vertical wavenumbers; Enter as -kz min max step',
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
    Nr     = args.Neig #Number of gridpoints
    halfNr = args.Neig // 2
    
    Nt  = 40
    kps = np.arange(args.k_phi[0], args.k_phi[1], args.k_phi[2])
    kzs = np.arange(args.k_z[0], args.k_z[1], args.k_z[2])
        
    nmodes   = args.modes
    printout = args.PrintOutputs

    def display(self):
        print('Ro = {}'.format(self.Ro))
        print('Bu = {}'.format(self.Bu))
        print('Lr = {}'.format(self.Lr))
        print('Nr = {}'.format(self.Nr))
        print('halfNr = {}'.format(self.halfNr))
        print('Nt = {}'.format(self.Nt))
        print('kps = {}'.format(self.kps))
        print('kzs = {}'.format(self.kzs))
        print('nmodes = {}'.format(self.nmodes))

class Geometry:

    def __init__(self, method, params):
        
        self.method = method
        
        if method == 'cheb':
           
            #Compute differentiation matrix (Dr) and Chebyshev-spaced grid (r)
            Dr, r = Chebyshev(params.Nr)
            
            #Scale gridpoints and variable of differentiation to fit domain
            self.r, self.Dr = r * params.Lr, Dr / params.Lr
            
            #Second-order differentiation matrix
            self.Dr2 = np.matmul(self.Dr, self.Dr)
        
        elif method == 'FD':

            dr = 2 * params.Lr / params.Nr #Uniform r-interval size

            #Discretized domain with gridpoints listed in descending order 
            self.r = np.arange(params.Lr, (-params.Lr - dr), -dr)
            
            #Nonzero entries of (sparse) differentiation matrix using 8th-order 
            # stencil. Size of matrix is len(r) x len(r).
            self.Dr = FiniteDiff(self.r, 8, True, True)

            #Second-order differentiation matrix
            self.Dr2 = np.dot(self.Dr, self.Dr)

def Build_Laplacian(params, geom):
   
    halfNr, Nr = params.halfNr, params.Nr

    #Note: we do not explicitly set zeroth indices because we impose zero
    # Dirichlet BCs

    #Quadrants of 2nd-order r-derivative matrix to be retained
    D1 = geom.Dr2[1:halfNr+1, :][:, 1:halfNr+1] #(pos, pos)
    D2 = geom.Dr2[1:halfNr+1, :][:, np.arange(Nr-1, halfNr, -1)] #(pos, neg)

    #Quadrants of 1st-order r-derivatives matrix to be retained
    E1 = geom.Dr[1:halfNr+1, :][:, 1:halfNr+1] #(pos, pos)
    E2 = geom.Dr[1:halfNr+1, :][:, np.arange(Nr-1, halfNr, -1)] #(pos, neg)

    #Build diagonal matrix from reciprocals of r_j for 1 <= j <= halfNr
    if sp.issparse(geom.Dr):
        R = sp.spdiags(np.transpose(1 / geom.r[1:(halfNr+1)]), 
                       np.array([0]), halfNr, halfNr)
    else:
        R = np.diag(1 / np.ravel(geom.r[1:(halfNr+1)]))

    #Build discretized Laplacian as done in Ch.11 of Trefethen
    Laplacian = D1 + D2 + np.dot(R, E1 + E2)

    return Laplacian

def Build_Bkgd_Operators(params, geom):

    #Array of those r-values at interior gridpoints that lie in physical space
    rInterior = geom.r[1:(params.halfNr + 1)]

    #Build array of (1/r) * (dPsi/dr) evaluated at gridpoints
    # and array of (1/r) * (dQ/dr) evaluated at gridpoints

    if params.bkgd == "GM":
        Ψ_op = np.ravel(-0.5 * np.exp(-rInterior**2))
        Q_op = np.ravel(-2 * np.exp(-rInterior**2) * (rInterior**2 - 2))

    elif params.bkgd == "BG":
        Ψ_op = np.ravel(np.sqrt(2*e) * np.exp(-(rInterior**2)))
        Q_op = np.ravel(np.sqrt(32*e) * (rInterior**2 - 2) 
                        * np.exp(-rInterior**2))

    return Ψ_op, Q_op

def Print_npArray(fp, arr):
    for ii in xrange(0, arr.shape[0]):
        for jj in xrange(0, arr.shape[1]):
            if jj == (arr.shape[1] - 1):
                fp.write('{0:+2.2e}'.format(arr[ii, jj]))
            else:
                fp.write('{0:+2.2e}, '.format(arr[ii, jj]))
        fp.write('\n')
            
def QG_Vortex_Stability():
 
    #For Chebyshev solver

    #Initialize parameters and set up geometry
    paramsCheb   = Parameters()
    GeomCheb     = Geometry("cheb", paramsCheb)
    GeomCheb.Lap = Build_Laplacian(paramsCheb, GeomCheb)
    
    #Set up background-state-flow operators
    GeomCheb.Ψ_op, GeomCheb.Q_op = Build_Bkgd_Operators(paramsCheb, GeomCheb)
    
    kps    = paramsCheb.kps
    kzs    = paramsCheb.kzs
    nmodes = paramsCheb.nmodes

    growthCheb = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    propCheb   = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    modesCheb  = np.zeros([kzs.shape[0], kps.shape[0], nmodes, 
                           paramsCheb.halfNr], 
                          dtype=complex)

    #For finite-difference solver

    #Initialize parameters and set up geometry
    paramsFD   = Parameters()
    GeomFD     = Geometry('FD', paramsFD)
    GeomFD.Lap = Build_Laplacian(paramsFD, GeomFD)

    #Set up background-state-flow operators

    #Array of those r-values at interior gridpoints that lie in physical space
    rInterior = GeomFD.r[1:(paramsFD.halfNr + 1)]

    #Array of (1/r) * (dPsi/dr) evaluated at gridpoints
    Ψ_op_FD = np.ravel(np.sqrt(2*e) * np.exp(-(rInterior**2)))

    #Array of (1/r) * (dQ/dr) evaluated at gridpoints
    Q_op_FD = np.ravel(np.sqrt(32*e) * (rInterior**2 - 2) 
                       * np.exp(-rInterior**2))
    """
    kps    = paramsFD.kps
    kzs    = paramsFD.kzs
    nmodes = paramsFD.nmodes
    """
    growthFD = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    propFD   = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    
    #Start solving

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
 
            recipR2_Cheb = np.diag(
                          np.ravel(1 / GeomCheb.r[1:(paramsCheb.halfNr + 1)]**2)
                                  )

            B_Cheb = (GeomCheb.Lap - (kp2 * recipR2_Cheb) 
                      - (kz2 * (1/paramsCheb.Bu)
                         * np.eye(paramsCheb.halfNr, paramsCheb.halfNr)))
            A_Cheb = np.dot(np.diag(GeomCheb.Ψ_op), B_Cheb) - np.diag(GeomCheb.Q_op)

            #For finite-difference solver

            recipR2_FD = np.diag(
                              np.ravel(1 / GeomFD.r[1:(paramsFD.halfNr + 1)]**2)
                                )
            B_FD = (GeomFD.Lap - (kp2 * recipR2_FD) 
                    - (kz2 * (1/paramsFD.Bu)
                       * np.eye(paramsFD.halfNr, paramsFD.halfNr)))
            A_FD = np.dot(np.diag(Ψ_op_FD), B_FD) - np.diag(Q_op_FD)
            
            ############################
            # FIND EIGENSPACE DIRECTLY #
            ############################
            
            t0 = timeit.timeit()
        
            #Compute eigvals c and eigvecs psi with direct solver ('eig')
            eigValCheb, eigVecCheb = spalg.eig(A_Cheb, B_Cheb)
             
            timeCheb = timeit.timeit() - t0 #Time for direct solve
            
            #Indexing that sorts eigvals by DESCENDING Im(c)
            ind        = (-eigValCheb.imag).argsort()

            eigValCheb = eigValCheb[ind] #Sort eigvals
            eigVecCheb = eigVecCheb[:, ind] #Sort eigvecs in the same order
            omegaCheb  = eigValCheb * kp #Corresp. omegas for this k_phi
            
            growthCheb[kz_idx, kp_idx, :] = omegaCheb[0:nmodes].imag
            propCheb[kz_idx, kp_idx, :]   = omegaCheb[0:nmodes].real
            modesCheb[kz_idx, kp_idx, :, :] = eigVecCheb[0:nmodes]
            
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
    
    dim_growthFD = growthFD * paramsFD.Ro * paramsFD.f0
    dim_propFD   = propFD * paramsFD.Ro * paramsFD.f0

    dim_growthCheb = growthCheb * paramsCheb.Ro * paramsCheb.f0
    dim_propCheb   = propCheb * paramsCheb.Ro * paramsCheb.f0
    
    nkp = (np.ravel(kps)).shape[0]
    nkz = (np.ravel(kzs)).shape[0]

    for jj in range(0, nmodes):
        
        if nkp < 4:
        
            fig, axes = plt.subplots(nkp, 2, figsize = (100, 70), sharex = "col")
            
            for ii in range(0, nkp):
                
                ax_growth = axes[ii, 0]
                ax_growth.plot(kzs, 
                                   4 * np.ravel(dim_growthCheb[:, ii, jj]), '-o')#,
                                   #kzs, 
                                   #4 * np.ravel(dim_growthFD[:, ii, jj]), 
                                   #'-*')
                ax_growth.set(title = 
                              f'Growth rate; $k_{{\phi}}$ = {kps[ii]}', 
                              ylabel = 'Growth rate ($s^{-1}$)')
                
                ax_prop = axes[ii, 1]
                ax_prop.plot(kzs, 
                                 4 * np.ravel(dim_propCheb[:, ii, jj]), 
                                 '-o')#,
                                 #kzs, 
                                 #4 * np.ravel(dim_freqFD[:, ii, jj]), 
                                 #'-*')
                ax_prop.set(title = 
                        f'Propagation speed; $k_{{\phi}}$ = {kps[ii]}',
                            ylabel = 'Azimuthal speed ($s^{-1}$)')

            ax_growth.set(xlabel = r'Vertical wavenumber ($\times \sigma_z$)')
            ax_prop.set(xlabel = r'Vertical wavenumber ($\times \sigma_z$)')

            plt.show()
            fig.savefig(f"omega_vs_m_mode{jj}.png")
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

        #Plot spatial structure of streamfunction

        kz_idx, kp_idx = 0, 1

        eigvec = modesCheb[kz_idx, kp_idx, jj, :] 
        psi    = GetStreamfunc(eigvec, k=2, phi=0)
        psiAmp = np.sqrt(psi.real**2 + psi.imag**2)

        fig, ax = plt.subplots(figsize = (20, 20)) 
         
        ax.scatter(GeomCheb.r[(paramsCheb.halfNr + 1):1:-1], 
                (psi.real)/max(psiAmp), color = "red")
        ax.scatter(GeomCheb.r[(paramsCheb.halfNr + 1):1:-1],
                (psi.imag)/max(psiAmp), color = "indigo")
        #ax.plot(GeomCheb.r[(paramsCheb.halfNr+1):1:-1], amp, color="blue")
        
        ax.set(xlabel = r"$r/\sigma_r$", 
               ylabel = r"Components of $\hat{\psi}$",
              title = f"Real and imaginary components of mode-{jj} eigenvector")
        plt.show()
        fig.savefig(f"eigvec_structure_mode{jj}.png")

if __name__ == '__main__': #For testing
   QG_Vortex_Stability()
