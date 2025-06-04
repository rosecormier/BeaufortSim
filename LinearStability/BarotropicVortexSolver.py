"""
Modification of Storer's code "Linear Stability of a Barotropic QG Vortex".

Some of the notation follows "Spectral Methods in MATLAB" by L. Trefethen.
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

from scipy.sparse.linalg import eigs
from scipy.interpolate import interp1d
from scipy.special import factorial

from Chebyshev import Chebyshev
from FiniteDiff import FiniteDiff

#############################
# PARSE COMMAND-LINE INPUTS #
#############################

parser = argparse.ArgumentParser()
parser.add_argument('--Neig', 
                    help = 'Number of grid points for eig computations',
                    type = int, default = 1001)
parser.add_argument('--Neigs', 
                    help = 'Number of grid points for eigs computations',
                    type = int, default = 1001)
parser.add_argument('-H', '--depth', 
                    help = 'Fluid depth parameter (DOES NOTHING)',
                    type = float, default = 1e3)
parser.add_argument('-L', '--width', 
                    help = 'Radius of the domain (DOES NOTHING)',
                    type = float, default = 1e6)
parser.add_argument('-f0', '--coriolis', 
                    help = 'Coriolis f0 value (DOES NOTHING)',
                    type = float, default = 1e-4)
parser.add_argument('-g', '--gravity', 
                    help = 'Acceleration due to gravity (DOES NOTHING)',
                    type = float, default = 9.81)
parser.add_argument('-p', '--PrintOutputs', 
                    help = 'Flag to turn on display for each computation',
                    action = 'store_true')
parser.add_argument('-N', '--buoyancy', 
                    help = 'Buoyancy frequency (DOES NOTHING)',
                    type = float, default = np.sqrt(5)*1e-3)
parser.add_argument('-kp', '--k_phi', 
                    help = 'Azimuthal wavenumbers; Enter as -kp min max step',
                    type = float, default = [1, 3, 1], nargs = 3)
parser.add_argument('-kz', '--k_z', 
                    help = 'Vertical wavenumbers; Enter as -kz min max step',
                    type = float, default = [0, 2, 0.1], nargs = 3)
parser.add_argument('--modes', 
                    help = 'Number of modes of instability to be considered',
                    type = int, default = 1)
args = parser.parse_args()
                    
class Parameters:
    
    H        = args.depth
    L        = args.width
    f0       = args.coriolis
    g        = args.gravity
    N        = args.buoyancy

    Lr       = 6.25 #Max. r-value in physical space; half computational domain
    Nr       = args.Neig #Number of gridpoints
    halfNr   = args.Neig // 2
    
    Nt       = 40
    kps      = np.arange(args.k_phi[0], args.k_phi[1], args.k_phi[2])
    kzs      = np.arange(args.k_z[0], args.k_z[1], args.k_z[2])
    
    nmodes   = args.modes
    printout = args.PrintOutputs

    def display(self):
        print('H = {}'.format(self.H))
        print('L = {}'.format(self.L))
        print('f0 = {}'.format(self.f0))
        print('g = {}'.format(self.g))
        print('N = {}'.format(self.N))
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
            
            #Scale gridpoints and variable of differentiation to fit 
            # desired domain
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

    #Quadrants of 2nd-order r-derivative matrix to be retained
    D1 = geom.Dr2[1:halfNr+1, :][:, 1:halfNr+1] #(pos, pos)
    D2 = geom.Dr2[1:halfNr+1, :][:, np.arange(Nr-1, halfNr, -1)] #(pos, neg)

    #Quadrants of 1st-order r-derivatives matrix to be retained
    E1 = geom.Dr[1:halfNr+1, :][:, 1:halfNr+1] #(pos, pos)
    E2 = geom.Dr[1:halfNr+1, :][:, np.arange(Nr-1, halfNr, -1)] #(pos, neg)

    if sp.issparse(geom.Dr):
        R = sp.spdiags(np.transpose(1.0/geom.r[1:halfNr+1]), np.array([0]), halfNr, halfNr)
    else:
        R = np.diag(1.0/np.ravel(geom.r[1:halfNr+1]))

    Lap = D1 + D2 + np.dot(R, E1 + E2)

    return Lap

def Print_npArray(fp, arr):
    for ii in xrange(0, arr.shape[0]):
        for jj in xrange(0, arr.shape[1]):
            if jj == (arr.shape[1] - 1):
                fp.write('{0:+2.2e}'.format(arr[ii, jj]))
            else:
                fp.write('{0:+2.2e}, '.format(arr[ii, jj]))
        fp.write('\n')
            
def QG_Vortex_Stability():

    #Initialize parameters and set up geometry

    paramsCheb = Parameters()
    paramsFD   = Parameters()

    GeomCheb = Geometry('cheb', paramsCheb)
    GeomFD   = Geometry('FD', paramsFD)

    GeomCheb.Lap = Build_Laplacian(paramsCheb, GeomCheb)
    GeomFD.Lap   = Build_Laplacian(paramsFD, GeomFD)

    #Set up background-state flow profile

    rin    = GeomCheb.r[1:(paramsCheb.halfNr + 1)]
    Prsp   = np.ravel(-0.5 * np.exp(-rin**2))            # 1/r*Psi_r
    Qrsp   = np.ravel(-2 * np.exp(-rin**2) * (rin**2 - 2))   # 1/r*Q_r
    
    rin    = GeomFD.r[1:(paramsFD.halfNr + 1)]
    Prfd   = np.ravel(-0.5 * np.exp(-rin**2))            # 1/r*Psi_r
    Qrfd   = np.ravel(-2 * np.exp(-rin**2) * (rin**2 - 2))   # 1/r*Q_r

    kps    = paramsCheb.kps
    kzs    = paramsCheb.kzs
    nmodes = paramsCheb.nmodes
 
    growthsp = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    frequysp = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    growthfd = np.zeros([kzs.shape[0], kps.shape[0], nmodes])
    frequyfd = np.zeros([kzs.shape[0], kps.shape[0], nmodes])

    #Start solving

    for cntz in range(0, kzs.shape[0]):

        kz  = kzs[cntz]
        kz2 = kz**2
  
        for cntp in range(0, kps.shape[0]):

            kp  = kps[cntp]
            kp2 = kp**2
    
            #Build A and B for eigen-analysis

            R2invC = np.diag(
                      np.ravel(1 / GeomCheb.r[1:paramsCheb.halfNr+1]**2))
            Bcheb  = (GeomCheb.Lap - (kp2 * R2invC) 
                      - (kz2 * np.eye(paramsCheb.halfNr, paramsCheb.halfNr)))
            Acheb  = np.dot(np.diag(Prsp), Bcheb) - np.diag(Qrsp)

            R2invF = np.diag(np.ravel(1 / GeomFD.r[1:paramsFD.halfNr+1]**2))
            Bfd    = (GeomFD.Lap - (kp2 * R2invF) 
                      - (kz2 * np.eye(paramsFD.halfNr, paramsFD.halfNr)))
            Afd    = np.dot(np.diag(Prfd), Bfd) - np.diag(Qrfd)
            
            #Find eigenspace directly
            
            t0 = timeit.timeit()

            eigValCheb, eigVecCheb = spalg.eig(Acheb, Bcheb)
            
            timesp = timeit.timeit() - t0
            
            ind        = (-eigValCheb.imag).argsort()
            eigVecCheb = eigVecCheb[:, ind]
            eigValCheb = eigValCheb[ind]

            omegaCheb               = eigValCheb * kp
            growthsp[cntz, cntp, :] = omegaCheb[0:nmodes].imag
            frequysp[cntz, cntp, :] = omegaCheb[0:nmodes].real
            
            for ii in range(0, nmodes): #Loop over modes
            
                grow = omegaCheb[ii].imag
                freq = omegaCheb[ii].real
      
                #Find eigenvalues indirectly

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

                    sig1, vec1 = eigs(np.dot(Afd, Tinv), 1, np.dot(Bfd, Tinv),
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
                
                timefd = timeit.timeit() - t0
                
                omegafd = kp * sig1[0]
                
                growfd = omegafd.imag
                freqfd = omegafd.real

                growthfd[cntz, cntp, ii] = growfd
                frequyfd[cntz, cntp, ii] = freqfd
        
                if paramsCheb.printout: #Display results
                    print('----------')
                    print('kz = {0:4f}, kp = {1:2f}'.format(kz, kp))
                    print('eig : growth rate = {0:+4e}, frequency = {1:+4e}, cputime = {2:+4e}'.format(grow, freq, timesp))
                    print('eigs: growth rate = {0:+4e}, frequency = {1:+4e}, cputime = {2:+4e}'.format(growfd, freqfd, timefd))
                    sys.stdout.flush()

    #Plot eigenvalue results
    
    nkp = (np.ravel(kps)).shape[0]
    nkz = (np.ravel(kzs)).shape[0]
    
    for jj in range(0, nmodes):
        
        plt.figure(jj)
        
        if nkp < 4:
            
            for ii in range(0, nkp):

                plt.subplot(nkp, 2, (1 + 2*ii))
                plt.plot(kzs, 4 * np.ravel(growthfd[:, ii, jj]), '-o',
                         kzs, 4 * np.ravel(growthsp[:, ii, jj]), '-*')
                plt.title('Growth rate for azimuthal wavenumber = {}'.format(ii))
                plt.xlabel('Vertical wavenumber (units?)')
                plt.ylabel('Growth rate (units?)')
        
                plt.subplot(nkp, 2, (2 + 2*ii))
                plt.plot(kzs, 4 * np.ravel(frequyfd[:, ii, jj]), '-o',
                         kzs, 4 * np.ravel(frequysp[:, ii, jj]), '-*')
                plt.title('Propagation speed for azimuthal wavenumber = {}'.format(ii))
                plt.xlabel('Vertical wavenumber (units?)')
                plt.ylabel('Speed (units?)')

        elif nkz < 4:

            for ii in range(0, nkz):
                
                plt.subplot(nkz, 2, (1 + 2*ii))
                plt.plot(np.ravel(kps), 4 * np.ravel(growthfd[ii, :, jj]), '-o',
                         np.ravel(kps), 4 * np.ravel(growthsp[ii, :, jj]), '-*')
                plt.title('Growth rate for vertical wavenumber = {} (units?)'.format(ii))
                plt.xlabel('Azimuthal wavenumber')
                plt.ylabel('Growth rate (units?)')

                plt.subplot(nkz, 2, (2 + 2*ii))
                plt.plot(np.ravel(kps), 4 * np.ravel(frequyfd[ii, :, jj]), '-o',
                         np.ravel(kps), 4 * np.ravel(frequysp[ii, :, jj]), '-*')
                plt.title('Propagation speed for vertical wavenumber = {} (units?)'.format(ii))
                plt.xlabel('Azimuthal wavenumber')
                plt.ylabel('Propagation speed (units?)')

        else:

            plt.subplot(2, 2, 1)
            plt.contour(np.ravel(kps), np.ravel(kzs), 4 * growthfd[:, :, jj])
            plt.title('Growth rate (eigs)')

            plt.subplot(2, 2, 2)
            plt.contour(np.ravel(kps), np.ravel(kzs), 4 * frequyfd[:, :, jj])
            plt.title('Propagation speed (eigs)')

            plt.subplot(2, 2, 3)
            plt.contour(np.ravel(kps), np.ravel(kzs), 4 * growthfd[:, :, jj])
            plt.title('Growth rate (eig)')

            plt.subplot(2, 2, 4)
            plt.contour(np.ravel(kps), np.ravel(kzs), 4 * frequyfd[:, :, jj])
            plt.title('Propagation speed (eig)')

        plt.show()
    
if __name__ == '__main__': #For testing
   QG_Vortex_Stability()
