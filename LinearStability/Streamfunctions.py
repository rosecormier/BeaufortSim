import numpy as np
import scipy.sparse as ssp

from math import cos, sin

from FiniteDiff import FiniteDiff

def Streamfunction(eigvec, **kwargs):
    """
    Evaluates streamfunction at a specified (discrete) eigenvector and any of
     φ, z, and t, where k, m, and omega, respectively, must be provided.
    """

    k, m, ω = kwargs.get("k"), kwargs.get("m"), kwargs.get("ω")
    φ, z, t = kwargs.get("φ"), kwargs.get("z"), kwargs.get("t")

    psi = eigvec

    if k is not None:
        psi = psi * (cos(k*φ) + 1j * sin(k*φ))

    if m is not None:
        psi = psi * (cos(m*z) + 1j * sin(m*z))

    if ω is not None:
        psi = psi * (cos(ω*t) - 1j * sin(ω*t))

    return psi

def EigvelFrom1DEigvec(params, geom, eigMode, kφ_idx, **kwargs):
    """
    Evaluates QG velocities corresponding to eigenvector(s)
     and azimuthal wavenumber, as well as any of φ, z, and t, (m and 
     omega, respectively, must be provided for the latter two).
    """
    
    kφ = params.kps[kφ_idx]

    m, ω    = kwargs.get("m"), kwargs.get("ω")
    φ, z, t = kwargs.get("φ"), kwargs.get("z"), kwargs.get("t")

    if (params.discretizeRadial and not params.discretizeVertical):

        r  = geom.r[0:(params.halfNr + 1)]
        Dr = geom.Dr[1:(params.halfNr + 1), 1:(params.halfNr + 1)]
        
        eigvec = eigMode
        
        #Compute velocity components for single mode 'eigvec'
        
        ur = np.zeros((params.halfNr + 1), dtype = complex)
        uφ = np.zeros((params.halfNr + 1), dtype = complex)
        
        ur[1:] = 1j * kφ * (eigvec[1:] / r[1:])
        uφ[1:] = -np.matmul(Dr, eigvec[1:])
        
        if m is not None:
            ur = ur * (np.cos(m * z) + 1j * np.sin(m * z))
            uφ = uφ * (np.cos(m * z) + 1j * np.sin(m * z))
        
    elif (params.discretizeVertical and not params.discretizeRadial):
        
        r  = params.rs
        Dr = FiniteDiff(r, 1, sparse = False)
        
        ur = np.zeros((len(geom.z), len(r)), dtype = complex)
        uφ = np.zeros((len(geom.z), len(r)), dtype = complex)
        
        for z_idx in range(len(geom.z)):
        
            #Compute eigen-velocities at this z-level
        
            mode_constant_z = eigMode[:, z_idx]
            ur_constant_z   = 1j * kφ * (mode_constant_z / r)
            uφ_constant_z   = -np.matmul(Dr, mode_constant_z)
            
            ur[z_idx, :] = ur_constant_z
            uφ[z_idx, :] = uφ_constant_z

    if φ is not None:
        ur = ur * (np.cos(kφ * φ) + 1j * np.sin(kφ * φ))
        uφ = uφ * (np.cos(kφ * φ) + 1j * np.sin(kφ * φ))

    if ω is not None:
        ur = ur * (np.cos(ω * t) - 1j * np.sin(ω * t))
        uφ = uφ * (np.cos(ω * t) - 1j * np.sin(ω * t))

    return ur, uφ

def EigvelFrom2DEigmode(params, geom, eigMode, k, **kwargs):
    """
    Evaluates QG velocities corresponding to specified (discrete) eigenvector
     and azimuthal wavenumber, as well as either of φ and t, (omega must be
     provided for the latter).
    """

    halfNr = params.halfNr
    Nz     = params.Nz
   
    iz = np.ones(params.Nz + 1) 
    Iz = ssp.eye(params.Nz + 1, format = "csr")

    r     = geom.r[:(halfNr + 1)]
    Dr_2D = ssp.kron(geom.Dr[:(halfNr + 1), :(halfNr + 1)], Iz, format = "csr")
    
    ω, φ, t = kwargs.get("ω"), kwargs.get("φ"), kwargs.get("t")

    ur = 1j * k * (eigMode / np.kron(r, iz))
    uφ = -(Dr_2D @ ssp.csr_array(eigMode)).toarray()

    if φ is not None:
        ur = ur * (np.cos(k * φ) + 1j * np.sin(k * φ))
        uφ = uφ * (np.cos(k * φ) + 1j * np.sin(k * φ))

    if ω is not None:
        ur = ur * (np.cos(ω * t) - 1j * np.sin(ω * t))
        uφ = uφ * (np.cos(ω * t) - 1j * np.sin(ω * t))

    return ur, uφ