import numpy as np

from math import cos, sin

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

def EigenvelocityFrom1DEigvec(params, geom, eigvec, k, **kwargs):
    """
    Evaluates QG velocities corresponding to specified (discrete) eigenvector
     and azimuthal wavenumber, as well as any of φ, z, and t, (m and 
     omega, respectively, must be provided for the latter two).
    """

    halfNr = params.halfNr
    r      = geom.r[0:(halfNr + 1)]
    Dr     = geom.Dr[1:(halfNr + 1), 1:(halfNr + 1)] #I need to double-check whether I should be adding quadrants instead

    m, ω     = kwargs.get("m"), kwargs.get("ω")
    φ, z, t  = kwargs.get("φ"), kwargs.get("z"), kwargs.get("t")

    ur = np.zeros((halfNr + 1), dtype = complex)
    uφ = np.zeros((halfNr + 1), dtype = complex)

    ur[1:] = 1j * k * (eigvec[1:] / r[1:])
    uφ[1:] = -np.matmul(Dr, eigvec[1:])

    if φ is not None:
        ur = ur * (np.cos(k * φ) + 1j * np.sin(k * φ))
        uφ = uφ * (np.cos(k * φ) + 1j * np.sin(k * φ))

    if m is not None:
        ur = ur * (np.cos(m * z) + 1j * np.sin(m * z))
        uφ = uφ * (np.cos(m * z) + 1j * np.sin(m * z))

    if ω is not None:
        ur = ur * (np.cos(ω * t) - 1j * np.sin(ω * t))
        uφ = uφ * (np.cos(ω * t) - 1j * np.sin(ω * t))

    return ur, uφ

def EigenvelocityFrom2DEigmode(params, geom, eigmode, k, **kwargs):
    """
    Evaluates QG velocities corresponding to specified (discrete) eigenvector
     and azimuthal wavenumber, as well as either of φ and t, (omega must be
     provided for the latter).
    """

    halfNr = params.halfNr
    Nz     = params.Nz
   
    iz = np.ones(params.Nz + 1) 
    Iz = np.eye(params.Nz + 1)
    
    r     = geom.r[0:(halfNr + 1)]
    Dr_2D = np.kron(geom.Dr[:(halfNr + 1), :(halfNr + 1)], Iz)
    
    ω, φ, t = kwargs.get("ω"), kwargs.get("φ"), kwargs.get("t")
    
    ur = np.zeros((Nz + 1) * (halfNr + 1), dtype = complex)
    uφ = np.zeros((Nz + 1) * (halfNr + 1), dtype = complex)
   
    ur = 1j * k * (eigmode / np.kron(r, iz))
    uφ = -np.matmul(Dr_2D, eigmode)

    if φ is not None:
        ur = ur * (np.cos(k * φ) + 1j * np.sin(k * φ))
        uφ = uφ * (np.cos(k * φ) + 1j * np.sin(k * φ))

    if ω is not None:
        ur = ur * (np.cos(ω * t) - 1j * np.sin(ω * t))
        uφ = uφ * (np.cos(ω * t) - 1j * np.sin(ω * t))

    return ur, uφ