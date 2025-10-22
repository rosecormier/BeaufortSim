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

def EigenvelocityFromEigvec(params, geom, eigvec, k, **kwargs):
    """
    Evaluates QG velocities corresponding to specified (discrete) eigenvector
     and azimuthal wavenumber, as well as any of φ, z, and t, (m and 
     omega, respectively, must be provided for the latter two).
    """

    halfNr = params.halfNr
    r      = geom.r[1:(halfNr + 1)]
    Dr     = geom.Dr[1:(halfNr + 1), 1:(halfNr + 1)]

    m, ω     = kwargs.get("m"), kwargs.get("ω")
    φ, z, t  = kwargs.get("φ"), kwargs.get("z"), kwargs.get("t")

    ur = 1j * k * (eigvec / r)
    uφ = -np.dot(Dr, eigvec)

    if φ is not None:
        ur = ur * (np.cos(k*φ) + 1j * np.sin(k*φ))
        uφ = uφ * (np.cos(k*φ) + 1j * np.sin(k*φ))

    if m is not None:
        ur = ur * (np.cos(m*z) + 1j * np.sin(m*z))
        uφ = uφ * (np.cos(m*z) + 1j * np.sin(m*z))

    if ω is not None:
        ur = ur * (np.cos(ω*t) - 1j * np.sin(ω*t))
        uφ = uφ * (np.cos(ω*t) - 1j * np.sin(ω*t))

    return ur, uφ
