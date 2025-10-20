import numpy as np

from math import cos, sin

def GetStreamfunc(eigvec, **kwargs):
    """
    Evaluates streamfunction at a specified (discrete) eigenvector and any of
     phi, z, and t, where k, m, and omega, respectively, must be provided.
    """

    k, m, omega = kwargs.get("k"), kwargs.get("m"), kwargs.get("omega")
    phi, z, t   = kwargs.get("phi"), kwargs.get("z"), kwargs.get("t")

    psi = eigvec

    if k is not None:
        psi = psi * (cos(k*phi) + 1j * sin(k*phi))

    if m is not None:
        psi = psi * (cos(m*z) + 1j * sin(m*z))

    if omega is not None:
        psi = psi * (cos(omega*t) - 1j * sin(omega*t))

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

    m, omega = kwargs.get("m"), kwargs.get("omega")
    φ, z, t  = kwargs.get("φ"), kwargs.get("z"), kwargs.get("t")

    ur = np.dot((1j * k / r), eigvec)
    uφ = -np.dot(Dr, eigvec)

    if φ is not None:
        ur = ur * (cos(k*φ) + 1j * sin(k*φ))
        uφ = uφ * (cos(k*φ) + 1j * sin(k*φ))

    if m is not None:
        ur = ur * (cos(m*z) + 1j * sin(m*z))
        uφ = uφ * (cos(m*z) + 1j * sin(m*z))

    if omega is not None:
        ur = ur * (cos(omega*t) - 1j * sin(omega*t))
        uφ = uφ * (cos(omega*t) - 1j * sin(omega*t))

    return ur, uφ
