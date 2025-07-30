from math import cos, sin

def GetStreamfunc(eigenvec, **kwargs):
    """
    Evaluates streamfunction at a specified (discrete) eigenvector and any of
     phi, z, and t, where k, m, and omega, respectively, must be provided.
    """

    k, m, omega = kwargs.get("k"), kwargs.get("m"), kwargs.get("omega")
    phi, z, t   = kwargs.get("phi"), kwargs.get("z"), kwargs.get("t")

    psi = eigenvec

    if k is not None:
        psi = psi * (cos(k*phi) + 1j * sin(k*phi))

    if m is not None:
        psi = psi * (cos(m*z) + 1j * sin(m*z))

    if omega is not None:
        psi = psi * (cos(omega*t) - 1j * sin(omega*t))

    return psi
