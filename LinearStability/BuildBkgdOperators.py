import numpy as np
import scipy.sparse as sp

def BuildBkgdOperators(params, geom):
    """
    Build array of (1/r) * (dPsi/dr) and array of (1/r) * (dQ/dr), each
     evaluated at gridpoints.
    """

    halfNr, Np       = params.halfNr, params.Np
    [rInterior, phi] = np.meshgrid(geom.r[1:(halfNr + 1)],
                                   np.arange(1, (Np + 1)))

    #Array of those r-values at interior gridpoints that lie in physical space
    rInterior = np.hstack(np.stack(rInterior[:], axis = -1))

    if params.bkgd == "GM":
        Ψ_op = np.ravel(-0.5 * np.exp(-rInterior**2))
        Q_op = np.ravel(-2 * np.exp(-rInterior**2) * (rInterior**2 - 2))

        if sp.issparse(geom.Dr):
            Ψ_op = sp.spdiags(Ψ_op, 0, (halfNr * Np), (halfNr * Np))
            Q_op = sp.spdiags(Q_op, 0, (halfNr * Np), (halfNr * Np))

        else:
            Ψ_op, Q_op = np.diag(Ψ_op), np.diag(Q_op)

    elif params.bkgd == "BG":
        Ψ_op = np.ravel(np.sqrt(2*e) * np.exp(-(rInterior**2)))
        Q_op = np.ravel(np.sqrt(32*e) * (rInterior**2 - 2)
                        * np.exp(-rInterior**2))

    return rInterior, Ψ_op, Q_op

