import numpy as np
import scipy.sparse as sp

from math import e

def rInterior(params, geom, discretize2D = False):
    """
    Build array of those r-values at interior gridpoints that are physical.
    """

    if discretize2D:

        halfNr, Np       = params.halfNr, params.Np
        [rInterior, phi] = np.meshgrid(geom.r[1:(halfNr + 1)],
                                       np.arange(1, (Np + 1)))
        rInterior        = np.hstack(np.stack(rInterior[:], axis = -1))

    else:

        halfNr    = params.halfNr
        rInterior = np.ravel(geom.r[1:(halfNr + 1)])

    return rInterior

def BuildBkgdOperators(params, geom, discretize2D = False, 
                       dimensional_U = 1, dimensional_σr = 1):
    """
    Build array of (1/r) * (dPsi/dr) and array of (1/r) * (dQ/dr), each
     evaluated at gridpoints.
    """

    rInterior = geom.rInterior

    if params.bkgd == "GM":
        Ψ_op = np.ravel((-0.5 * np.exp(-rInterior**2 / dimensional_σr**2)) 
                        * (dimensional_σr * dimensional_U))
        Q_op = np.ravel((-2 * np.exp(-rInterior**2 / dimensional_σr**2) 
                         * ((rInterior / dimensional_σr)**2 - 2)) 
                        * (dimensional_U / dimensional_σr))

    elif params.bkgd == "BG":
        Ψ_op = np.ravel((np.sqrt(2*e) * np.exp(-(rInterior**2 / dimensional_σr**2))) 
                        * (dimensional_σr * dimensional_U))
        Q_op = np.ravel((np.sqrt(32*e) * ((rInterior / dimensional_σr)**2 - 2)
                        * np.exp(-rInterior**2 / dimensional_σr**2))
                        * (dimensional_U / dimensional_σr))

    if sp.issparse(geom.Dr):
            
        if discretize2D:
            opLength = params.halfNr * params.Np
        else:
            opLength = params.halfNr
    
        Ψ_op = sp.spdiags(Ψ_op, 0, opLength, opLength)
        Q_op = sp.spdiags(Q_op, 0, opLength, opLength)

    else:
        Ψ_op, Q_op = np.diag(Ψ_op), np.diag(Q_op)

    return Ψ_op, Q_op

