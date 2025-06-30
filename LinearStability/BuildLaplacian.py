import numpy as np
import scipy.sparse as sp

def BuildLaplacian(params, geom):
    """
    Build discretized Laplacian operator.
    Note: we do not explicitly set zeroth indices because we impose zero
     Dirichlet BCs.
    """

    Nr, halfNr, Np = params.Nr, params.halfNr, params.Np

    I, Z = np.eye(Np // 2), np.zeros(((Np // 2), (Np // 2)))

    #Quadrants of 2nd-order r-derivative matrix to be retained
    D1 = geom.Dr2[1:(halfNr + 1), 1:(halfNr + 1)] #(pos, pos)
    D2 = geom.Dr2[1:(halfNr + 1), np.arange(Nr-1, halfNr, -1)] #(pos, neg)

    #Quadrants of 1st-order r-derivative matrix to be retained
    E1 = geom.Dr[1:(halfNr + 1), 1:(halfNr + 1)] #(pos, pos)
    E2 = geom.Dr[1:(halfNr + 1), np.arange(Nr-1, halfNr, -1)] #(pos, neg)

    #Build diagonal matrix from reciprocals of r_j for 1 <= j <= halfNr
    if sp.issparse(geom.Dr):
        R = sp.spdiags(np.transpose(1 / geom.r[1:(halfNr + 1)]),
                       np.array([0]), halfNr, halfNr)
    else:
        R = np.diag(1 / np.ravel(geom.r[1:(halfNr + 1)]))

    #Build discretized Laplacian as done in Ch.11 of Trefethen
    Lap = (np.kron((D1 + np.dot(R, E1)), np.eye(Np))
           + np.kron((D2 + np.dot(R, E2)), np.block([[Z, I], [I, Z]])))

    return Lap
