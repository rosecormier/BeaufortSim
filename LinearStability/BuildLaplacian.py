import numpy as np
import scipy.sparse as sp

def BuildLaplacian(params, geom, discretizeVertical = False):
    """
    Build discretized Laplacian operator.
    Note: we do not explicitly set zeroth indices because we impose zero
     Dirichlet BCs.
    """
    
    Nr, halfNr = params.Nr, params.halfNr

    if discretizeVertical:
        Nz   = params.Nz
        I, Z = np.eye(Nz // 2), np.zeros(((Nz // 2), (Nz // 2)))

    #Quadrants of 2nd-order r-derivative matrix to be retained
    D1 = geom.Dr2[1:(halfNr + 1), 1:(halfNr + 1)]              #(pos, pos)
    D2 = geom.Dr2[1:(halfNr + 1), np.arange(Nr-1, halfNr, -1)] #(pos, neg)

    #Quadrants of 1st-order r-derivative matrix to be retained
    #E1 = geom.Dr[1:(halfNr + 1), 1:(halfNr + 1)]              #(pos, pos)
    #E2 = geom.Dr[1:(halfNr + 1), np.arange(Nr-1, halfNr, -1)] #(pos, neg)

    if sp.issparse(geom.Dr):
        
        #Build diagonal matrix from reciprocals of r_j for 1 <= j <= halfNr
        R = sp.spdiags(np.transpose(1 / geom.r[1:(halfNr + 1)]),
                       np.array([0]), halfNr, halfNr)
        
        #Build discretized Laplacian as done in Ch.11 of Trefethen
        if discretizeVertical:
            
            E = sp.kron((R @ geom.Dr[:, :]), sp.eye(Nz))
            
            Lap = sp.kron(D1, sp.eye(Nz)) + E[1:(halfNr + 1), :] + sp.kron(D2, np.block([[Z, I], [I, Z]]))

            #Lap = (sp.kron((D1 + R.dot(E1)), sp.eye(Nz))
            #       + sp.kron((D2 + R.dot(E2)), np.block([[Z, I], [I, Z]])))
        else:
            
            #Quadrants of 1st-order r-derivative matrix to be retained
            E1 = geom.Dr[1:(halfNr + 1), 1:(halfNr + 1)]              #(pos, pos)
            E2 = geom.Dr[1:(halfNr + 1), np.arange(Nr-1, halfNr, -1)] #(pos, neg)
            
            Lap = (D1 + R.dot(E1)) + (D2 + R.dot(E2))

    else:
        
        #Build diagonal matrix from reciprocals of r_j for 1 <= j <= halfNr
        R = np.diag(1 / np.ravel(geom.r[1:(halfNr + 1)]))

        #Build discretized Laplacian as done in Ch.11 of Trefethen
        if discretizeVertical:
            
            E = np.kron(np.dot(R, geom.Dr[1:halfNr+1, :]), np.eye(Nz))

            Lap = np.kron(D1, np.eye(Nz))# + E[:, :] + np.kron(D2, np.block([[Z, I], [I, Z]]))

            #Lap = (np.kron((D1 + np.dot(R, E1)), np.eye(Nz))
            #       + np.kron((D2 + np.dot(R, E2)), np.block([[Z, I], [I, Z]])))
        else:
            Lap = (D1 + np.dot(R, E1)) + (D2 + np.dot(R, E2))
    print("all good")
    return Lap
