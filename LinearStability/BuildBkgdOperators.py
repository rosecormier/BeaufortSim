import scipy.sparse as sp

from numpy import arange, diag, exp, hstack, meshgrid, ravel, sqrt, stack

def rInterior(params, geom, discretizeAzim = False):
    """
    Build array of those r-values at interior gridpoints that are physical.
    """

    if discretizeAzim:
        halfNr, Nφ     = params.halfNr, params.Nφ
        [rInterior, φ] = meshgrid(geom.r[1:(halfNr + 1)], arange(1, (Nφ + 1)))
        rInterior      = hstack(stack(rInterior[:], axis = -1))

    else:
        halfNr    = params.halfNr
        rInterior = ravel(geom.r[1:(halfNr + 1)])

    return rInterior

def BuildBkgdOperators(params, geom, 
                       discretizeAzimuth = False, 
                       discretizeVertical = False,
                       dimensional_U = 1, dimensional_σr = 1):
    """
    Build array of (1/r) * (dPsi/dr) and array of (1/r) * (dQ/dr), each
     evaluated at gridpoints.
    """

    rInterior = geom.rInterior
    
    if params.bkgd == "GM":
        Ψ_op = ravel(-0.5
                     * exp(-rInterior**2 / dimensional_σr**2)
                     * (dimensional_U / dimensional_σr))
        Q_op = ravel(-2
                     * exp(-rInterior**2 / dimensional_σr**2) 
                     * ((rInterior / dimensional_σr)**2 - 2)
                     * (dimensional_U / dimensional_σr**3))

    elif params.bkgd == "BG":

        if discretizeVertical: #Discretize in both r and z
            Ψ_op = ravel(sqrt(2) 
                         * exp(0.5 - (rInterior**2 / dimensional_σr**2)
                               - (zInterior**2 / dimensional_σz**2))
                         * (dimensional_U / dimensional_σr))
            Q_op = ravel(sqrt(8) 
                         * ((2 * (rInterior / dimensional_σr)**2 - 2)
                            (1 / Bu) * ())
                         * exp(0.5 - (rInterior**2 / dimensional_σr**2)
                               - (zInterior**2 / dimensional_σz**2))
                         * (dimensional_U / dimensional_σr**3))
        
        else: #Discretize in r only
            Ψ_op = ravel(sqrt(2) 
                         * exp(0.5 - (rInterior**2 / dimensional_σr**2))
                         * (dimensional_U / dimensional_σr))
            Q_op = ravel(sqrt(32) 
                         * ((rInterior / dimensional_σr)**2 - 2)
                         * exp(0.5 - rInterior**2 / dimensional_σr**2)
                         * (dimensional_U / dimensional_σr**3))

    if sp.issparse(geom.Dr):
            
        if discretizeAzimuth:
            opLength = params.halfNr * params.Nφ
        else:
            opLength = params.halfNr
    
        Ψ_op = sp.spdiags(Ψ_op, 0, opLength, opLength)
        Q_op = sp.spdiags(Q_op, 0, opLength, opLength)

    else:
        Ψ_op, Q_op = diag(Ψ_op), diag(Q_op)

    return Ψ_op, Q_op
