import scipy.sparse as sp

from numpy import arange, block, diag, exp, eye, hstack, kron, matmul, meshgrid, ravel, sqrt, stack, zeros

from Stratification import N2_profile

def GridInterior(params, geom, discretizeVertical = False):
    """
    Build array of those r- or (r,z)-values at PHYSICAL interior gridpoints.
    If grid is 2D (r, z), also build array of N^2-values at physical interior
     points.
    """

    if discretizeVertical:
        
        halfNr, Nz  = params.halfNr, params.Nz
        N2_function = N2_profile(params.stratification_kw)

        #[rInterior, zInterior]  = meshgrid(geom.r[1:(halfNr + 1)], 
        #                                   geom.z[1:(Nz + 1)])

        halfNz_I, halfNz_Z = eye(Nz // 2), zeros(((Nz // 2), (Nz // 2)))
        rInterior = kron(diag(geom.r[1:(halfNr+1)]), eye(Nz)) + kron(diag(geom.r[1:(halfNr+1)]), block([[halfNz_Z, halfNz_I], [halfNz_I, halfNz_Z]]))
        r2RecipInterior = kron(diag(1/geom.r[1:(halfNr+1)]**2), eye(Nz)) + kron(diag(1/geom.r[1:(halfNr+1)]**2), block([[halfNz_Z, halfNz_I], [halfNz_I, halfNz_Z]]))

        quarterNr_I, quarterNr_Z = eye(halfNr // 2), zeros(((halfNr //2), (halfNr // 2)))
        zInterior = kron(eye(halfNr), diag(geom.z[1:(Nz+1)])) + kron(block([[quarterNr_Z, quarterNr_I], [quarterNr_I, quarterNr_Z]]), diag(geom.z[1:(Nz+1)]))

        #rInterior  = hstack(stack(rInterior[:], axis = -1))
        #zInterior  = hstack(stack(zInterior[:], axis = 0))
        
        N2Recip = diag(1 / ravel(N2_function(geom.z[1:(Nz+1)])))
        N2RecipInterior = kron(eye(halfNr), N2Recip) + kron(block([[quarterNr_Z, quarterNr_I], [quarterNr_I, quarterNr_Z]]), N2Recip)
        
        #Update Dz: take Kron. product of (Nr x Nr) identity matrix and Dz
        geom.Dz = kron(eye(halfNr), geom.Dz[1:(Nz + 1), 1:(Nz + 1)])

        geom.rInterior, geom.zInterior = rInterior, zInterior
        geom.r2Recip = r2RecipInterior
        geom.N2Recip = N2RecipInterior

    else:
       
        halfNr    = params.halfNr
        rInterior = ravel(geom.r[1:(halfNr + 1)])

        return rInterior

def BuildBkgdOperators(params, geom, 
                       discretizeVertical = False,
                       dimensional_U = 1, dimensional_N2 = 1,
                       dimensional_σr = 1, dimensional_σz = 1):
    """
    Build array of (1/r) * (dPsi/dr) and array of (1/r) * (dQ/dr), each
     evaluated at gridpoints.
    """

    rTilde = geom.rInterior / dimensional_σr

    if discretizeVertical:

        Dz = geom.Dz
        Bu = params.Bu
        
        #Nondimensionalize z and N^{-2}
        zTilde       = geom.zInterior / dimensional_σz
        N2RecipTilde = geom.N2Recip * dimensional_N2

    if params.bkgd == "GM":
        Ψ_op = ravel(-0.5 * exp(-rTilde**2)
                     * dimensional_U / dimensional_σr
                    )
        Q_op = ravel(-2 * exp(-rTilde**2) * (rTilde**2 - 2)
                     * dimensional_U / dimensional_σr**3
                    )

    elif params.bkgd == "BG":

        if discretizeVertical: #Discretize in both r and z
            Ψ_op = (sqrt(2) * exp(0.5 - rTilde**2 - zTilde**2)
                         * dimensional_U / dimensional_σr
                   )
            Q_op = (sqrt(8) 
                         * (2 * (rTilde**2 - 2)
                            - (1 / Bu) 
                              * (zTilde * matmul(Dz, N2RecipTilde) 
                                 + (1 - 2 * zTilde**2) * (N2RecipTilde**2)
                                )
                           )
                         * exp(0.5 - rTilde**2 - zTilde**2)
                         * dimensional_U / dimensional_σr**3
                    )
        
        else: #Discretize in r only
            Ψ_op = ravel(sqrt(2) 
                         * exp(0.5 - rTilde**2)
                         * dimensional_U / dimensional_σr
                        )
            Q_op = ravel(sqrt(32) 
                         * ((rTilde**2 - 2)
                            * exp(0.5 - rTilde**2)
                            * dimensional_U / dimensional_σr**3
                           )
                         )
            Ψ_op, Q_op = diag(Ψ_op), diag(Q_op)
    
            #if sp.issparse(geom.Dr):
            #    opLength = params.halfNr
            #    Ψ_op = sp.spdiags(Ψ_op, 0, opLength, opLength)
            #    Q_op = sp.spdiags(Q_op, 0, opLength, opLength)
            #else:
                #Ψ_op, Q_op = diag(Ψ_op), diag(Q_op)

    return Ψ_op, Q_op
