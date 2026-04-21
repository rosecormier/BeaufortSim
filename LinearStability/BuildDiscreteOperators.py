import numpy as np

from math import e, pi

def N2_profile(stratification_kw, dimensional_N2 = 1):

    if stratification_kw == "constant":
        N2_function = lambda z : dimensional_N2 * np.ones_like(z)

    #elif stratification_kw == "baroclinicthermalwind":

    return N2_function

def ComputeRecips(params, geom, discretizeVertical = False,
                  nondimensional = False):

    geom.rRecip = np.diag(1 / geom.r[1:(params.halfNr + 1)])
    
    if discretizeVertical:
    
        if nondimensional:
            dimensional_N2 = 1
        else:
            dimensional_N2 = params.Nmax**2

        N2_function  = N2_profile(params.stratification_kw, 
                                  dimensional_N2 = dimensional_N2)
        geom.N2      = N2_function(geom.z)
        geom.N2Recip = 1 / geom.N2

def BuildBkgdOperators(params, geom, discretizeVertical = False,
                       nondimensional = False):

    if nondimensional:
        dimensional_U, dimensional_σr = 1, 1
    else:
        dimensional_U, dimensional_σr = params.Umax, params.σr

    if not discretizeVertical:
        rTilde = np.ravel(geom.r[1:(params.halfNr + 1)]) / dimensional_σr

    if params.bkgd == "GM":

        Ψ_op = -0.5 * np.exp(-rTilde**2) * dimensional_U / dimensional_σr
        Q_op = (-2 * np.exp(-rTilde**2) * (rTilde**2 - 2) * dimensional_U
                / dimensional_σr**3)

        geom.Ψ_op, geom.Q_op = np.diag(Ψ_op), np.diag(Q_op)

    elif params.bkgd == "BG":

        if not discretizeVertical:
            
            Ψ_op = (np.sqrt(2 * e) * np.exp(-rTilde**2) 
                    * dimensional_U / dimensional_σr)
                    
            Q_op = (np.sqrt(32 * e) * (rTilde**2 - 2) * np.exp(-rTilde**2)
                    * dimensional_U / dimensional_σr**3)
            #This is the barotropic limit, so we are implicitly neglecting a
            # factor of 1/Bu.

            geom.Ψ_op, geom.Q_op = np.diag(Ψ_op), np.diag(Q_op)
                    
        elif discretizeVertical:
        
            f0, dimensional_σz = params.f0, params.σz
        
            r       = geom.r[1:-1]
            rTilde  = np.ravel(r) / dimensional_σr
            r2Tilde = rTilde**2
            
            z       = geom.z[1:-1]
            zTilde  = np.ravel(z) / dimensional_σz
            z2Tilde = zTilde**2
            
            Dz = geom.Dz[1:-1, 1:-1]
            
            N2Recip = np.ravel(geom.N2Recip[1:-1])

            Ψ_opRadialFactor   = np.diag(np.sqrt(2 * e) * np.exp(-r2Tilde)
                                         * dimensional_U / dimensional_σr)
            Ψ_opVerticalFactor = np.diag(np.exp(-z2Tilde))
            
            Ir = np.eye(params.Nr - 1)
            Iz = np.eye(params.Nz - 1)

            Ψ_op = np.matmul(np.kron(Ψ_opRadialFactor, Iz),
                             np.kron(Ir, Ψ_opVerticalFactor))
                             
            Q_opScaleFactor = (8 * e)**0.5 * (dimensional_U / dimensional_σr**3)
            
            Q_opFactor1 = np.diag(np.exp(-r2Tilde))
            Q_opFactor2 = np.diag(np.exp(-z2Tilde))
            
            Q_opFactor3RadialTerm   = np.diag(2 * (rTilde**2 - 2))
            Q_opFactor3VerticalTerm = np.diag(-(f0 * dimensional_σr 
                                                / dimensional_σz)**2 
                                              * (N2Recip * (1 - 2 * z2Tilde) 
                                                 + z * np.matmul(Dz, N2Recip)
                                                )
                                             )
            
            Q_opFactor3 = (np.kron(Q_opFactor3RadialTerm, Iz)
                           + np.kron(Ir, Q_opFactor3VerticalTerm))

            Q_op = np.matmul(np.matmul(np.kron(Q_opScaleFactor * Q_opFactor1,
                                               Iz),
                                       np.kron(Ir, Q_opFactor2)),
                             Q_opFactor3
                            )
            #N.b., in the limit of very large dimensional_σz, this construction
            # of Q_op agrees with the 1-dimensionsal Q_op constructed in the 
            # case of the other disjunct (!discretizeVertical).

            geom.Ψ_op, geom.Q_op = Ψ_op, Q_op

def ConvertQuadsToBlock(Q1, Q2, Q3, Q4):
    """
    Given 4 blocks of a matrix, each indexed from outside to inside of its 
     respective quadrant of the computational domain, re-index and assemble as
     a block matrix, with the result indexed in the global ordering of the 
     computational domain.
    """
    
    block1 = Q1[:, :]
    block2 = Q2[:, ::-1]
    block3 = Q3[::-1, :]
    block4 = Q4[::-1, ::-1]
    
    return np.block([[block1, block2], [block3, block4]])

def BuildHorizontalLaplacian(params, geom):

    halfNr, Nr = params.halfNr, params.Nr

    testfunction         = np.cos(geom.r * pi/(2*params.Lr))
    testfunction1stDeriv = -(pi/(2*params.Lr)) * np.sin(geom.r * pi/(2*params.Lr))
    testfunction2ndDeriv = -(pi/(2*params.Lr))**2 * np.cos(geom.r * pi/(2*params.Lr))
    
    print("Maximum error in discretized first-order r-derivative applied to test function on computational domain:", np.max(np.abs(np.matmul(geom.Dr[1:-1, 1:-1], testfunction[1:-1]) - testfunction1stDeriv[1:-1])))
    
    #Quadrants of 1st-order r-derivative matrix
    Dr_quad1 = geom.Dr[1:(halfNr + 1), 1:(halfNr + 1)]
    Dr_quad2 = geom.Dr[1:(halfNr + 1), (Nr - 1):halfNr:-1]
    Dr_quad3 = geom.Dr[(Nr - 1):halfNr:-1, 1:(halfNr+1)]
    Dr_quad4 = geom.Dr[(Nr - 1):halfNr:-1, (Nr - 1):halfNr:-1]

    print("Maximum error in discretized first-order r-derivative applied to test function on physical domain:", np.max(np.abs(np.matmul(Dr_quad1 + Dr_quad2, testfunction[1:(halfNr+1)]) - testfunction1stDeriv[1:(halfNr+1)])))

    #Quadrants of 2nd-order r-derivative matrix
    Dr2_quad1 = geom.Dr2[1:(halfNr + 1), 1:(halfNr + 1)]
    Dr2_quad2 = geom.Dr2[1:(halfNr + 1), (Nr - 1):halfNr:-1]
    Dr2_quad3 = geom.Dr2[(Nr - 1):halfNr:-1, 1:(halfNr+1)]
    Dr2_quad4 = geom.Dr2[(Nr - 1):halfNr:-1, (Nr - 1):halfNr:-1]
    
    print("Maximum error in discretized second-order r-derivative applied to test function on physical domain:", np.max(np.abs(np.matmul(Dr2_quad1 + Dr2_quad2, testfunction[1:halfNr+1]) - testfunction2ndDeriv[1:halfNr+1])))

    #Quadrants of the full horizontal Laplacian
    geom.Lap_quad1 = (Dr2_quad1 + np.matmul(geom.rRecip, Dr_quad1))
    geom.Lap_quad2 = (Dr2_quad2 + np.matmul(geom.rRecip, Dr_quad2))
    geom.Lap_quad3 = (Dr2_quad3 + np.matmul(geom.rRecip, Dr_quad3))
    geom.Lap_quad4 = (Dr2_quad4 + np.matmul(geom.rRecip, Dr_quad4))

    print("Maximum error in discretized horizontal Laplacian applied to test function on computational domain:", np.max(np.abs(np.matmul(ConvertQuadsToBlock(geom.Lap_quad1, geom.Lap_quad2, geom.Lap_quad3, geom.Lap_quad4), testfunction[1:-1])[:halfNr] - (testfunction2ndDeriv[1:(halfNr+1)] + np.matmul(geom.rRecip, testfunction1stDeriv[1:(halfNr+1)])))))
    
def BuildMatrixB(params, geom, kφ, kz = None, discretizeVertical = False):

    kφ2 = kφ**2
    
    Nr, halfNr, f0 = params.Nr, params.halfNr, params.f0
    
    r2Recip = geom.rRecip**2
    
    if not discretizeVertical:
    
        kz2 = kz**2
    
        N2max = params.Nmax**2
    
        geom.B_quad1 = ((kφ2 * r2Recip) + kz2 * (f0**2 / N2max) * np.eye(halfNr)
                        - geom.Lap_quad1)
        geom.B_quad2 = -geom.Lap_quad2
        geom.B_quad3 = -geom.Lap_quad3
        geom.B_quad4 = ((kφ2 * r2Recip) + kz2 * (f0**2 / N2max) * np.eye(halfNr)
                        - geom.Lap_quad4)

        #Assemble square matrix to be used in gen. eig. solver
        geom.B = geom.B_quad1 + geom.B_quad2
        
        return geom.B

def BuildMatrixA(params, geom, discretizeVertical = False):

    geom.A = np.matmul(geom.Ψ_op, geom.B) + geom.Q_op

    return geom.A