import numpy as np

from math import e, pi

def N2_profile(stratification_kw, dimensional_N2 = 1):
    """
    Return a function to compute N^2 at discrete values of z.
    """

    if stratification_kw == "constant":
        N2_function = lambda z : dimensional_N2 * np.ones_like(z)

    #elif stratification_kw == "baroclinicthermalwind":

    return N2_function

def ComputeRecips(params, geom, discretizeVertical = False,
                  nondimensional = False):
    """
    Compute and save values of 1/r and 1/N^2 at grid points.
    """

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
    """
    Build discrete representations of operators Ψ_op := (1/r) * (∂Ψ/∂r) and
     Q_op := (1/r) * (∂Q/∂r) for the prescribed background flow.
    """

    if nondimensional:
        dimensional_U, dimensional_σr = 1, 1
    else:
        dimensional_U, dimensional_σr = params.Umax, params.sigmar

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
        
            f0, dimensional_σz = params.f0, params.sigmaz
        
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

def BuildHorizontalLaplacian(params, geom, discretizeVertical = False):
    """
    Build discrete representation of horizontal Laplacian in cylindrical
     coordinates, assuming azimuthal symmetry (i.e., neglecting derivatives 
     w.r.t. φ).
    """

    halfNr, Nr = params.halfNr, params.Nr
    
    #Quadrants of 1st-order r-derivative matrix
    geom.Dr_Q1 = geom.Dr[1:(halfNr + 1), 1:(halfNr + 1)]
    geom.Dr_Q2 = geom.Dr[1:(halfNr + 1), (Nr - 1):halfNr:-1]
    geom.Dr_Q3 = geom.Dr[(Nr - 1):halfNr:-1, 1:(halfNr+1)]
    geom.Dr_Q4 = geom.Dr[(Nr - 1):halfNr:-1, (Nr - 1):halfNr:-1]

    #Quadrants of 2nd-order r-derivative matrix
    geom.Dr2_Q1 = geom.Dr2[1:(halfNr + 1), 1:(halfNr + 1)]
    geom.Dr2_Q2 = geom.Dr2[1:(halfNr + 1), (Nr - 1):halfNr:-1]
    geom.Dr2_Q3 = geom.Dr2[(Nr - 1):halfNr:-1, 1:(halfNr+1)]
    geom.Dr2_Q4 = geom.Dr2[(Nr - 1):halfNr:-1, (Nr - 1):halfNr:-1]

    #Quadrants of the full horizontal Laplacian
    geom.LapH_Q1 = geom.Dr2_Q1 + np.matmul(geom.rRecip, geom.Dr_Q1)
    geom.LapH_Q2 = geom.Dr2_Q2 + np.matmul(geom.rRecip, geom.Dr_Q2)
    geom.LapH_Q3 = geom.Dr2_Q3 + np.matmul(geom.rRecip, geom.Dr_Q3)
    geom.LapH_Q4 = geom.Dr2_Q4 + np.matmul(geom.rRecip, geom.Dr_Q4)
    
    geom.LapH = ConvertQuadsToBlock(geom.LapH_Q1, geom.LapH_Q2, 
                                    geom.LapH_Q3, geom.LapH_Q4)
      
    if discretizeVertical:
    
        Iz = np.eye(params.Nz - 1)

        geom.LapH_2D = np.kron(geom.LapH, Iz)
    
def BuildMatrixB(params, geom, kφ, kz = None, discretizeVertical = False):
    """
    Build discrete representation of 'B' operator in generalized eigenvalue 
     problem.
    """
    
    Nr, halfNr, f0 = params.Nr, params.halfNr, params.f0
    
    kφ2     = kφ**2
    r2Recip = geom.rRecip**2
    
    if not discretizeVertical:
    
        kz2 = kz**2
    
        N2max = params.Nmax**2
    
        geom.B_Q1 = ((kφ2 * r2Recip) + kz2 * (f0**2 / N2max) * np.eye(halfNr)
                        - geom.LapH_Q1)
        geom.B_Q2 = -geom.LapH_Q2
        geom.B_Q3 = -geom.LapH_Q3
        geom.B_Q4 = ((kφ2 * r2Recip) + kz2 * (f0**2 / N2max) * np.eye(halfNr)
                        - geom.LapH_Q4)

        #Assemble square matrix to be used in gen. eig. solver
        geom.B = geom.B_Q1 + geom.B_Q2
        
        return geom.B
        
    #elif discretizeVertical:

def BuildMatrixA(params, geom, discretizeVertical = False):
    """
    Build discrete representation of 'A' operator in generalized eigenvalue 
     problem.
    """

    geom.A = np.matmul(geom.Ψ_op, geom.B) + geom.Q_op

    return geom.A