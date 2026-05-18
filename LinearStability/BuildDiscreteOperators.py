import numpy as np
import scipy.sparse as ssp

from math import e, pi

def N2_profile(stratification_kw, dimensional_N2 = 1):
    """
    Return a function to compute N^2 at discrete values of z.
    """

    if stratification_kw == "constant":
        N2_function = lambda z : dimensional_N2 * np.ones_like(z)

    #elif stratification_kw == "baroclinicthermalwind":

    return N2_function

def ComputeRecips(params, geom):
    """
    Compute and save values of 1/r and 1/N^2 at grid points.
    """

    if params.discretizeVertical:

        if params.nondimensional:
            dimensional_N2 = 1
        else:
            dimensional_N2 = params.Nmax**2

        N2_function  = N2_profile(params.stratification_kw, 
                                  dimensional_N2 = dimensional_N2)
        geom.N2      = N2_function(geom.z)
        geom.N2Recip = 1 / geom.N2
        
        geom.rRecip = ssp.diags_array(1 / geom.r[1:(params.halfNr + 1)], format = "csr")
        
    elif not params.discretizeVertical:
        geom.rRecip = np.diag(1 / geom.r[1:(params.halfNr + 1)])

def BuildBkgdOperators(params, geom):
    """
    Build discrete representations of operators Ψ_op := (1/r) * (∂Ψ/∂r) and
     Q_op := (1/r) * (∂Q/∂r) for the prescribed background flow.
    """

    if params.nondimensional:
        dimensional_U, dimensional_σr = 1, 1
    else:
        dimensional_U, dimensional_σr = params.Umax, params.sigmar

    if not params.discretizeVertical:
    
        rTilde = np.ravel(geom.r[1:(params.halfNr + 1)]) / dimensional_σr

        if params.bkgd == "GM":
    
            Ψ_op = -0.5 * np.exp(-rTilde**2) * dimensional_U / dimensional_σr
            Q_op = (-2 * np.exp(-rTilde**2) * (rTilde**2 - 2) * dimensional_U
                    / dimensional_σr**3)
    
            geom.Ψ_op, geom.Q_op = np.diag(Ψ_op), np.diag(Q_op)

        elif params.bkgd == "BG":

            Ψ_op = (np.sqrt(2 * e) * np.exp(-rTilde**2) 
                    * dimensional_U / dimensional_σr)
                    
            Q_op = (np.sqrt(32 * e) * (rTilde**2 - 2) * np.exp(-rTilde**2)
                    * dimensional_U / dimensional_σr**3)
            #This is the barotropic limit, so we are implicitly neglecting a
            # factor of 1/Bu.

            geom.Ψ_op, geom.Q_op = np.diag(Ψ_op), np.diag(Q_op)
                    
    elif params.discretizeVertical:
    
        if params.bkgd == "BG":
        
            #Where to start indexing in z to match size of matrix B
            zStartIdx = int((np.size(geom.z) - params.DzSize) / 2)

            if zStartIdx == 0:
                zEndIdx = np.size(geom.z) + 1
            elif zStartIdx > 0:
                zEndIdx = -zStartIdx
        
            f0, dimensional_σz = params.f0, params.sigmaz
        
            r       = geom.r[1:-1]
            rTilde  = np.ravel(r) / dimensional_σr
            r2Tilde = rTilde**2
            
            z       = geom.z[zStartIdx:zEndIdx]
            zTilde  = np.ravel(z) / dimensional_σz
            z2Tilde = zTilde**2
            
            Dz = geom.Dz[zStartIdx:zEndIdx, zStartIdx:zEndIdx]
            
            #N2Recip = np.ravel(geom.N2Recip[zStartIdx:zEndIdx])
            #N2Recip = ssp.csr_array(geom.N2Recip[zStartIdx:zEndIdx])
            N2Recip = geom.N2Recip[zStartIdx:zEndIdx]

            Ψ_opRadialFactor   = np.diag(np.sqrt(2 * e) * np.exp(-r2Tilde)
                                         * dimensional_U / dimensional_σr)
            Ψ_opVerticalFactor = np.diag(np.exp(-z2Tilde))
            
            Ir = ssp.eye_array(params.Nr - 1, format = "csr") #np.eye(params.Nr - 1)
            Iz = ssp.eye_array(params.DzSize, format = "csr") #np.eye(params.DzSize)

            Ψ_op = (ssp.kron(Ψ_opRadialFactor, Iz, format = "csr") @
                    ssp.kron(Ir, Ψ_opVerticalFactor, format = "csr"))
            #Ψ_op = np.matmul(np.kron(Ψ_opRadialFactor, Iz), np.kron(Ir, Ψ_opVerticalFactor))

            Q_opScaleFactor = (8 * e)**0.5 * (dimensional_U / dimensional_σr**3)
            
            Q_opFactor1 = np.diag(np.exp(-r2Tilde))
            Q_opFactor2 = np.diag(np.exp(-z2Tilde))
            
            #Q_opFactor3RadialTerm   = ssp.diags_array(2 * (rTilde**2 - 2), format = "csr")
            Q_opFactor3RadialTerm = np.diag(2 * (rTilde**2 - 2))
            #Q_opFactor3VerticalTerm = np.diag(-(f0 * dimensional_σr 
            #                                    / dimensional_σz)**2 
            #                                  * (N2Recip * (1 - 2 * z2Tilde) 
            #                                     + z * np.matmul(Dz, N2Recip)
            #                                    )
            #                                 )
            Q_opFactor3VerticalTerm = np.diag(-(f0 * dimensional_σr 
                                                / dimensional_σz)**2 
                                              * (N2Recip * (1 - 2 * z2Tilde) 
                                                 + z * np.matmul(Dz.toarray(), N2Recip)
                                                )
                                             )
            #Q_opFactor3VerticalTerm = ssp.diags_array(
            #                 -(f0 * dimensional_σr / dimensional_σz)**2 
            #                 * (N2Recip * (1 - 2 * z2Tilde) 
            #                    + z * (Dz @ N2Recip)),
            #                    format = "csr")
                                
            Q_opFactor3 = (ssp.kron(Q_opFactor3RadialTerm, Iz, format = "csr")
                           + ssp.kron(Ir, Q_opFactor3VerticalTerm, format = "csr"))
            #Q_opFactor3 = (np.kron(Q_opFactor3RadialTerm, Iz) + np.kron(Ir, Q_opFactor3VerticalTerm))
  
            #Q_op = np.matmul(np.matmul(np.kron(Q_opScaleFactor * Q_opFactor1,
            #                                   Iz),
            #                           np.kron(Ir, Q_opFactor2)),
            #                 Q_opFactor3
            #                )
            Q_op = ((ssp.kron(Q_opScaleFactor * Q_opFactor1, Iz, format = "csr") @
                                       ssp.kron(Ir, Q_opFactor2, format = "csr")) @ Q_opFactor3)
            #N.b., in the limit of very large dimensional_σz, this construction
            # of Q_op agrees with the 1-dimensional Q_op constructed in the 
            # case of the other disjunct (!discretizeVertical).

            geom.Ψ_op, geom.Q_op = Ψ_op, Q_op

def ConvertQuadsToBlock(geom, Q1, Q2, Q3, Q4):
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
    
    if geom.sparse:
        return ssp.block_array([[block1, block3], [block2, block4]], format = "csr")
    elif not geom.sparse:
        return np.block([[block1, block2], [block3, block4]])

def BuildHorizontalLaplacian(params, geom):
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
    geom.LapH_Q1 = geom.Dr2_Q1 + (geom.rRecip @ geom.Dr_Q1)
    geom.LapH_Q2 = geom.Dr2_Q2 + (geom.rRecip @ geom.Dr_Q2)
    geom.LapH_Q3 = geom.Dr2_Q3 + (geom.rRecip @ geom.Dr_Q3)
    geom.LapH_Q4 = geom.Dr2_Q4 + (geom.rRecip @ geom.Dr_Q4)
    #geom.LapH_Q1 = geom.Dr2_Q1 + np.matmul(geom.rRecip, geom.Dr_Q1)
    #geom.LapH_Q2 = geom.Dr2_Q2 + np.matmul(geom.rRecip, geom.Dr_Q2)
    #geom.LapH_Q3 = geom.Dr2_Q3 + np.matmul(geom.rRecip, geom.Dr_Q3)
    #geom.LapH_Q4 = geom.Dr2_Q4 + np.matmul(geom.rRecip, geom.Dr_Q4)
    
    geom.LapH = ConvertQuadsToBlock(geom, geom.LapH_Q1, geom.LapH_Q2, 
                                    geom.LapH_Q3, geom.LapH_Q4)
    
    if params.discretizeVertical:
    
        Iz = ssp.eye_array(params.Nz - 1, format = "csr")
        
        geom.LapH_2D = ssp.kron(geom.LapH, Iz, format = "csr")
    
def BuildMatrixB(params, geom, kφ, kz = None):
    """
    Build discrete representation of 'B' operator in generalized eigenvalue 
     problem.
    """
    
    Nr, halfNr = params.Nr, params.halfNr
    
    kφ2     = kφ**2
    r2Recip = geom.rRecip**2
    
    if not params.discretizeVertical:
    
        kz2 = kz**2
        
        geom.B_Q2 = -geom.LapH_Q2
        geom.B_Q3 = -geom.LapH_Q3
        
        if params.nondimensional:
            geom.B_Q1 = (kφ2 * r2Recip + kz2 * (1 / params.Bu) * np.eye(halfNr)
                         - geom.LapH_Q1)
            geom.B_Q4 = (kφ2 * r2Recip + kz2 * (1 / params.Bu) * np.eye(halfNr)
                         - geom.LapH_Q4)
            
        else:
            geom.B_Q1 = (kφ2 * r2Recip
                         + kz2 * (params.f0 / params.Nmax)**2 * np.eye(halfNr)
                         - geom.LapH_Q1)
            geom.B_Q4 = (kφ2 * r2Recip
                         + kz2 * (params.f0 / params.Nmax)**2 * np.eye(halfNr)
                         - geom.LapH_Q4)

        #Assemble square matrix to be used in gen. eig. solver
        geom.B = geom.B_Q1 + geom.B_Q2
        
        return geom.B
        
    elif params.discretizeVertical:
    
        #Quadrants of terms depending on r, discretized on r-grid
        horizontalB_Q1 = (kφ2 * r2Recip) - geom.LapH_Q1
        horizontalB_Q2 = -geom.LapH_Q2
        horizontalB_Q3 = -geom.LapH_Q3
        horizontalB_Q4 = (kφ2 * r2Recip) - geom.LapH_Q4
        
        #Terms depending on r, discretized on r-grid
        horizontalB = ConvertQuadsToBlock(geom, horizontalB_Q1, horizontalB_Q2,
                                          horizontalB_Q3, horizontalB_Q4)
        
        Iz = ssp.eye(params.DzSize, format = "csr") #np.eye(params.DzSize)
        
        #Deal with vertical boundary conditions
        
        if params.verticalBCs == "homogeneous":
        
            #Terms depending on r, discretized on rz-grid
            horizontalB_2D = ssp.kron(horizontalB, Iz, format = "csr")
            
            #Retain interior values (eigfunction will vanish at boundary pts)
            N2Recip = np.diag(geom.N2Recip[1:-1])
            Dz      = geom.Dz[1:-1, 1:-1]
            Dz2     = geom.Dz2[1:-1, 1:-1]
        
        elif params.verticalBCs == "continuousBuoyancy":
        
            zz = np.zeros((1, (params.Nz + 1)))
            
            #Terms depending on r, discretized on rz-grid
            horizontalB_2D = ssp.kron(horizontalB, Iz, format = "csr") #np.kron(horizontalB, Iz)
            
            #Load these operators in full (w.r.t. z)
            #N2Recip = ssp.diags_array(geom.N2Recip, format = "csr")
            N2Recip = np.diag(geom.N2Recip)
            Dz      = geom.Dz
            Dz2Full = geom.Dz2[1:-1, :].toarray()
            
            #Impose BCs by zeroing first and last rows of the full Dz2
            #Dz2 = ssp.vstack([zz, Dz2Full, zz], format = "csr") 
            Dz2 = np.vstack((zz, Dz2Full, zz))

        #Terms depending on z, discretized on z-grid
        #-(np.matmul((params.f0**2 * N2Recip), Dz2)
        #              + np.matmul(np.matmul(params.f0**2 * Dz, N2Recip), Dz)
        #             )
        #verticalB = -(((params.f0**2 * N2Recip) @ Dz2)
        #              + (((params.f0**2 * Dz) @ N2Recip) @ Dz)
        #             )
        verticalB = -(np.matmul((params.f0**2 * N2Recip), Dz2)
                      + np.matmul(np.matmul(params.f0**2 * Dz.toarray(), N2Recip), Dz.toarray())
                     )
        
        Ir = ssp.eye(params.Nr - 1, format = "csr") #np.eye(params.Nr - 1)
            
        #Terms depending on z, discretized on rz-grid
        verticalB_2D = ssp.kron(Ir, verticalB, format = "csr") #np.kron(Ir, verticalB)
            
        geom.B_Q1 = (horizontalB_2D + verticalB_2D)[:(params.halfNr * params.DzSize),
                                                    :(params.halfNr * params.DzSize)]
        geom.B_Q2 = (horizontalB_2D + verticalB_2D)[(params.halfNr * params.DzSize):,
                                                    :(params.halfNr * params.DzSize)]
        geom.B_Q3 = (horizontalB_2D + verticalB_2D)[:(params.halfNr * params.DzSize),
                                                    (params.halfNr * params.DzSize):]
        geom.B_Q4 = (horizontalB_2D + verticalB_2D)[(params.halfNr * params.DzSize):,
                                                    (params.halfNr * params.DzSize):]

        #Assemble square matrix to be used in gen. eig. solver
        geom.B = geom.B_Q1 + geom.B_Q2

        return geom.B
        
def BuildMatrixA(params, geom):
    """
    Build discrete representation of 'A' operator in generalized eigenvalue 
     problem.
    """

    if not params.discretizeVertical:
        geom.A = np.matmul(geom.Ψ_op, geom.B) + geom.Q_op
    
    elif params.discretizeVertical:
        #(np.matmul(geom.Ψ_op[:(params.halfNr * params.DzSize),
        #                             :(params.halfNr * params.DzSize)],
        #                    geom.B)
        #          + geom.Q_op[:(params.halfNr * params.DzSize),
        #                      :(params.halfNr * params.DzSize)]
        #         )

        geom.A = ((geom.Ψ_op[:(params.halfNr * params.DzSize),
                             :(params.halfNr * params.DzSize)] @ geom.B_Q1 + geom.B_Q2)
                  + geom.Q_op[:(params.halfNr * params.DzSize),
                              :(params.halfNr * params.DzSize)]
                 )

    return geom.A