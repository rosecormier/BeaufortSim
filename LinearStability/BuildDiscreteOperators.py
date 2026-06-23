import numpy as np
import scipy.sparse as ssp

from math import e, pi

def N2_profile(params, dimensional_N2_far = 1):
    """
    Return a function to compute N^2 at discrete values of z.
    """

    def TWB_N2(r, z):

        f0, dimensional_U              = params.f0, params.Umax
        dimensional_σr, dimensional_σz = params.sigmar, params.sigmaz
        
        N2 = ((np.sqrt(2) * dimensional_σr * f0 * dimensional_U 
               / (dimensional_σz**2))
              * (1 - np.exp(-(r / dimensional_σr)**2)) 
              * (1 - 2 * z / (dimensional_σz**2)) 
              * np.exp(0.5 - (z / dimensional_σz)**2)
             )
        return N2

    def fromDoubleTanhN2(z, doubleTanhParams):
    
        g, rho_0, A_s, z_s, C_s, A_d, z_d, C_d = params.doubleTanhParams
        
        N2 = -(g / rho_0) * ((A_s / C_s) 
                             * (1 / (np.cosh((z - z_s) / C_s))**2) 
                             + (A_d / C_d) 
                             * (1 / (np.cosh((z - z_d) / C_d))**2)
                            )
        return N2
       
    def constantN2(z):
        return dimensional_N2_far * np.ones_like(z)

    if params.stratification_kw == "constant":
        totalN2function = lambda r, z : constantN2(z)
        print("Warning: to ensure thermal-wind balance of background state, ensure U is set to 0 (for no background flow) or sigma_z is set very large (to approximate a barotropic background). \n")
        
    elif params.stratification_kw == "TWB":
        totalN2function = lambda r, z : constantN2(z) + TWB_N2(r, z)

    elif params.stratification_kw == "doubleTanh":
        totalN2function = lambda r, z : (constantN2(z) 
                                 + fromDoubleTanhN2(z, params.doubleTanhParams))
        print("Warning: to ensure thermal-wind balance of background state, ensure U is set to 0 (for no background flow) or sigma_z is set very large (to approximate a barotropic background). \n")
        
    elif params.stratification_kw == "doubleTanhTWB":
        totalN2function = lambda r, z : (constantN2(z) 
                  + TWB_N2(r, z) + fromDoubleTanhN2(z, params.doubleTanhParams))

    return totalN2function

def ComputeRecips(params, geom, r = None):
    """
    Compute and save values of 1/r and 1/N^2 at grid points.
    """

    if params.discretizeVertical:

        if params.nondimensional:
            dimensional_N2 = 1
        else:
            dimensional_N2 = params.N2_far

        N2_function  = N2_profile(params, 
                                  dimensional_N2_far = dimensional_N2)
        geom.N2      = N2_function(r, geom.z)
        geom.N2Recip = 1 / geom.N2
        
        #2D eigenvalue problem
        if params.discretizeRadial:
            geom.rRecip = ssp.diags_array(1 / geom.r[1:(params.halfNr + 1)],
                                          format = "csr")
        
    #1D (r) eigenvalue problem
    elif (params.discretizeRadial and not params.discretizeVertical):
        geom.rRecip = np.diag(1 / geom.r[1:(params.halfNr + 1)])

    #1D (z) eigenvalue problem - TESTING
    if (params.discretizeVertical and not params.discretizeRadial):
        geom.rRecip = 1 / r

def BuildBkgdOperators(params, geom, r_idx = None):
    """
    Build discrete representations of operators Ψ_op := (1/r) * (∂Ψ/∂r) and
     Q_op := (1/r) * (∂Q/∂r) for the prescribed background flow.
    """

    if params.nondimensional:
        dimensional_U, dimensional_σr = 1, 1
    else:
        dimensional_U, dimensional_σr = params.Umax, params.sigmar

    #1D (r) eigenvalue problem
    if (params.discretizeRadial and not params.discretizeVertical):
    
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
        
            #2D eigenvalue problem
            if params.discretizeRadial:
                rTilde = np.ravel(geom.r[1:-1]) / dimensional_σr
                
            #1D (z) eigenvalue problem
            elif not params.discretizeRadial:
                rTilde = np.array([params.rs[r_idx] / dimensional_σr])
                
            r2Tilde = rTilde**2
                
            z       = geom.z[zStartIdx:zEndIdx]
            zTilde  = np.ravel(z) / dimensional_σz
            z2Tilde = zTilde**2
            Dz      = geom.Dz[zStartIdx:zEndIdx, zStartIdx:zEndIdx]
            N2Recip = geom.N2Recip[zStartIdx:zEndIdx]

            Ψ_opRadialFactor   = np.diag(np.sqrt(2 * e) * np.exp(-r2Tilde)
                                         * dimensional_U / dimensional_σr)
            Ψ_opVerticalFactor = np.diag(np.exp(-z2Tilde))
            
            Q_opScaleFactor = (8 * e)**0.5 * (dimensional_U / dimensional_σr**3)
            
            Q_opFactor1 = np.diag(np.exp(-r2Tilde))
            Q_opFactor2 = np.diag(np.exp(-z2Tilde))

            Q_opFactor3RadialTerm = np.diag(2 * (rTilde**2 - 2))
            
            #2D eigenvalue problem
            if params.discretizeRadial:
            
                Ir = ssp.eye_array(params.Nr - 1, format = "csr")
                Iz = ssp.eye_array(params.DzSize, format = "csr")

                Ψ_op = (ssp.kron(Ψ_opRadialFactor, Iz, format = "csr") @
                        ssp.kron(Ir, Ψ_opVerticalFactor, format = "csr"))

            
                Q_opFactor3VerticalTerm = np.diag(-(f0 * dimensional_σr 
                                                    / dimensional_σz)**2 
                                                  * (N2Recip * (1 - 2 * z2Tilde)
                                                     + z 
                                                       * np.matmul(Dz.toarray(),
                                                                   N2Recip)
                                                    )
                                                 )
  
                Q_opFactor3 = (ssp.kron(Q_opFactor3RadialTerm, Iz, 
                                        format = "csr")
                               + ssp.kron(Ir, Q_opFactor3VerticalTerm, 
                                          format = "csr")
                              )

                Q_op = ((ssp.kron(Q_opScaleFactor * Q_opFactor1, Iz, 
                                  format = "csr") 
                         @ ssp.kron(Ir, Q_opFactor2, format = "csr")
                        ) 
                        @ Q_opFactor3
                       )
                    
            #1D (z) eigenvalue problem   
            elif not params.discretizeRadial:
                
                Ψ_op = Ψ_opRadialFactor * Ψ_opVerticalFactor
                
                Q_opFactor3VerticalTerm = np.diag(-(f0 * dimensional_σr 
                                                    / dimensional_σz)**2 
                                                  * (N2Recip * (1 - 2 * z2Tilde)
                                                     + z 
                                                       * np.matmul(Dz, N2Recip)
                                                    )
                                                 )
                
                Q_op = (Q_opScaleFactor * Q_opFactor1 
                        * np.matmul(Q_opFactor2, 
                                    (Q_opFactor3RadialTerm 
                                     + Q_opFactor3VerticalTerm)
                                   )
                       )

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
    
    geom.LapH = ConvertQuadsToBlock(geom, geom.LapH_Q1, geom.LapH_Q2, 
                                    geom.LapH_Q3, geom.LapH_Q4)
    
    #2D eigenvalue problem
    if params.discretizeVertical:
        Iz           = ssp.eye_array(params.Nz - 1, format = "csr")
        geom.LapH_2D = ssp.kron(geom.LapH, Iz, format = "csr")
    
def BuildMatrixB(params, geom, kφ, kz = None, r_idx = None):
    """
    Build discrete representation of 'B' operator in generalized eigenvalue 
     problem.
    """
    
    kφ2 = kφ**2
    
    if params.discretizeRadial:
        Nr, halfNr = params.Nr, params.halfNr
        r2Recip    = geom.rRecip**2

    #1D (r) eigenvalue problem
    if (params.discretizeRadial and not params.discretizeVertical):
    
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
        
    #1D (z) eigenvalue problem
    elif (params.discretizeVertical and not params.discretizeRadial):
    
        #Deal with vertical boundary conditions
        
        if params.verticalBCs == "homogeneous":
        
            #Retain interior values (eigfunction will vanish at boundary pts)
            N2Recip = np.diag(geom.N2Recip[1:-1])
            Dz      = geom.Dz[1:-1, 1:-1]
            Dz2     = geom.Dz2[1:-1, 1:-1]
            
        elif params.verticalBCs == "constantStreamfunction":
            
            zz = np.zeros((1, (params.Nz + 1)))

            N2Recip = np.diag(geom.N2Recip)
            DzFull  = geom.Dz[1:-1, :]

            #Impose BCs by zeroing first and last rows of the full Dz
            Dz  = np.vstack((zz, DzFull, zz))
            Dz2 = geom.Dz2
            
        elif params.verticalBCs == "continuousBuoyancy":
            
            zz = np.zeros((1, (params.Nz + 1)))

            N2Recip = np.diag(geom.N2Recip)
            DzFull  = geom.Dz
            Dz2Full = geom.Dz2[1:-1, :]
            
            #Impose BCs by zeroing first and last rows of the full Dz2
            Dz2 = np.vstack((zz, Dz2Full, zz))

        geom.B = -(np.matmul((params.f0**2 * N2Recip), Dz2)
                   + np.matmul(np.matmul(params.f0**2 * Dz, N2Recip), Dz)
                  ) + kφ2 * geom.rRecip**2
    
        return geom.B
        
    #2D eigenvalue problem
    elif (params.discretizeRadial and params.discretizeVertical):
    
        #Quadrants of terms depending on r, discretized on r-grid
        horizontalB_Q1 = (kφ2 * r2Recip) - geom.LapH_Q1
        horizontalB_Q2 = -geom.LapH_Q2
        horizontalB_Q3 = -geom.LapH_Q3
        horizontalB_Q4 = (kφ2 * r2Recip) - geom.LapH_Q4
        
        #Terms depending on r, discretized on r-grid
        horizontalB = ConvertQuadsToBlock(geom, horizontalB_Q1, horizontalB_Q2,
                                          horizontalB_Q3, horizontalB_Q4)
        
        Iz = ssp.eye_array(params.DzSize, format = "csr")
        
        #Deal with vertical boundary conditions
        
        if params.verticalBCs == "homogeneous":
        
            #Terms depending on r, discretized on rz-grid
            horizontalB_2D = ssp.kron(horizontalB, Iz, format = "csr")
            
            #Retain interior values (eigfunction will vanish at boundary pts)
            N2Recip = np.diag(geom.N2Recip[1:-1])
            Dz      = geom.Dz[1:-1, 1:-1].toarray()
            Dz2     = geom.Dz2[1:-1, 1:-1].toarray()
        
        elif params.verticalBCs == "continuousBuoyancy":
        
            zz = np.zeros((1, (params.Nz + 1)))
            
            #Terms depending on r, discretized on rz-grid
            horizontalB_2D = ssp.kron(horizontalB, Iz, format = "csr")
            
            #Load these operators in full (w.r.t. z)
            N2Recip = np.diag(geom.N2Recip)
            Dz      = geom.Dz.toarray()
            Dz2Full = geom.Dz2[1:-1, :].toarray()
            
            #Impose BCs by zeroing first and last rows of the full Dz2
            Dz2 = np.vstack((zz, Dz2Full, zz))

        #Terms depending on z, discretized on z-grid
        verticalB = -(np.matmul((params.f0**2 * N2Recip), Dz2)
                      + np.matmul(np.matmul(params.f0**2 * Dz, N2Recip), Dz)
                     )
        
        Ir = ssp.eye_array(params.Nr - 1, format = "csr")
            
        #Terms depending on z, discretized on rz-grid
        verticalB_2D = ssp.kron(Ir, verticalB, format = "csr")
            
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

    #Any 1D eigenvalue problem
    if ((params.discretizeRadial and not params.discretizeVertical) or
        (params.discretizeVertical and not params.discretizeRadial)):
        geom.A = np.matmul(geom.Ψ_op, geom.B) + geom.Q_op  
    
    #2D eigenvalue problem
    elif (params.discretizeRadial and params.discretizeVertical):
        geom.A = ((geom.Ψ_op[:(params.halfNr * params.DzSize),
                             :(params.halfNr * params.DzSize)] @ geom.B_Q1 
                   + geom.B_Q2)
                  + geom.Q_op[:(params.halfNr * params.DzSize),
                              :(params.halfNr * params.DzSize)]
                 )

    return geom.A