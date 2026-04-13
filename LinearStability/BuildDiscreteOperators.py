import numpy as np

from math import e, pi

def N2_profile(stratification_kw, dimensional_N2 = 1):

    if stratification_kw == "constant":
        N2_function = lambda z : dimensional_N2 * np.ones_like(z)

    #elif stratification_kw == "baroclinicthermalwind":

    return N2_function

def DiscretizeGrid(params, geom, discretizeVertical = False):
    geom.rRecip = np.diag(1 / geom.r[1:params.halfNr+1])

def BuildBkgdOperators(params, geom, discretizeVertical = False,
                       nondimensional = False):

    if nondimensional:
        dimensional_U, dimensional_σr = 1, 1
    else:
        dimensional_U, dimensional_σr = params.Umax, params.σr

    rTilde = np.ravel(geom.r[1:params.halfNr+1]) / dimensional_σr

    if params.bkgd == "GM":

        Ψ_op = -0.5 * np.exp(-rTilde**2) * dimensional_U / dimensional_σr
        Q_op = (-2 * np.exp(-rTilde**2) * (rTilde**2 - 2)
                        * dimensional_U / dimensional_σr**3)

        geom.Ψ_op, geom.Q_op = np.diag(Ψ_op), np.diag(Q_op)

    elif params.bkgd == "BG":
    
        Ψ_opRadialFactor = (np.sqrt(2 * e) * np.exp(-rTilde**2) 
                            * dimensional_U / dimensional_σr)

        if not discretizeVertical:
            
            Q_op = -6*(2*e)**0.5 * (dimensional_U / dimensional_σr**3) * np.exp(-rTilde**2)
            ###(np.sqrt(32 * e) * (dimensional_U / dimensional_σr**3) * np.exp(-rTilde**2) * (rTilde**2 - 2))
                    
            geom.Ψ_op, geom.Q_op = Ψ_opRadialFactor, Q_op
                    
        elif discretizeVertical:
            
            f0, dimensional_σz = params.f0, params.σz

            r2Tilde = rTilde**2
            z       = np.diag(geom.z)
            zTilde  = np.diag(geom.z / dimensional_σz)
            z2Tilde = zTilde**2
            N2Recip = np.diag(geom.N2Recip)
            
            Ψ_opVerticalFactor = np.exp(-z2Tilde)
            
            Ir = np.eye(params.Nr + 1)
            Iz = np.eye(params.Nz + 1)
            
            Q_op = (-np.sqrt(8 * e) * (dimensional_U / dimensional_σr)
                    * np.kron(np.exp(-r2Tilde), Iz)
                    * np.kron(Ir, np.exp(-z2Tilde))
                    * (2 * (2 - np.kron(r2Tilde, Iz)) / dimensional_σr**2)
                       + (f0**2 / dimensional_σz**2)
                          * np.kron(Ir, 
                                    N2Recip - 2 * np.matmul(z2Tilde, N2Recip) 
                                    + np.matmul(z, np.matmul(geom.Dz, N2Recip))
                                   )
                   )
            
            Ψ_op = (np.kron(Ψ_opRadialFactor, Iz) 
                    + np.kron(Ir, Ψ_opVerticalFactor))

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