import numpy as np

from math import e

def N2_profile(stratification_kw, dimensional_N2 = 1):

    if stratification_kw == "constant":
        N2_function = lambda z : dimensional_N2 * np.ones_like(z)

    #elif stratification_kw == "baroclinicthermalwind":

    return N2_function

def DiscretizeGrid(params, geom, discretizeVertical = False):

    #if not discretizeVertical:
    geom.rRecip  = np.diag(1 / geom.r)
    geom.r2Recip = np.diag(1 / geom.r**2)

    if discretizeVertical:
    
        #Nr, Nz = params.Nr, params.Nz
        
        N2_function  = N2_profile(params.stratification_kw, 
                                 dimensional_N2 = params.Nmax**2)
        geom.N2Recip = (N2_function(geom.z))**(-2)

        #r_2D       = kron(diag(geom.r), eye(Nz + 1))
        #r2Recip_2D = kron(diag(geom.r**(-2)), eye(Nz + 1))

        #z_2D  = kron(eye(Nr + 1), diag(geom.z))

        #N2Recip_2D = kron(eye(Nr + 1), N2Recip)

        #geom.r_2D       = r_2D
        #geom.r2Recip_2D = r2Recip_2D

def BuildFullBkgdOperators(params, geom, discretizeVertical = False,
                           nondimensional = False):

    if nondimensional:
        dimensional_U, dimensional_σr = 1, 1
    else:
        dimensional_U, dimensional_σr = params.Umax, params.σr
        
    if params.bkgd == "GM":
    
        Ψ_op = np.ravel(-0.5 * np.exp(-rTilde**2) * dimensional_U
                        / dimensional_σr)
        Q_op = np.ravel(-2 * np.exp(-rTilde**2) * (rTilde**2 - 2)
                        * dimensional_U / dimensional_σr**3)
                        
        geom.Ψ_op, geom.Q_op = np.diag(Ψ_op), np.diag(Q_op)

    elif params.bkgd == "BG":
    
        rTilde = np.diag(geom.r / dimensional_σr)
    
        Ψ_opRadialFactor = (np.sqrt(2 * e) * np.exp(-rTilde**2) 
                            * dimensional_U / dimensional_σr)

        if not discretizeVertical:
            
            Q_op = (-np.sqrt(8 * e) * (dimensional_U / dimensional_σr)
                    * np.exp(-rTilde**2) 
                    * (2 * (2 - rTilde**2) / dimensional_σr**2))
                    
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

def BuildFullLaplacian(params, geom, discretizeVertical = False):

    halfNr, Nr = params.halfNr, params.Nr

    if not discretizeVertical:

        #Quadrants of 2nd-order r-derivative matrix to retain
        Dr2_quad1 = geom.Dr2[1:(halfNr + 1), 1:(halfNr + 1)] #(pos, pos)
        Dr2_quad2 = geom.Dr2[1:(halfNr + 1), 
                            np.arange(Nr-1, halfNr, -1)]     #(pos, neg)

        #Quadrants of 1st-order r-derivative matrix to retain
        Dr_quad1 = geom.Dr[1:(halfNr + 1), 1:(halfNr + 1)] #(pos, pos)
        Dr_quad2 = geom.Dr[1:(halfNr + 1),
                           np.arange(Nr-1, halfNr, -1)]    #(pos, neg)

        rRecip_quad1 = geom.rRecip[1:(halfNr + 1), 1:(halfNr + 1)]
        rRecip_quad2 = geom.rRecip[1:(halfNr + 1),
                np.arange(Nr-1, halfNr, -1)]
        #Quad 2 is just zeros, but I build it explicitly to be consistent

        geom.Lap_quad1 = (Dr2_quad1 + np.matmul(rRecip_quad1, Dr_quad1))
        geom.Lap_quad2 = (Dr2_quad2 + np.matmul(rRecip_quad2, Dr_quad2))
    
def BuildMatrixB(params, geom, kφ, kz, discretizeVertical = False):

    kφ2, kz2 = kφ**2, kz**2
    
    halfNr, f0, N2max = params.halfNr, params.f0, params.Nmax**2
    
    B_quad1 = (geom.Lap_quad1 
               - (kφ2
                  * geom.r2Recip[1:(halfNr + 1), 1:(halfNr + 1)])
               - kz2 * (f0 / N2max) * np.eye(halfNr))
    B_quad2 = (geom.Lap_quad2 
               - (kφ2
                  * geom.r2Recip[1:(halfNr + 1), (halfNr + 1):-1])
               - kz2 * (f0 / N2max) * np.eye(halfNr))

    geom.B_quad1, geom.B_quad2 = B_quad1, B_quad2
    
    Z = np.zeros((halfNr, halfNr))

    #Assemble and return square matrix
    return (np.block([[B_quad1, Z], [Z, B_quad1[::-1, ::-1]]]) 
            + np.block([[Z, B_quad2], [B_quad2[::-1, ::-1], Z]]))

def BuildMatrixA(params, geom, discretizeVertical = False):

    halfNr = params.halfNr

    A_quad1 = (np.matmul(geom.Ψ_op[1:(halfNr + 1), 1:(halfNr + 1)],
                        geom.B_quad1) 
               - geom.Q_op[1:(halfNr + 1), 1:(halfNr + 1)])
    A_quad2 = (np.matmul(geom.Ψ_op[1:(halfNr + 1), (halfNr + 1):-1], 
                         geom.B_quad1)
               - geom.Q_op[1:(halfNr + 1), (halfNr + 1):-1])
    
    Z = np.zeros((halfNr, halfNr))
    
    #Assemble and return square matrix
    return (np.block([[A_quad1, Z], [Z, A_quad1[::-1, ::-1]]]) 
            + np.block([[Z, A_quad2], [A_quad2[::-1, ::-1], Z]]))