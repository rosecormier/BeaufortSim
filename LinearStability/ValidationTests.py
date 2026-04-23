import numpy as np

from math import pi
from numpy import inf

from Chebyshev import Parameters, ChebyshevGeometry

def ConvertQuadsToBlock(Q1, Q2, Q3, Q4):
    """
    Helper function.
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


def Validate1DHorizontalLaplacian():

    from numpy.linalg import norm

    from BuildDiscreteOperators import BuildHorizontalLaplacian, ComputeRecips
    
    testArgs = {"Nr": 201, "Lr": 1, "bkgd": "BG",
                "buoyancyfreq": 1, "Coriolis": None, "bkgdU": None, 
                "sigmar": None, "Np": None, "k_phi": [0, 1, 1], 
                "k_z": [0, 1, 1], "nmodes": None}
        
    params = Parameters(args = testArgs)
    geom   = ChebyshevGeometry(params)

    cosineTestFunction = lambda r : np.cos(r * pi / (2 * params.Lr))

    cosineTest1stDeriv = lambda r : (-(pi / (2 * params.Lr))
                                     * np.sin(r * pi / (2 * params.Lr)))
    
    cosineTest2ndDeriv = lambda r : (-(pi / (2 * params.Lr))**2
                                     * np.cos(r * pi / (2 * params.Lr)))
    
    testFunctions = [cosineTestFunction]
    
    test1stDerivs = {cosineTestFunction: cosineTest1stDeriv}
    
    test2ndDerivs = {cosineTestFunction: cosineTest2ndDeriv}
    
    for testFunctionExpression in testFunctions:
    
        test1stDerivExpression = test1stDerivs[testFunctionExpression]
        test2ndDerivExpression = test2ndDerivs[testFunctionExpression]
      
        ComputeRecips(params, geom)
      
        testFunction      = testFunctionExpression(geom.r)
        test1stDerivExact = test1stDerivExpression(geom.r)
        test2ndDerivExact = test2ndDerivExpression(geom.r)
        
        testHorizontalLapExact = (test2ndDerivExact[1:(params.halfNr + 1)]
                                  + np.matmul(geom.rRecip,
                                              test1stDerivExact[1:(params.halfNr
                                                                   + 1)]
                                             )
                                 )
        
        test1stDerivCompDomain = np.matmul(geom.Dr[1:-1, 1:-1],
                                           testFunction[1:-1])
        
        print("Max. error in 1st-order r-derivative applied to test function on computational 1D domain:",
              norm(test1stDerivCompDomain - test1stDerivExact[1:-1], ord = inf
                  ),
              "\nL2-error in 1st-order r-derivative applied to test function on computational 1D domain:",
              norm(test1stDerivCompDomain - test1stDerivExact[1:-1]),
              "\n"
             )
             
        BuildHorizontalLaplacian(params, geom)
        
        test1stDerivPhysDomain = np.matmul(geom.Dr_Q1 + geom.Dr_Q2,
                                           testFunction[1:(params.halfNr + 1)])
             
        print("Max. error in 1st-order r-derivative applied to test function on physical 1D domain:", 
              norm(test1stDerivPhysDomain
                   - test1stDerivExact[1:(params.halfNr + 1)], ord = inf
                  ),
              "\nL2-error in 1st-order r-derivative applied to test function on computational 1D domain:",
              norm(test1stDerivPhysDomain
                   - test1stDerivExact[1:(params.halfNr + 1)]
                  ),
              "\n"
             )
             
        test2ndDerivPhysDomain = np.matmul(geom.Dr2_Q1 + geom.Dr2_Q2,
                                           testFunction[1:(params.halfNr + 1)])
             
        print("Max. error in 2nd-order r-derivative applied to test function on physical 1D domain:", 
              norm(test2ndDerivPhysDomain
                   - test2ndDerivExact[1:(params.halfNr + 1)], ord = inf
                  ),
              "\nL2-error in 2nd-order r-derivative applied to test function on physical 1D domain:",
              norm(test2ndDerivPhysDomain
                   - test2ndDerivExact[1:(params.halfNr + 1)]
                  ),
              "\n"
             )
        
        blockLap = ConvertQuadsToBlock(geom.LapH_Q1, geom.LapH_Q2,
                                       geom.LapH_Q3, geom.LapH_Q4)
        testLap  = np.matmul(blockLap, testFunction[1:-1])[:params.halfNr]
        
        print("Max. error in horizontal Laplacian applied to test function on physical 1D domain:", 
              norm(testLap - testHorizontalLapExact, ord = inf),
              "\nL2-error in horizontal Laplacian applied to test function on physical 1D domain:",
              norm(testLap - testHorizontalLapExact)
             )
             
        testArgs.update({"Nz": 20, "Lz": 1, "strat_shape": "constant",
                         "sigmaz": None})
             
        params2D = Parameters(args = testArgs, discretizeVertical = True)
        geom2D   = ChebyshevGeometry(params2D)
        
        ComputeRecips(params2D, geom2D, discretizeVertical = True)
        BuildHorizontalLaplacian(params2D, geom2D, discretizeVertical = True)

        Iz = np.eye(params2D.Nz - 1)
        
        test2DFunction = np.kron(np.diag(testFunction[1:-1]), Iz)
        
        test2DLap = np.matmul(geom2D.LapH, 
                              test2DFunction)[:(params2D.halfNr
                                                * (params2D.Nz - 1)),
                                              :(params2D.halfNr
                                                * (params2D.Nz - 1))]
                                              
        print(test2DLap - np.kron((np.diag(test2ndDerivExact[1:(params2D.halfNr+1)]) + np.matmul(geom2D.rRecip, np.diag(test1stDerivExact[1:(params2D.halfNr+1)]))), Iz))

#def Validate2DBkgdOpsInBarotropicLimit():

if __name__ == '__main__': #For testing
    Validate1DHorizontalLaplacian()