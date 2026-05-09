import numpy as np

from math import pi
from numpy import inf
from numpy.linalg import norm

from BuildDiscreteOperators import *
from Chebyshev import Parameters, ChebyshevGeometry

def ErrorsInDiscreteHorizontalLaplacians():
    
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
        
        print("Max. fractional error in 1st-order r-derivative applied to test function on computational 1D domain:",
              norm(((test1stDerivCompDomain - test1stDerivExact[1:-1]) 
                    / test1stDerivExact[1:-1]
                   ),
                   ord = inf
                  ),
              "\nL2 fractional error in 1st-order r-derivative applied to test function on computational 1D domain:",
              norm((test1stDerivCompDomain - test1stDerivExact[1:-1])
                    / test1stDerivExact[1:-1]
                  ),
              "\n"
             )
             
        BuildHorizontalLaplacian(params, geom)
        
        test1stDerivPhysDomain = np.matmul(geom.Dr_Q1 + geom.Dr_Q2,
                                           testFunction[1:(params.halfNr + 1)])
             
        print("Max. fractional error in 1st-order r-derivative applied to test function on physical 1D domain:", 
              norm(((test1stDerivPhysDomain
                     - test1stDerivExact[1:(params.halfNr + 1)])
                    / test1stDerivExact[1:(params.halfNr + 1)]
                   ),
                   ord = inf
                  ),
              "\nL2 fractional error in 1st-order r-derivative applied to test function on computational 1D domain:",
              norm((test1stDerivPhysDomain
                    - test1stDerivExact[1:(params.halfNr + 1)])
                   / test1stDerivExact[1:(params.halfNr + 1)]
                  ),
              "\n"
             )
             
        test2ndDerivPhysDomain = np.matmul(geom.Dr2_Q1 + geom.Dr2_Q2,
                                           testFunction[1:(params.halfNr + 1)])
             
        print("Max. fractional error in 2nd-order r-derivative applied to test function on physical 1D domain:", 
              norm(((test2ndDerivPhysDomain
                     - test2ndDerivExact[1:(params.halfNr + 1)])
                    / test2ndDerivExact[1:(params.halfNr + 1)]
                   ),
                   ord = inf
                  ),
              "\nL2 fractional error in 2nd-order r-derivative applied to test function on physical 1D domain:",
              norm((test2ndDerivPhysDomain
                    - test2ndDerivExact[1:(params.halfNr + 1)])
                   / test2ndDerivExact[1:(params.halfNr + 1)]
                  ),
              "\n"
             )
        
        testLapH = np.matmul(geom.LapH, testFunction[1:-1])[:params.halfNr]
        
        print("Max. fractional error in horizontal Laplacian applied to test function on physical 1D domain:", 
              norm(((testLapH - testHorizontalLapExact) 
                    / testHorizontalLapExact),
                   ord = inf
                  ),
              "\nL2 fractional error in horizontal Laplacian applied to test function on physical 1D domain:",
              norm((testLapH - testHorizontalLapExact) 
                   / testHorizontalLapExact
                  ),
              "\n"
             )
             
        testArgs.update({"Nz": 20, "Lz": 1, "strat_shape": "constant",
                         "sigmaz": None})
             
        params2D = Parameters(args = testArgs, discretizeVertical = True)
        geom2D   = ChebyshevGeometry(params2D)
        
        ComputeRecips(params2D, geom2D, discretizeVertical = True)
        BuildHorizontalLaplacian(params2D, geom2D, discretizeVertical = True)

        iz = np.ones(params2D.Nz - 1)
        
        test2DFunction = np.kron(testFunction[1:-1], iz)

        test2DLapH = np.matmul(geom2D.LapH_2D, 
                               test2DFunction)[:(params2D.halfNr
                                                 * (params2D.Nz - 1))]

        print("Max. fractional error in horizontal Laplacian applied to test function on physical 2D domain:", 
              norm(((test2DLapH - np.kron(testLapH, iz))
                    / np.kron(testLapH, iz)
                   ),
                   ord = inf
                  ),
              "\nL2 fractional error in horizontal Laplacian applied to test function on physical 2D domain:",
              norm((test2DLapH - np.kron(testLapH, iz)) 
                   / np.kron(testLapH, iz)
                  )
             )
             
#def Validate2DBkgdOpsInBarotropicLimit():

if __name__ == '__main__': #For testing
    ErrorsInDiscreteHorizontalLaplacians()