
import numpy as np
import sys

# secant method
def secant(f, x0, x1, TOL, imax):
    """
    searches for roots of f(x) using the secant method
    
    INPUT:
    f = function whose roots are being sought 
    x0 = initial guess 1
    x1 = initial guess 2
    TOL = allowed tolerance
    imax = maximum number of iterations
    
    OUTPUT:
    location of the root to within the allowed tolerance, or failure message
    """
    
    # test initial bracket 
    if np.sign(f(x0))*np.sign(f(x1)) > 0:
        # if sgn > 0, problem may not be well-defined
        print('WARNING: function has the same sign at both ends of the interval.')
        print('This may lead to convergence problems. Proceed with caution!')

    # initialize iteration
    xOld = x0
    xNew = x1
    
    # iterate search using bisection method
    i = 1  # reset iteration number
    while i <= imax:
            
        # rotate previous approximate locations of root
        xOldOld = xOld # xOld --> xOldOld
        xOld = xNew    # xNew --> xOld
        
        # update approximate location of root
        fOld = f(xOld)
        xNew = xOld - fOld*(xOld - xOldOld)/(f(xOld) - f(xOldOld))
        print('iteration',i,': approximate location of root at',xNew)
    
        # calculate errors
        xErr = np.abs((xNew - xOld)/xNew)
        fErr = f(xNew)

        # check if errors are within the allowed tolerance
        if (xErr <= TOL and fErr <= TOL):
            best = xNew #best estimate
            delta = np.abs(xNew-xOld) #uncertainty
            print('SUCCESS! Root has been located to within the specified tolerance after',i,'iterations.')
            print('Root is located at',best,'+/-',delta)
            return

        # print message if max iteration has been reached
        if i == imax:
            print('FAIL! Max number of iterations has been reached. Stopping.')
            return
        
        # increment iteration number
        i = i + 1