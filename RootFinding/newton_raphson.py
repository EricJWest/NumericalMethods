
import numpy as np
import sys

# Newton-Raphson method
def newton_raphson(f, df, x0, TOL, imax):
    """
    searches for roots of f(x) using the Newton-Raphson method
    
    INPUT:
    f = function whose roots are being sought 
    df = first derivative of function whose roots are being sought
    x0 = initial guess for root
    TOL = allowed tolerance
    imax = maximum number of iterations
    
    OUTPUT:
    location of the root to within the allowed tolerance, or failure message
    """
    
    # initialize iteration
    xOld = x0 # set initial approximant
    print('iteration',0,': initial guess for location of root at',x0)
    
    # iterate search using method of false position
    i = 0  # reset iteration number
    while i <= imax:
    
        # increment iteration number
        i = i + 1
                
        # update approximate location of root
        xNew = xOld - f(xOld)/df(xOld)
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
        
        # rotate approximants
        xOld = xNew
        
    # print message that max iteration has been reached
    print('FAIL! Max number of iterations has been reached. Stopping.')
    return