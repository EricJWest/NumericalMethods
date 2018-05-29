
import numpy as np
import sys

# method of false position
def false_position(f, a, b, TOL, imax):
    """
    searches for roots of f(x) using the method of false position
    
    INPUT:
    f = function whose roots are being sought 
    a = initial left end point
    b = initial right end point
    TOL = allowed tolerance
    imax = maximum number of iterations
    
    OUTPUT:
    location of the root to within the allowed tolerance, or failure message
    """
    
    # test initial bracket 
    if np.sign(f(a))*np.sign(f(b)) > 0:
        # if sgn > 0, problem may not be well-defined
        print('WARNING: function has the same sign at both ends of the interval.')
        print('This may lead to convergence problems. Proceed with caution!')

    # initialize iteration
    xLeft = a   # set left bound
    xRight = b  # set right bound
    pOld = b # set zeroth approximant
    
    # iterate search using method of false position
    i = 0  # reset iteration number
    while i <= imax:
    
        # increment iteration number
        i = i + 1
                
        # update approximate location of root
        pNew = xRight - f(xRight)*(xRight - xLeft)/(f(xRight) - f(xLeft))
        print('iteration',i,': approximate location of root at',pNew)
    
        # calculate errors
        xErr = np.abs((pNew - pOld)/pNew)
        fErr = f(pNew)

        # check if errors are within the allowed tolerance
        if (xErr <= TOL and fErr <= TOL):
            best = pNew #best estimate
            delta = np.abs(pNew-pOld) #uncertainty
            print('SUCCESS! Root has been located to within the specified tolerance after',i,'iterations.')
            print('Root is located at',best,'+/-',delta)
            return

        # print message if max iteration has been reached
        if i == imax:
            print('FAIL! Max number of iterations has been reached. Stopping.')
            return

        # determine which subinterval the root lies in, and update end points      
        sgn = np.sign(f(xLeft))*np.sign(f(pNew))
        if sgn <= 0:
            # if sgn < 0, root is in left subinterval,
            # make midpt the new right end point
            xRight = pNew
            print('root is now brackted by [',xLeft,',',xRight,']')
        elif sgn > 0:
            # if sgn > 0, root is in right subinterval,
            # make midpt the new left end point
            xLeft = pNew 
            print('root is now brackted by [',xLeft,',',xRight,']')
        elif sgn == 0:
            # if sgn = 0, root is located at the midpoint
            print('SUCCESS! Root has been located to within machine precision after',i,'iterations.')
            print('Root is located at',pNew)
            return
        
        # rotate approximants
        pOld = pNew