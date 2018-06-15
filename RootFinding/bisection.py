
import numpy as np
import sys

# bisection method
def bisect(f, a, b, TOL, imax):
    """
    searches for roots of f(x) using the bisection method
    
    INPUT:
    f = function whose roots are being sought 
    a = initial left bracket
    b = initial right bracket
    TOL = allowed tolerance
    imax = maximum number of iterations
    
    OUTPUT:
    location of the root to within the allowed tolerance, or failure message
    """
    
    # test initial bracket 
    if np.sign(f(a))*np.sign(f(b)) > 0:
        # if sgn > 0, problem may not be well-defined
        print('ERROR: function has the same sign at both ends of the interval.')
        print('Bisection method is not well-suited to this problem, with the given input.')
        print('Try a different initial interval, or choose a different method. Stopping.')
        return

    # initialize iteration
    xLeft = a   # set left bound
    xRight = b  # set right bound
    xMid = xLeft + (xRight - xLeft)/2. # calculate midpoint
    
    # iterate search using bisection method
    i = 0  # reset iteration number
    while i <= imax:
    
        # increment iteration number
        i = i + 1

        # determine which subinterval the root lies in, and update end points      
        sgn = np.sign(f(xLeft))*np.sign(f(xMid))
        if sgn <= 0:
            # if sgn < 0, root is in left subinterval,
            # make midpt the new right end point
            print('iteration',i,': root is located in left subinterval, between',xLeft,'and',xMid)
            xRight = xMid
        elif sgn > 0:
            # if sgn > 0, root is in right subinterval,
            # make midpt the new left end point
            print('iteration',i,': root is located in right subinterval, between',xMid,'and',xRight)
            xLeft = xMid 
        elif sgn == 0:
            # if sgn = 0, root is located at the midpoint
            print('SUCCESS! Root has been located to within machine precision after',i,'iterations.')
            print('Root is located at',xMid)
            return

        # calculate new midpoint
        xMid = xLeft + (xRight - xLeft)/2.
    
        # calculate errors
        xErr = np.abs((xRight - xLeft)/xMid)
        fErr = f(xMid)

        # check if errors are within the allowed tolerance
        if (xErr <= TOL and fErr <= TOL):
            best = xMid #best estimate
            delta = xRight-xLeft #uncertainty
            print('SUCCESS! Root has been located to within the specified tolerance after',i,'iterations.')
            print('Root is located at',best,'+/-',delta)
            return

    # print message that max iteration has been reached
    print('FAIL! Max number of iterations has been reached. Stopping.')
    return