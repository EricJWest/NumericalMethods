
import numpy as np

def get_brackets(f, a, b, n):
    """
    scans an interval for bracketed roots by performing a simple divide and conquer strategy
    
    INPUT:
    f = function to scan for bracketed roots
    a = left end point of search interval
    b = right end point of search interval
    n = number of subdivisions to consider
    
    OUTPUT:
    list of potential bracketed roots
    """
    
    # create grid
    x = np.linspace(a,b,n)
    
    # initialize output grid, adding dummy row (to be removed later)
    keep = np.array([[a,b]]) 
    
    # test for bracketed roots (sign changes)
    for i in range(0,n-1):
        phi = np.sign(f(x[i]))*np.sign(f(x[i+1]))
        if phi < 0:
            print('sign change detected between',x[i],'and',x[i+1])
            keep = np.concatenate((keep, [[x[i],x[i+1]]]))
    
    # remove first row (dummy row)
    keep = np.delete(keep, 0, axis=0)

    # check whether bracketed roots have been detected
    if keep.size == 0:
        print('WARNING: no sign changes have been detected between',a,'and',b)
    
    return keep