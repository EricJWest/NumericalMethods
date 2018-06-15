
import numpy as np

# back substitution
def back_sub(n, U, y):

    # check that input matrix has the correct dimensions, shape, and size
    if U.ndim != 2:
        print('ERROR: input array must have ndim=2. Stopping.')
        return
    if U.shape[0] != U.shape[1]:
        print('ERROR: input array must have shape=(n,n). Stopping.')
        return
    if U.size != n**2:
        print('ERROR: input array must have size=nxn. Stopping.')
        return

    # check that input matrix is upper-triangular

    
    # check that input vector has the correct dimensions, shape, and size
    if y.ndim != 1:
        print('ERROR: input vector must have ndim=1. Stopping.')
        return
    if y.shape[0] != n:
        print('ERROR: input vector must have shape=(n,). Stopping.')
        return
    if y.size != n:
        print('ERROR: input vector must have size=n. Stopping.')
        return

    # create output array
    x = np.zeros(n)
    
    # back substitution loop
    for i in range (n-1, -1, -1): #loop over rows, in reverse order
        
        # calculate the sum using previously found solutions
        xsum = 0 #reset sum to zero
        for j in range (i+1, n): #loop over columns to the right of the diagonal
            xsum = xsum + U[i,j]*x[j]
        
        # calculate the unknown variable in the current row
        x[i] = (y[i] - xsum)/U[i,i]
        
    return x