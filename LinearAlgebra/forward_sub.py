import numpy as np

# forward substitution
def forward_sub(n, L, b):

    # check that input matrix has the correct dimensions, shape, and size
    if L.ndim != 2:
        print('ERROR: input array must have ndim=2. Stopping.')
        return
    if L.shape[0] != L.shape[1]:
        print('ERROR: input array must have shape=(n,n). Stopping.')
        return
    if L.size != n**2:
        print('ERROR: input array must have size=nxn. Stopping.')
        return

    # check that input matrix is lower-triangular

    
    # check that input vector has the correct dimensions, shape, and size
    if b.ndim != 1:
        print('ERROR: input vector must have ndim=1. Stopping.')
        return
    if b.shape[0] != n:
        print('ERROR: input vector must have shape=(n,). Stopping.')
        return
    if b.size != n:
        print('ERROR: input vector must have size=n. Stopping.')
        return

    # create output array
    y = np.zeros(n)
    
    # forward substitution loop
    for i in range (0, n): #loop over rows
        
        # calculate the sum using previously found solutions
        ysum = 0 #reset sum to zero
        for j in range (0, i): #loop over columns to the left of the diagonal
            ysum = ysum + L[i,j]*y[j]
        
        # calculate the unknown variable in the current row
        y[i] = (b[i] - ysum)/L[i,i]
        
    return y