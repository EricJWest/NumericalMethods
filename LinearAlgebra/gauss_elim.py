
import numpy as np

# gaussian elimination
def gauss_elim(n, A, b):
    
    # check if input matrix has the correct dimensions, shape, and size
    if A.ndim != 2:
        print('ERROR: input array must have ndim=2. Stopping.')
        return
    if A.shape[0] != A.shape[1]:
        print('ERROR: input array must have shape=(n,n). Stopping.')
        return
    if A.size != n**2:
        print('ERROR: input array must have size=nxn. Stopping.')
        return

    # check if input vector has the correct dimensions, shape, and size
    if b.ndim != 1:
        print('ERROR: input vector must have ndim=1. Stopping.')
        return
    if b.shape[0] != n:
        print('ERROR: input vector must have shape=(n,). Stopping.')
        return
    if b.size != n:
        print('ERROR: input vector must have size=n. Stopping.')
        return
    
    # calculate the scale factor for each row as the magnitude of
    # the largest element in that row
    s = np.zeros(n) #initialize all scale factors to zero
    for i in range (0, n): #loop over rows
        s[i] = np.amax(A[i,:])
        #print('scale factor for row',i,'is',s[i]) #debugging
        if s[i] == 0:
            print('ERROR: matrix is singular. Stopping.')
            return
    
    # initialize row pointer
    nrow = np.arange(0, n) # nrow i --> i
    
    # create augmented matrix
    b = b.reshape((n,1))
    A = np.concatenate((A,b), axis=1)

    # forward elimination loop
    for j in range (0, n): #loop over columns
        
        #print('forward elimination for column',j,':\n') #debugging
        
        # identify row number of the best pivot candidate
        pjj = np.abs(A[nrow[j],j]/s[nrow[j]]) # relative magnitude of first pivot candidate
        pmax = pjj                            # set first candidate to be the pivot element, so far
        k = j                                 # set first row index to be pivot index, so far
        for i in range (j+1, n): # loop over rows below the diagonal
            pij = np.abs(A[nrow[i],j]/s[nrow[i]]) # relative magnitude of next pivot candidate
            if pij > pmax:
                # if relative magnitude of the current pivot candidate is greater than 
                # the relative magnitude of the previous best pivot candidate...
                pmax = pij  #set the current candidate to be the pivot element, so far
                k = nrow[i] #set the current row index to be the pivot index, so far

        if pmax == 0:
            print('ERROR: matrix is singular. Stopping.')
            return
    
        # simulate row swap by swapping row pointers
        if k != j:
            #print('Pivoting at column',j,': row',k,'<--> row',j) #debugging
            ncopy = nrow[j]   # nrow j --> copy
            nrow[j] = nrow[k] # nrow k --> nrow j
            nrow[k] = ncopy   # copy --> nrow k
        else:
            #print('Pivoting at column',j,': no re-ordering needed') #debugging
            continue
        
        for i in range (j+1, n): #loop over rows below the diagonal        
            # perform row operation on ith row, to eliminate element a_{ij}
            A[nrow[i],:] = A[nrow[i],:] - A[nrow[i],j]/A[nrow[j],j]*A[nrow[j],:]
                
    # create output arrays
    AA = np.zeros((n,n))
    bb = np.zeros(n)
    for i in range (0,n):
        AA[i,:] = A[nrow[i],0:-1]
        bb[i] = A[nrow[i],-1]
    
    return AA, bb