
import numpy as np

# LU factorization
def lu_factor(n, A):

    # check that input matrix has the correct dimensions, shape, and size
    if A.ndim != 2:
        print('ERROR: input array must have ndim=2. Stopping.')
        return
    if A.shape[0] != A.shape[1]:
        print('ERROR: input array must have shape=(n,n). Stopping.')
        return
    if A.size != n**2:
        print('ERROR: input array must have size=nxn. Stopping.')
        return
    
    # ensure that elements of A are floats
    A = A.astype(float)

    # calculate the scale factor for each row as the magnitude of
    # the largest element in that row
    s = np.zeros(n) #initialize all scale factors to zero
    for i in range (0, n): #loop over rows
        s[i] = np.amax(A[i,:]) #caclulate scale factor
        print('scale factor for row',i,'is',s[i])                      #debugging
        if s[i] == 0:
            print('ERROR: matrix is singular. Stopping.')
            return
    
    # initialize row pointer to keep track of row ordering
    # (pointer indices are swapped during pivoting to simulate row swapping)
    nrow = np.arange(0, n)
    
    # LU factorization loop
    for alpha in range (0, n): #loop over row/column super index
        
        ### final loop ###
        if alpha == n-1:
            
            # final diagonal element
            print('Calculating final element of U:')                      #debugging

            # calculate the sum
            lusum = 0. #reset sum
            for k in range (0, alpha):
                lusum = lusum + A[nrow[i],k]*A[nrow[k],alpha]
        
            # calculate diagonal element
            A[nrow[alpha],alpha] = A[nrow[alpha],alpha] - lusum
            print('  u[',alpha,',',alpha,'] =',A[nrow[alpha],alpha]) #debugging
            
            #then quit
            break
        
        ### calculate diagonal element, with pivoting ###

        # initialize pivoting
        print('Pivoting at column',alpha,':')                         #debugging
        pmax = 0            #max relative magnitude
        prow = alpha        #virtual pivot row number
        nprow = nrow[alpha] #actual pivot row number in stored array

        # calculate pivot candidates in row alpha
        print('  calculating pivot candidates in column',alpha,'...') #debugging
        for i in range (alpha, n): #loop over rows
        
            # calculate the sum using previously found solutions
            lusum = 0. #reset sum
            for k in range (0, alpha):
                lusum = lusum + A[nrow[i],k]*A[nrow[k],alpha]
        
            # calculate the pivot candidate in the current row
            A[nrow[i],alpha] = A[nrow[i],alpha] - lusum
            print('  row',i,': u[',i,',',alpha,'] =',A[nrow[i],alpha]) #debugging

            # check if pivot candidate is the best so far...
            pij = np.abs(A[nrow[i],alpha]/s[nrow[i]]) #calculate relative magnitude
            if pij > pmax:
                # if relative magnitude of the current pivot candidate is greater than 
                # the relative magnitude of the previous best pivot candidate...
                pmax = pij      #current max relative magnitude, so far
                prow = i        #current virtual pivot row number, so far
                nprow = nrow[i] #current actual pivot row number in stored array, so far   
                
        # pivot
        print('  all pivot candidates have been calculated...')      #debugging
        print('  pmax for column',alpha,'is',pmax)                   #debugging
        print('  pivot row for column',alpha,'is',prow)              #debugging
        if pmax == 0:
            print('ERROR: matrix is singular. Stopping.')
            return
    
        # simulate row swap by swapping row pointers
        if prow != nrow[alpha]:
            print('  row swap: row',prow,'<--> row',alpha)           #debugging
            ncopy = nrow[alpha]
            nrow[alpha] = nrow[prow]
            nrow[prow] = ncopy
            print('  new U element, after pivoting:')                #debugging
            print('  u[',alpha,',',alpha,'] =',A[nrow[alpha],alpha]) #debugging 
        else:
            print('  pivoting is not needed, skipping')              #debugging
            continue


        ### calculate column elements of L ###
        print('Calculating elements of L along column',alpha,':')    #debugging
        for i in range (alpha+1, n): #loop over rows
                
            # calculate the element of L in the current row
            A[nrow[i],alpha] = A[nrow[i],alpha]/A[nrow[alpha],alpha]
            print('  l[',i,',',alpha,'] =',A[nrow[i],alpha])         #debugging
            
            
        ### calculate row elements of U ###
        print('Calculating elements of U along row',alpha,':')       #debugging
        for j in range (alpha+1, n): #loop over columns
        
            # calculate the sum using previously found solutions
            lusum = 0. #reset sum
            for k in range (0, alpha):
                lusum = lusum + A[nrow[alpha],k]*A[nrow[k],j]
        
            # calculate the element of U in the current column
            A[nrow[alpha],j] = A[nrow[alpha],j] - lusum
            print('  u[',alpha,',',j,'] =',A[nrow[alpha],j])         #debugging

            
    # create output arrays
    P = np.zeros((n,n)) #initialize P as a zero array
    L = np.eye(n)       #initialize L as diag(1,...,1)
    U = np.zeros((n,n)) #initialize U as a zero array
    for i in range (0, n):
        P[i, nrow[i]] = 1             #fill permutation matrix 
        U[i, i:] = A[nrow[i], i:]     #fill upper-triangular elements of U
    for i in range (1, n):
        L[i, :i] = A[nrow[i], :i]     #fill lower-triangular elements of L

    return P, L, U