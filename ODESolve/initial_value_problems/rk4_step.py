
import numpy as np

def rk4_step (x, y, f, h):
    """
    Fourth-order Runge-Kutta method to solve a system of first-order 
    differential equations subject to initial conditions.
    
    INPUT:
    x = current x position
    y = current y values (1D array)
    f = rhs function handle (takes x and y as arguments)
    h = current stepsize
    
    OUTPUT:
    YY = new y value (1D array)
    """
        
    # Runge-Kutta weights
    b1 = 1/6
    b2 = 1/3
    b3 = 1/3
    b4 = 1/6
    
    # Runge-Kutta nodes
    c1 = 0
    c2 = 1/2
    c3 = 1/2
    c4 = 1
    
    # Runge-Kutta matrix elements
    a21 = 1/2;
    a31 = 0; a32 = 1/2 ;
    a41 = 0; a42 = 0 ; a43 = 1;

    # Runge-Kutta stages
    k1 = f(x + h*c1, y)
    k2 = f(x + h*c2, y + h*(a21*k1))
    k3 = f(x + h*c3, y + h*(a31*k1 + a32*k2))
    k4 = f(x + h*c4, y + h*(a41*k1 + a42*k2 + a43*k3))
       
    # advance the solution
    YY = y + h*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
        
    return YY 