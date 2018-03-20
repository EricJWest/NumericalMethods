
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
    
    # Runge-Kutta coefficients
    k1 = h*f(x,y)
    k2 = h*f(x + h/2, y + (1/2)*k1)
    k3 = h*f(x + h/2, y + (1/2)*k2)
    k4 = h*f(x + h, y + k3)
        
    # advance the solution
    YY = y + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        
    return YY 