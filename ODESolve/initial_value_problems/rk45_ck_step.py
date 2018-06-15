
import numpy as np

def rk45_ck_step (x, y, f, h):
    """
    Fourth-order Runge-Kutta-Fehlberg method with fifth-order local truncation error.
    RK weights, nodes, matrix elements are those of Cash and Karp, as recommended by Numerical Recipes.
    
    INPUT:
    x = current x position
    y = current y values (1D array)
    f = rhs function handle (takes x and y as arguments)
    h = current stepsize
    
    OUTPUT:
    YY4 = new y value (1D array), based on lower 4th-order step
    LTE = estimated local truncation error
    """
    
    # Runge-Kutta-Fehlberg weights
    b1 = 37/378
    b2 = 0
    b3 = 250/621
    b4 = 125/594
    b5 = 0
    b6 = 512/1771

    e1 = 2825/27648
    e2 = 0
    e3 = 18575/48384
    e4 = 13525/55296
    e5 = 277/14336
    e6 = 1/4

    # Runge-Kutta-Fehlberg nodes
    c1 = 0
    c2 = 1/5
    c3 = 3/10
    c4 = 3/5
    c5 = 1
    c6 = 7/8
    
    # Runge-Kutta-Fehlberg matrix elements
    a21 = 1/5;
    a31 = 3/40;       a32 = 9/40;
    a41 = 3/10;       a42 = -9/10;   a43 = 6/5;
    a51 = -11/54;     a52 = 5/2;     a53 = -70/27;    a54 = 35/27;
    a61 = 1631/55296; a62 = 175/512; a63 = 575/13824; a64 = 44275/110592; a65 = 253/4096;

    # Runge-Kutta-Fehlberg stages    
    k1 = f(x + h*c1, y)
    k2 = f(x + h*c2, y + h*(a21*k1))
    k3 = f(x + h*c3, y + h*(a31*k1 + a32*k2))
    k4 = f(x + h*c4, y + h*(a41*k1 + a42*k2 + a43*k3))
    k5 = f(x + h*c5, y + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4))
    k6 = f(x + h*c6, y + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5))
       
    # advance the solution, calculate local truncation error
    YY4 = y + h*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
    YY5 = y + h*(e1*k1 + e2*k2 + e3*k3 + e4*k4 + e5*k5 + e6*k6)
    LTE = np.abs(YY4 - YY5)
        
    return YY4, LTE