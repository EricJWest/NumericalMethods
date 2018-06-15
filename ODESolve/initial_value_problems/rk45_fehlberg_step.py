
import numpy as np

def rk45_fehlberg_step (x, y, f, h):
    """
    Fourth-order Runge-Kutta-Fehlberg method with fifth-order local truncation error.
    RK weights, nodes, matrix elements are those of Fehlberg's original paper.
    
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
    b1 = 25/216
    b2 = 0
    b3 = 1408/2565
    b4 = 2197/4104
    b5 = -1/5
    b6 = 0

    e1 = 16/135
    e2 = 0
    e3 = 6656/12825
    e4 = 28561/56430
    e5 = -9/50
    e6 = 2/55

    # Runge-Kutta-Fehlberg nodes
    c1 = 0
    c2 = 1/4
    c3 = 3/8
    c4 = 12/13
    c5 = 1
    c6 = 1/2
    
    # Runge-Kutta-Fehlberg matrix elements
    a21 = 1/4;
    a31 = 3/32;      a32 = 9/32;
    a41 = 1932/2197; a42 = -7200/2197; a43 = 7296/2197;
    a51 = 439/216;   a52 = -8;         a53 = 3680/513;   a54 = -845/4104;
    a61 = -8/27;     a62 = 2;          a63 = -3544/2565; a64 = 1859/4104; a65 = -11/40;

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