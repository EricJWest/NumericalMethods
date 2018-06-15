# midpt euler method
def euler_midpt_step(x, y, f, h):
    """ 
    Evolve solution by one step using the midpoint Euler method
    (a.k.a., second-order Runge-Kutta with a1=0, a2=1, b2=1/2, c2=1/2)
    
    INPUT:
    x = current x position
    y = current y values (array)
    f = rhs function handle (takes x and y as arguments)
    h = current stepsize
    
    OUTPUT:
    y_new = new y value
    """
    x_mid = x + 0.5*h
    y_mid = y + 0.5*h*f(x,y)
    y_new = y + h*f(x_mid, y_mid)
    return y_new