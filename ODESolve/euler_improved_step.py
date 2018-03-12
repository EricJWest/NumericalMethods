# improved euler method
def euler_improved_step(x, y, f, h):
    """ 
    Evolve solution by one step using the improved Euler method
    (a.k.a., second-order Runge-Kutta with a1=1/2, b1=1, c1=1, a2=1/2, b2=1, c2=1)
    
    INPUT:
    x = current x position
    y = current y values (array)
    f = rhs function handle (takes x and y as arguments)
    h = current stepsize
    
    OUTPUT:
    y_new = new y value
    """
    f_i = f(x,y)
    f_i_plus_1 = f(x + h, y + h*f_i)
    y_new = y + 0.5*h*(f_i + f_i_plus_1)
    return y_new