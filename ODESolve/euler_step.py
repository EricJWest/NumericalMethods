# simple euler method
def euler_step(x, y, f, h):
    """ 
    Evolve solution by one step using the simple Euler method
    
    INPUT:
    x = current x position
    y = current y values (array)
    f = rhs function handle (takes x and y as arguments)
    h = current stepsize
    
    OUTPUT:
    y_new = new y value
    """
    y_new = y + h*f(x,y)
    return y_new