def mandelbrot(c, max_iter=100):
    """This function determines whether the given c belongs to the Mandelbrot set, 
    i.e. whether the equation z = z ** 2 + c diverges or not when iterated from z = 0.
    
    c diverges if abs(z) ** 2 > 1e10 before max_iter is reached. In this 
    case, the function stops iterating and return the iteration number at which c diverges.
    
    Otherwise, c is bounded, and the function should return max_iter.
    
    Parameters:
    --------------
    c: complex
        A given complex number.
    
    Returns:
    --------------
    int
        i: the number of iterations until c diverges OR
        max_iter: if c does not diverge
    """
    # Initialize z to be 0
    z = 0

    # Iterate the equation z = z ** 2 + c
    for i in range(max_iter):
        z = z ** 2 + c
        
        # Check if abs(z) ** 2 > 1e10, i.e. c diverges
        if abs(z) ** 2 > 1e10:
            return(i)
    
    # c does not diverge
    return(max_iter)
