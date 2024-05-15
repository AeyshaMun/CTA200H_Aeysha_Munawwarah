def lorenz_equations(t, w, sig, r, b):
    """This function defines three Lorenz equations.
    
    Parameters:
    --------------
    t: array-like
        Time parameter (not used as the system is autonomous).
    w: array-like
        Array of three elements represetning the state variables x, y, z.
    sig: array-like
        The Prandtl number (dimensionless).
    r: array-like
        The Rayleigh number (dimensionless).
    b: array-like
        Length scale (dimensionless).
    
    Returns:
    --------------
    dwdt: list
        List of the time derivatives of the state variables.
    """
    # Unpack the state variables
    x, y, z = w
    
    # Define Lorenz equations
    dxdt = -sig * (x-y)
    dydt = r * x - y - x * z
    dzdt = -b * z + x * y
    
    dwdt = [dxdt, dydt, dzdt]
    
    return dwdt
