import numpy as np

def boucWen(y,params):
    """
    Bouc-Wen model, see [1]
    [1] J. Song and A. Der Kiureghian, Generalized Bouc–Wen Model for
    Highly Asymmetric Hysteresis, Journal of Engineering Mechanics,
    June 2006 DOI: 10.1061/(ASCE)0733-9399(2006)132:6(610)
    
    Parameters:
    y (ndarray): Input array containing [u, u_dot, z]. System states.
    u = position, u_dot=velocity, z = hysteretic variable

    params (dict): Dictionary containing model parameters.
    - gamma, beta = parameters that control the shape of the hysteresis loops
    - A, n        = parameters that control the scale and sharpness of the hysteresis loops
    - c0, k0      = linear stiffness and damping
    - a           = parameter that controls post- to pre-yielding stiffness ratio
    
    Returns:
    F (double): Output force.
    dz(double): Nonlinear stiffness.

    Implemented by: V. Denoël, 2024
    """
    u  = y[0]
    du = y[1]
    z  = y[2]

    # original Bouc-Wen model (Wen 1976)
    psi = params['gamma'] + params['beta'] * np.sign(du * z)

    # original Wang-Wen model (1998)
    #psi = psi +  params.phi * ( np.sign(du) + np.sign(z))

    F = params['c0'] * du + params['k0'] * (params['a'] * u + (1 - params['a']) * z)
    dz = (params['A'] - np.abs(z)**params['n'] * psi) * du

    return F, dz