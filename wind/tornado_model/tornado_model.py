import numpy as np

def wind_field(r, z, p):
    """
    Returns the dimensionless components of the tornado velocity (radial,
    circumferential, vertical) as a function of the dimensionless position
    with respect to the core of the tornado, (r,z). Model by Baker & Sterling.
    Parameters: S=swirl ratio, delta=tornado shape factor
    a,b,g refer to alpha, beta and gamma of the model (Default : (1,1,2))
    """
    
    S     = p['S']
    delta = p['delta']
    a     = p['alpha']
    b     = p['beta']
    g     = p['gamma']
    
    if a is None:
        a = 1
    if b is None:
        b = 1
    if g is None:
        g = 2
        
    K = S / (((g - 1) ** ((g - 1) / (a + 1))) * (np.log(2) ** (g / (a + 1))) / (g ** (g / (a + 1))))
    coef = a ** (a / (1 + a)) * b ** (b / (1 + b))

    U = -(1 + a) * (1 + b) / coef * \
        np.power(r, a) * np.power(z, b) / (1 + np.power(r, 1 + a)) / (1 + np.power(z, 1 + b))

    V = K * np.power(r, g - 1) * np.power(np.log(1 + np.power(z, 1 + b)), g / (1 + a)) / \
        np.power(1 + np.power(r, 1 + a), g / (1 + a))

    W = delta * np.power(1 + a, 2) * np.power(r, a - 1) * np.log(1 + np.power(z, 1 + b)) / \
        coef / np.power(1 + np.power(r, 1 + a), 2)

    return U, V, W



def flight_eq(t, y, p):
    """
    Debris flight equation in a tornado (after the model of Baker and Sterling)
    """

    r = y[0]  # radial position
    q = y[1]  # tangential position
    z = y[2]  # height above ground
    u = y[3]  # radial velocity
    v = y[4]  # tangential velocity
    w = y[5]  # vertical velocity

    U, V, W = wind_field(r, z, p)

    R = np.linalg.norm([U - u, V - v, W - w])

    dydt = np.zeros(6)

    dydt[0] = u
    dydt[1] = v
    dydt[2] = w
    dydt[3] = p['phi'] * R * (U - u) + v**2 / r
    dydt[4] = p['phi'] * R * (V - v)
    dydt[5] = p['phi'] * R * (W - w) - p['psi']

    # Set the time derivatives of position to zero if z <= 0
    if z <= 0:
        dydt[0] = u = 0
        dydt[1] = v = 0
        dydt[2] = w = 0

    return dydt
