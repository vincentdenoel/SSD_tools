import numpy as np

def BeamElasticFoundation(EI, k, L, p):
    """
    Solve a beam on elastic foundation problem.
    
    Parameters:
    EI : Flexural rigidity of the beam.
    k  : Stiffness of the elastic foundation.
    L  : Length of the beam.
    p : (array) Vector with applied distributed loads at n points uniformly spaced along the beam.
    
    Returns:
    x : ndarray, Positions along the beam.
    v : ndarray, Displacements at each node.
    M : ndarray, Bending moments at each node.
    """
    n = len(p)   # Number of nodes
    nel = n - 1  # Number of elements

    dx = L / nel
    x = np.linspace(0, L, n)

    K = np.zeros((2 * n, 2 * n))

    for iel in range(nel):
        Ke = np.array([
            [12,       6*dx,    -12,      6*dx],
            [6*dx,     4*dx**2, -6*dx,    2*dx**2],
            [-12,     -6*dx,     12,     -6*dx],
            [6*dx,     2*dx**2, -6*dx,    4*dx**2]
        ]) * (EI / dx**3)

        Ke += np.diag([0.5, 0, 0.5, 0]) * k * dx
        idof = slice(2*iel, 2*iel+4)
        K[idof, idof] += Ke

    P = np.zeros(2 * n)
    P[0::2] = p

    V = np.linalg.solve(K, P)
    v = V[0::2]

    M_left = np.zeros(nel)
    M_right = np.zeros(nel)

    for iel in range(nel):
        idof = slice(2*iel, 2*iel+4)
        vel = V[idof]

        K_left  = np.array([6*dx, 4*dx**2, -6*dx, 2*dx**2]) * (EI / dx**3)
        K_right = np.array([6*dx, 2*dx**2, -6*dx, 4*dx**2]) * (EI / dx**3)

        M_left[iel] = -np.dot(K_left, vel)  # Use "-" for Strength of Materials conventions
        M_right[iel] = np.dot(K_right, vel)

    M = np.zeros(n)
    M[0] = M_left[0]
    M[1:-1] = (M_left[1:] + M_right[:-1]) / 2
    M[-1] = M_right[-1]

    return x, v, M
