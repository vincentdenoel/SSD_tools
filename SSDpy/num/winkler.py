import numpy as np

def BeamElasticFoundation(EI, k, L, p, x_hinge=None):
    """
    Solve a beam on elastic foundation problem.
    
    Parameters:
    EI : Flexural rigidity of the beam.
    k  : Stiffness of the elastic foundation.
    L  : Length of the beam.
    p : (array) Vector with applied distributed loads at n points uniformly spaced along the beam.
    x_hinge : (1D array of integers) contains indices related to the position of hinges if any, otherwise None.
              Numbering refers to the definition of the load vector p.
              The first node is 0, the second node is 1, etc. The last node is n-1.
    
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

    if x_hinge is None:
        for iel in range(nel):
            Ke = np.array([
                [12,       6*dx,    -12,      6*dx],
                [6*dx,     4*dx**2, -6*dx,    2*dx**2],
                [-12,     -6*dx,     12,     -6*dx],
                [6*dx,     2*dx**2, -6*dx,    4*dx**2]
            ]) * (EI / dx**3)

            #Ke += np.diag([0.5, 0, 0.5, 0]) * k * dx
            Ke += np.diag([0.5, dx/12, 0.5, -dx/12]) * k * dx

            idof = slice(2*iel, 2*iel+4)
            K[idof, idof] += Ke
    else:
        x_hinge =np.sort(x_hinge)
        # remove values smaller than 1 and larger than n-2
        x_hinge = x_hinge[x_hinge > 0]
        x_hinge = x_hinge[x_hinge < n-1]

        for iel in range(nel):
            if iel+1 not in x_hinge:
                Ke = np.array([
                    [12,       6*dx,    -12,      6*dx],
                    [6*dx,     4*dx**2, -6*dx,    2*dx**2],
                    [-12,     -6*dx,     12,     -6*dx],
                    [6*dx,     2*dx**2, -6*dx,    4*dx**2]
                ]) * (EI / dx**3)
            else:
                Ke = np.array([
                    [3 ,       3*dx,     -3,     0],
                    [3*dx,     3*dx**2, -3*dx,   0],
                    [-3,      -3*dx,      3,     0],
                    [0,           0,      0,     0]
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

        if x_hinge is None or (iel + 1) not in x_hinge:

            K_left  = np.array([6*dx, 4*dx**2, -6*dx, 2*dx**2]) * (EI / dx**3)
            K_right = np.array([6*dx, 2*dx**2, -6*dx, 4*dx**2]) * (EI / dx**3)
        else:
            K_left  = np.array([3, 3*dx, -3, 0]) * (EI / dx**3)
            K_right = np.array([3*dx, 3*dx**2, -3*dx, 0]) * (EI / dx**3)

        M_left[iel] = -np.dot(K_left, vel)  # Use "-" for Strength of Materials conventions
        M_right[iel] = np.dot(K_right, vel)

    M = np.zeros(n)
    M[0] = M_left[0]
    M[1:-1] = (M_left[1:] + M_right[:-1]) / 2
    M[-1] = M_right[-1]

    return x, v, M
