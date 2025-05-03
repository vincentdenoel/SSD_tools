import numpy as np

def stiffness_12dof(a, b, D, nu):
    """
    Calculates the stiffness matrix for a 12-DOF system.
    
    Parameters:
    a (float): Half-dimension along x-axis.
    b (float): Half-dimension along y-axis.
    D (float): Bending modulus of the plate
    nu (float): Poisson's ratio.
    
    Returns:
    np.ndarray: Symmetric stiffness matrix (12x12).
    """

    # Define the matrix k_tmp
    sym = 1./2.
    k_tmp = np.zeros((12, 12))
    k_tmp[3, 3] = 16 * a * b * D * sym
    k_tmp[3, 5] = 16 * a * b * D * nu
    k_tmp[4, 4] = 8 * a * b * D * (1 - nu) * sym
    k_tmp[4, 10] = 8 * a**3 * b * D * (1 - nu)
    k_tmp[4, 11] = 8 * a * b**3 * D * (1 - nu)
    k_tmp[5, 5] = 16 * a * b * D * sym
    k_tmp[6, 6] = 48 * a**3 * b * D * sym
    k_tmp[6, 8] = 16 * a**3 * b * D * nu
    k_tmp[7, 7] = (16 / 3) * a * b**3 * D + (32 / 3) * a**3 * b * D * (1 - nu) * sym
    k_tmp[7, 9] = 16 * a * b**3 * D * nu
    k_tmp[8, 8] = ((16 / 3) * a**3 * b * D + (32 / 3) * a * b**3 * D * (1 - nu)) * sym
    k_tmp[9, 9] = 48 * a * b**3 * D * sym
    k_tmp[10, 10] = (16 * a**3 * b**3 * D + (72 / 5) * a**5 * b * D * (1 - nu)) * sym
    k_tmp[10, 11] = (4 / 9) * a**3 * b**3 * (18 * D * (1 - nu) + 36 * D * nu)
    k_tmp[11, 11] = (16 * a**3 * b**3 * D + (72 / 5) * a * b**5 * D * (1 - nu)) * sym

    # Ensure symmetry in k_tmp
    k_tmp = (k_tmp + k_tmp.T)

    # Calculate the stiffness matrix Ke
    C = const_12dof(a, b)
    C_inv = np.linalg.inv(C.T)
    Ke = C_inv @ k_tmp @ np.linalg.inv(C)
    Ke = (Ke + Ke.T) / 2  # Symmetrize Ke

    return Ke


def const_12dof(a, b):
    # Define the constant matrix C for a 12-DOF element
    # The matrix C is used to transform the local coordinates to global coordinates

    C = np.array([
        [1, -a, -b, a**2, a*b, b**2, -a**3, -a**2*b, -a*b**2, -b**3, a**3*b, a*b**3],
        [0,  0,  1,    0,  -a, -2*b, 0, a**2, 2*a*b, 3*b**2, -a**3, -3*a*b**2],
        [0,  1,  0, -2*a,  -b,    0, 3*a**2, 2*a*b, b**2, 0, -3*a**2*b, -b**3],
        [1, -a,  b, a**2,-a*b, b**2, -a**3, a**2*b, -a*b**2, b**3, -a**3*b, -a*b**3],
        [0,  0,  1,    0,  -a,  2*b, 0, a**2, -2*a*b, 3*b**2, -a**3, -3*a*b**2],
        [0,  1,  0, -2*a,   b,    0, 3*a**2, -2*a*b, b**2, 0, 3*a**2*b, b**3],
        [1,  a,  b, a**2, a*b, b**2, a**3, a**2*b, a*b**2, b**3, a**3*b, a*b**3],
        [0,  0,  1,    0,   a,  2*b, 0, a**2, 2*a*b, 3*b**2, a**3, 3*a*b**2],
        [0,  1,  0,  2*a,   b,    0, 3*a**2, 2*a*b, b**2, 0, 3*a**2*b, b**3],
        [1,  a, -b, a**2,-a*b, b**2, a**3, -a**2*b, a*b**2, -b**3, -a**3*b, -a*b**3],
        [0,  0,  1,    0,   a, -2*b, 0, a**2, -2*a*b, 3*b**2, a**3, 3*a*b**2],
        [0,  1,  0,  2*a,  -b,    0, 3*a**2, -2*a*b, b**2, 0, -3*a**2*b, -b**3]
    ])
    return C

def mass_12dof(a, b, mu):
    """
    Calculate the mass matrix for a 12-DOF element.

    Parameters:
        a (float): Half-dimension along x.
        b (float): Half-dimension along y.
        mu (float): Mass per unit surface.

    Returns:
        np.ndarray: A 12x12 mass matrix.
    """
    
    # Define the matrix m_tmp 
    m_tmp = mu * np.array([
        [4*a*b, 0, 0, (4*a**3*b)/3, 0, (4*a*b**3)/3, 0, 0, 0, 0, 0, 0],
        [0, (4*a**3*b)/3, 0, 0, 0, 0, (4*a**5*b)/5, 0, (4*a**3*b**3)/9, 0, 0, 0],
        [0, 0, (4*a*b**3)/3, 0, 0, 0, 0, (4*a**3*b**3)/9, 0, (4*a*b**5)/5, 0, 0],
        [(4*a**3*b)/3, 0, 0, (4*a**5*b)/5, 0, (4*a**3*b**3)/9, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, (4*a**3*b**3)/9, 0, 0, 0, 0, 0, (4*a**5*b**3)/15, (4*a**3*b**5)/15],
        [(4*a*b**3)/3, 0, 0, (4*a**3*b**3)/9, 0, (4*a*b**5)/5, 0, 0, 0, 0, 0, 0],
        [0, (4*a**5*b)/5, 0, 0, 0, 0, (4*a**7*b)/7, 0, (4*a**5*b**3)/15, 0, 0, 0],
        [0, 0, (4*a**3*b**3)/9, 0, 0, 0, 0, (4*a**5*b**3)/15, 0, (4*a**3*b**5)/15, 0, 0],
        [0, (4*a**3*b**3)/9, 0, 0, 0, 0, (4*a**5*b**3)/15, 0, (4*a**3*b**5)/15, 0, 0, 0],
        [0, 0, (4*a*b**5)/5, 0, 0, 0, 0, (4*a**3*b**5)/15, 0, (4*a*b**7)/7, 0, 0],
        [0, 0, 0, 0, (4*a**5*b**3)/15, 0, 0, 0, 0, 0, (4*a**7*b**3)/21, (4*a**5*b**5)/25],
        [0, 0, 0, 0, (4*a**3*b**5)/15, 0, 0, 0, 0, 0, (4*a**5*b**5)/25, (4*a**3*b**7)/21]
    ])
    
    C = const_12dof(a, b)

    # Compute the mass matrix Me
    Me = np.linalg.solve(C.T, m_tmp @ np.linalg.inv(C))

    # Symmetrize the result
    Me = (Me + Me.T) / 2

    return Me


def load_12dof(a, b, q):
    """
    Calculate the load vector for a 12-DOF element.

    Parameters:
        a (float): Dimension of the element (e.g., length).
        b (float): Dimension of the element (e.g., width).
        q (float): Load intensity.

    Returns:
        np.ndarray: A 12x1 load vector.
    """
    # Calculate the scaling factor
    fac = (1 / 3) * q * a * b
    
    # Define the load vector
    pe = fac * np.array([3, b, -a, 3, -b, -a, 3, b, a, 3, -b, a])
    
    return pe

def moment_coefs(a, b, nu, D):
    """
    Calculate the coefficients for the moment matrix. Returns results for 4 nodes.
    Strenght of material conventions.
    Returns:
        Coefficients to get moments from the nodal displacements:
        (Mx1, My1, Mxy1,  Mx2, My2, Mxy2, Mx3, My3, Mxy3, Mx4, My4, Mxy4) = Kmoment_coefs(a,b) * u_el
    """
    
    # Derivatives (D2x, D2y, D2xy) at each node
    matrix1 = [
        [0, 0, 0, 2, 0, 0, -6*a, -2*b, 0, 0, 6*a*b, 0],
        [0, 0, 0, 0, 0, 2, 0, 0, -2*a, -6*b, 0, 6*a*b],
        [0, 0, 0, 0, 1, 0, 0, -2*a, -2*b, 0, 3*a**2, 3*b**2]
    ]

    matrix2 = [
        [0, 0, 0, 2, 0, 0, -6*a, 2*b, 0, 0, -6*a*b, 0],
        [0, 0, 0, 0, 0, 2, 0, 0, -2*a, 6*b, 0, -6*a*b],
        [0, 0, 0, 0, 1, 0, 0, -2*a, 2*b, 0, 3*a**2, 3*b**2]
    ]

    matrix3 = [
        [0, 0, 0, 2, 0, 0, 6*a, 2*b, 0, 0, 6*a*b, 0],
        [0, 0, 0, 0, 0, 2, 0, 0, 2*a, 6*b, 0, 6*a*b],
        [0, 0, 0, 0, 1, 0, 0, 2*a, 2*b, 0, 3*a**2, 3*b**2]
    ]

    matrix4 = [
        [0, 0, 0, 2, 0, 0, 6*a, -2*b, 0, 0, -6*a*b, 0],
        [0, 0, 0, 0, 0, 2, 0, 0, 2*a, -6*b, 0, -6*a*b],
        [0, 0, 0, 0, 1, 0, 0, 2*a, -2*b, 0, 3*a**2, 3*b**2]
    ]
    
    derivatives = [matrix1, matrix2, matrix3, matrix4]

    C = const_12dof(a, b)
    C_inv = np.linalg.inv(C)

    coefMoments = []
    for matrix in derivatives:
        coefM = np.zeros((3, 12))
        coefM[0, :] = -D * (np.array(matrix[1]) + nu * np.array(matrix[0]))
        coefM[1, :] = -D * (np.array(matrix[0]) + nu * np.array(matrix[1]))
        coefM[2, :] = D * np.array(matrix[2]) * (1-nu)
        coefMoments.append(coefM @ C_inv)

    return coefMoments