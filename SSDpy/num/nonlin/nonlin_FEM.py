import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import SSDpy as ssd


def initialize(xnod, ynod, elemnoa, elemnob):
    """
    Initialize x0 and l0 based on the input nodal coordinates and element connectivity.

    Parameters:
    xnod (list or numpy array): x-coordinates of nodes.
    ynod (list or numpy array): y-coordinates of nodes.
    elemnoa (list or numpy array): Starting node indices of elements.
    elemnob (list or numpy array): Ending node indices of elements.

    Returns:
    x0 (list of numpy arrays): List of displacement vectors for each element.
    l0 (list of floats): List of lengths of each element.
    """
    nelem = len(elemnoa)
    x0 = []  # To store displacement vectors
    l0 = []  # To store lengths

    for i in range(nelem):
        dx = xnod[elemnob[i]-1] - xnod[elemnoa[i]-1]
        dy = ynod[elemnob[i]-1] - ynod[elemnoa[i]-1]
        displacement = np.array([dx, dy])
        x0.append(displacement)
        l0.append(np.linalg.norm(displacement))
    
    return x0, l0

def fhe(u, ndof, elemnoa, elemnob, EA, x0, l0, listfix, k_inf, lambda_, p):
    """
    Computes the residual vector at all degrees of freedom.

    Parameters:
    u (numpy.ndarray): Displacement vector.
    ndof (int): Number of degrees of freedom.
    elemnoa (list): Starting node indices for each element.
    elemnob (list): Ending node indices for each element.
    EA (list): Axial stiffness of elements.
    x0 (list): Initial position vectors of elements.
    l0 (list): Lengths of elements.
    listfix (list): Indices of fixed degrees of freedom.
    k_inf (float): Support stiffness value.
    lambda_ (float): Load multiplier.
    p (numpy.ndarray): External force vector.

    Returns:
    numpy.ndarray: Residual vector at all degrees of freedom.
    """
    # Compute the internal forces
    qtot, _, _, _ = internal_forces(u, ndof, elemnoa, elemnob, EA, x0, l0, listfix, k_inf)

    # Compute the residual vector
    fhe = lambda_ * p - qtot

    return fhe


def internal_forces(u, ndof, elemnoa, elemnob, EA, x0, l0, listfix, k_inf):
    """
    Computes the internal forces for a structure.

    Parameters:
    u (numpy.ndarray): Displacement vector.
    ndof (int): Number of degrees of freedom.
    elemnoa (list): Starting node indices for each element.
    elemnob (list): Ending node indices for each element.
    EA (list): Axial stiffness of elements.
    x0 (list): Initial position vectors of elements.
    l0 (list): Lengths of elements.
    listfix (list): Indices of fixed degrees of freedom.
    k_inf (float): Support stiffness value.

    Returns:
    tuple: (qtot, qa, qb, u)
        qtot: Total internal forces vector.
        qa: Internal forces at start nodes of elements.
        qb: Internal forces at end nodes of elements.
        u: Displacement vector (unchanged from input).
    """
    nelem = len(elemnoa)
    qa = []
    qb = []
    du = []
    x = []
    epsG = []

    # Compute internal forces for each element
    for i in range(nelem):
        # Extract displacements for nodes A and B
        ua = u[[elemnoa[i] * 2 - 2, elemnoa[i] * 2 - 1]]  # Node A displacements
        ub = u[[elemnob[i] * 2 - 2, elemnob[i] * 2 - 1]]  # Node B displacements

        du_i = ub - ua
        du.append(du_i)

        # Compute the new position vector
        x_i = x0[i] + du_i
        x.append(x_i)

        # Compute Green's strain
        epsG_i = (np.dot(x0[i], du_i) + 0.5 * np.dot(du_i, du_i)) / l0[i] ** 2
        epsG.append(epsG_i)

        # Internal forces at nodes A and B
        qa_i = -EA[i] / l0[i] * epsG_i * x_i
        qb_i = EA[i] / l0[i] * epsG_i * x_i
        qa.append(qa_i)
        qb.append(qb_i)

    # Assemble total internal forces vector
    qtot = np.zeros(ndof)
    for i in range(nelem):
        idofa = [elemnoa[i] * 2 - 2, elemnoa[i] * 2 - 1]
        idofb = [elemnob[i] * 2 - 2, elemnob[i] * 2 - 1]
        qtot[idofa] += qa[i]
        qtot[idofb] += qb[i]

    # Add reaction forces for large supports
    reac = np.zeros(ndof)
    reac[listfix-1] = k_inf * u[listfix-1]

    qtot += reac

    return qtot, qa, qb, u


def stiffness_matrix(u, ndof, listfix, elemnoa, elemnob, EA, x0, l0):
    """
    Computes the stiffness matrix and its components for a structure.

    Parameters:
    u (numpy.ndarray): Displacement vector.
    ndof (int): Number of degrees of freedom.
    listfix (list): Indices of fixed degrees of freedom.
    elemnoa (list): Starting node indices for each element.
    elemnob (list): Ending node indices for each element.
    EA (list): Axial stiffness of elements.
    x0 (list): Initial position vectors of elements.
    l0 (list): Lengths of elements.

    Returns:
    tuple: (Ktot, k_inf, K0, Ku, Ks, K0tot, Kutot, Kstot)
        Ktot: Total stiffness matrix.
        k_inf: Support stiffness value.
        K0: Initial stiffness matrix contributions.
        Ku: Stiffness matrix contributions from geometric changes.
        Ks: Stiffness matrix contributions from internal forces.
        K0tot, Kutot, Kstot: Assembled stiffness matrices for K0, Ku, and Ks.
    """

    nelem = len(elemnoa)
    K = []
    K0 = []
    Ku = []
    Ks = []
    du = []
    x = []
    epsG = []
    N = []

    # Compute internal forces and stiffness matrices for each element
    for i in range(nelem):
        # Extract displacements for nodes A and B
        ua = u[[elemnoa[i] * 2 - 2, elemnoa[i] * 2 - 1]]  # Node A displacements
        ub = u[[elemnob[i] * 2 - 2, elemnob[i] * 2 - 1]]  # Node B displacements

        du_i = ub - ua
        du.append(du_i)

        # Compute the new position vector
        x_i = x0[i] + du_i
        x.append(x_i)

        # Compute Green's strain
        epsG_i = (np.dot(x0[i], du_i) + 0.5 * np.dot(du_i, du_i)) / l0[i] ** 2
        epsG.append(epsG_i)

        # Compute axial force
        N_i = EA[i] * epsG_i
        N.append(N_i)

        # Initial stiffness matrix
        tmp = np.outer(x0[i], x0[i])
        K0_i = EA[i] / l0[i] ** 3 * np.block([[tmp, -tmp], [-tmp, tmp]])
        K0.append(K0_i)

        # Stiffness matrix contribution from geometric changes
        tmp = np.outer(x0[i], du_i) + np.outer(du_i, x0[i]) + np.outer(du_i, du_i)
        Ku_i = EA[i] / l0[i] ** 3 * np.block([[tmp, -tmp], [-tmp, tmp]])
        Ku.append(Ku_i)

        # Stiffness matrix contribution from internal forces
        tmp = np.eye(2)
        Ks_i = N_i / l0[i] * np.block([[tmp, -tmp], [-tmp, tmp]])
        Ks.append(Ks_i)

        # Total stiffness matrix for the element
        K_i = K0_i + Ku_i + Ks_i
        K.append(K_i)

    # Assemble global stiffness matrix
    Ktot = assemble(K, ndof, elemnoa, elemnob)

    # Create a support condition matrix (penalization method)
    k_inf = 10**6 * np.max(np.diag(Ktot))
    K_support = np.zeros((ndof, ndof))
    np.fill_diagonal(K_support, [k_inf if i+1 in listfix else 0 for i in range(ndof)])

    # Add support conditions to the total stiffness matrix
    Ktot += K_support

    # Assemble additional stiffness matrices if needed
    K0tot = assemble(K0, ndof, elemnoa, elemnob) + K_support
    Kutot = assemble(Ku, ndof, elemnoa, elemnob)
    Kstot = assemble(Ks, ndof, elemnoa, elemnob)

    return Ktot, k_inf, K0, Ku, Ks, K0tot, Kutot, Kstot


def assemble(K, ndof, elemnoa, elemnob):
    """
    Assembles the global stiffness matrix from element stiffness matrices.

    Parameters:
    K (list): Element stiffness matrices.
    ndof (int): Number of degrees of freedom.
    elemnoa (list): Starting node indices for each element.
    elemnob (list): Ending node indices for each element.

    Returns:
    numpy.ndarray: Global stiffness matrix.
    """
    Ktot = np.zeros((ndof, ndof))

    for i in range(len(K)):
        idofa = [elemnoa[i] * 2 - 2, elemnoa[i] * 2 - 1]
        idofb = [elemnob[i] * 2 - 2, elemnob[i] * 2 - 1]
        Ke = K[i]
        Ktot[np.ix_(idofa, idofa)] += Ke[:2, :2]
        Ktot[np.ix_(idofa, idofb)] += Ke[:2, 2:]
        Ktot[np.ix_(idofb, idofa)] += Ke[2:, :2]
        Ktot[np.ix_(idofb, idofb)] += Ke[2:, 2:]

    return Ktot


def plot_structure(xnod, ynod, elemnoa, elemnob, u=None, drawing=None):
    """
    Plots the initial configuration (if `u` is not specified) or the deformed configuration.

    Parameters:
    xnod (list or numpy.ndarray): x-coordinates of nodes.
    ynod (list or numpy.ndarray): y-coordinates of nodes.
    elemnoa (list or numpy.ndarray): Starting node indices of elements.
    elemnob (list or numpy.ndarray): Ending node indices of elements.
    u (list or numpy.ndarray, optional): Displacement vector (alternating x, y displacements).
    """
    # Use the default Matplotlib color cycle
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Plot initial nodes
    plt.plot(xnod, ynod, 'o', label='Original Nodes', color=colors[0])

    # Draw initial elements
    for i in range(len(elemnoa)):
        plt.plot(
            [xnod[elemnoa[i]-1], xnod[elemnob[i]-1]],
            [ynod[elemnoa[i]-1], ynod[elemnob[i]-1]],
            color=colors[0],
            label='Initial Configuration' if i == 0 else ""
        )

    # If displacements `u` are provided, plot deformed configuration    
    if u is not None:
        xnod_deformed = np.array(xnod) + u[0::2]  # Add x displacements
        ynod_deformed = np.array(ynod) + u[1::2]  # Add y displacements        
        if drawing is None:
            drawing = {}
            for i in range(len(elemnoa)):
                drawing[i] = plt.plot(
                    [xnod_deformed[elemnoa[i]-1], xnod_deformed[elemnob[i]-1]],
                    [ynod_deformed[elemnoa[i]-1], ynod_deformed[elemnob[i]-1]],
                    color=colors[1],
                    linestyle='--',
                    label='Deformed Configuration' if i == 0 else ""
                )
        else:
            for i in range(len(elemnoa)):
                drawing[i][0].set_xdata([xnod_deformed[elemnoa[i]-1], xnod_deformed[elemnob[i]-1]])
                drawing[i][0].set_ydata([ynod_deformed[elemnoa[i]-1], ynod_deformed[elemnob[i]-1]])

    # Format the plot
    plt.axis('equal')
    plt.xlabel('X-coordinate')
    plt.ylabel('Y-coordinate')
    plt.title('Structure Configuration')
    plt.grid(True)    

    return drawing

def run_static_analysis(opts, ndof, xnod, ynod, elemnoa, elemnob, EA, x0, l0, 
                 listfix, k_inf, p, dofload, dofobs, drawing):
    """
    Runs structural analysis using specified solver options.

    Parameters:
    opts (dict): Options for the solver and analysis.
    ndof (int): Number of degrees of freedom.
    xnod, ynod (list): Node coordinates.
    elemnoa, elemnob (list): Element connectivity.
    EA (list): Axial stiffness of elements.
    x0, l0 (list): Initial positions and lengths of elements.
    listfix (list): Fixed degrees of freedom.
    k_inf (float): Support stiffness penalty factor.
    p (numpy.ndarray): External force vector.
    dofload (int): Degree of freedom where load is applied.
    dofobs (int): Degree of freedom to observe during analysis.
    drawing (dict): Drawing objects for the deformed configuration.
    """
    u_init = np.zeros(ndof)

    # Linear analysis before first step
    lambda_ = 1  # Load multiplier
    Ktot, k_inf, K0, Ku, Ks, K0tot, _, _ = stiffness_matrix(u_init, ndof, listfix, elemnoa, elemnob, EA, x0, l0)
    u_lin = np.linalg.solve(K0tot, lambda_ * p)

    if opts['solver'] == 'Imposed force - Python fsolve':
        u_init = np.zeros(ndof)

        for lambda_ in opts['lambda']:
            # Solve using fsolve
            u_new = fsolve(
                lambda u: fhe(u, ndof, elemnoa, elemnob, EA, x0, l0, listfix, k_inf, lambda_, p),
                u_init,
                xtol=1e-6
            )

            # Plot deformed configuration and force-displacement curve
            plt.subplot(1, 2, 1)
            plot_structure(xnod, ynod, elemnoa, elemnob, u_new, drawing=drawing)
            plt.title("Deformed Configuration")
            plt.subplot(1, 2, 2)
            plt.plot(u_new[dofload], lambda_, 'o', color='blue')
            plt.xlabel("Displacement")
            plt.ylabel("Force")
            plt.grid(True)

            # Update initial guess
            u_init = u_new
            plt.pause(0.01)
            plt.clf()  # Clear the figure for the next iteration

    elif opts['solver'] == 'Continuation Arc length':
        prev_du = u_lin
        R = opts['R']

        # Current position of the center of the hypersphere
        u0 = u_init
        lambda0 = 0

        solution = []
        solution.append({'u': u0, 'lambda': lambda0})

        # Iterate over steps
        u_i = u0
        lambda_i = lambda0

        for istep in range(1, opts['nstep'] + 1):
            print(f"Step = {istep}")

            converged = False
            iter_count = 0

            while not converged:
                # Compute internal forces and stiffness matrix
                f_int, _, _, _ = internal_forces(u_i, ndof, elemnoa, elemnob, EA, x0, l0, listfix, k_inf)
                Ktot, _, _, _, _, _, _, _ = stiffness_matrix(u_i, ndof, listfix, elemnoa, elemnob, EA, x0, l0)

                y = np.linalg.solve(Ktot, p)
                z = np.linalg.solve(Ktot, f_int)
                w = u_i - u0 - z

                # Solve quadratic equation for lambda
                coeffs = [1 + np.dot(y, y), 2 * (np.dot(w, y) - lambda0), np.dot(w, w) + lambda0**2 - R**2]
                lambda_np1 = np.roots(coeffs)

                # Choose the correct solution for lambda
                u_a = u_i + lambda_np1[0] * y - z
                u_b = u_i + lambda_np1[1] * y - z

                if np.dot((u_a - u0), prev_du) > np.dot((u_b - u0), prev_du):
                    new_un = u_a
                    new_lambdan = lambda_np1[0]
                else:
                    new_un = u_b
                    new_lambdan = lambda_np1[1]

                iter_count += 1
                if iter_count > 40:
                    converged = True
                if abs(new_lambdan - lambda_i) / abs(new_lambdan) < 1e-6 and np.linalg.norm(new_un - u_i) / np.linalg.norm(u_i) < 1e-6:
                    converged = True

                print(f"  Iter {iter_count}: rel norm(du)={np.linalg.norm(new_un - u_i) / np.linalg.norm(new_un):.3f}, dlambda={new_lambdan - lambda_i:.3f}")

                lambda_i = new_lambdan
                u_i = new_un

            # Plot results
            plt.subplot(1, 2, 1)
            plot_structure(xnod, ynod, elemnoa, elemnob, u_i, drawing=drawing)
            plt.title("Deformed Configuration (Arc Length)")
            plt.subplot(1, 2, 2)
            for i in range(0, len(dofobs)):
                clr = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
                plt.plot(u_i[dofobs[i]-1], lambda_i, 'o', markersize=4, color=clr[i])
            plt.xlabel("Displacement")
            plt.ylabel("Force")
            plt.grid(True)

            plt.pause(0.001)
            
            # Store solution
            solution.append({'u': u_i, 'lambda': lambda_i})
            prev_du = u_i - u0

            # Update the center of the hypersphere
            u0 = u_i
            lambda0 = lambda_i


if __name__ == "__main__":
    
    plt.ion()

    # Load structure
    xnod, ynod, elemnoa, elemnob, ndof, EA, dofload, p, listfix, dofobs = ssd.example_structures.load_structure(4)

    # Plot initial configuration
    plt.figure()
    drawing = plot_structure(xnod, ynod, elemnoa, elemnob)

    # Initialize x0 and l0 for each element
    x0, l0 = initialize(xnod, ynod, elemnoa, elemnob)

    # Elastic step: run a 1st order linear elastic step
    plt.figure()
    u_init = np.zeros(ndof)
    lambda_val = 0.01  # Small load multiplier

    # Calculate the stiffness matrix and related parameters
    Ktot, k_inf, K0, Ku, Ks, K0tot, _, _ = stiffness_matrix(u_init, ndof, listfix, elemnoa, elemnob, EA, x0, l0)
    
    # Solve and compute nodal displacements
    u_lin = np.linalg.solve(K0tot, lambda_val * p)

    # Display deformed configuration
    plt.subplot(1, 2, 1)
    drawing = plot_structure(xnod, ynod, elemnoa, elemnob, u_lin)
    
    # Compute internal forces
    f_int, f_a, f_b, u = internal_forces(u_lin, ndof, elemnoa, elemnob, EA, x0, l0, listfix, k_inf)

    # Nonlinear static analysis
    opts = {
        'solver': 'Continuation Arc length',
        'R': 0.1,
        'nstep': 100
    }

    run_static_analysis(opts, ndof, xnod, ynod, elemnoa, elemnob, EA, x0, l0,
                        listfix, k_inf, p, dofload, dofobs, drawing)

    plt.show(block=True)

