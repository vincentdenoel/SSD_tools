import numpy as np

def load_structure(cas):
    # Define structure (geometry, material, supports)
    
    if cas == 1:
        # Node coordinates - there are 3 nodes (0,1), (1,1) and (1,0)
        # Nodes are used to number dofs. At node 1: 1&2, at node 2: 3&4, etc.
        xnod = np.array([0, 1, 2])
        ynod = np.array([0, 1, 0])

        # Elements - there are two elements (1>2) and (2>3). They link nodes 1&2 and nodes 2&3.
        elemnoa = np.array([1, 2])
        elemnob = np.array([2, 3])

        ndof = 2 * len(xnod)   # Total number of dofs (fixed + free)
        nelem = len(elemnoa)   # Total number of elements

        # Define axial stiffness for each element
        EA = np.array([1, 1])

        # There is one load, imposed at dof 3
        dofload = 4

        # Create vector of loading
        p = np.zeros(ndof)
        p[dofload - 1] = 1  # Adjust for 0-based indexing

        listfix = np.array([1, 2, 5, 6])  # List of fixed dofs
        listfree = np.setdiff1d(np.arange(1, ndof + 1), listfix)  # List of free dofs

        # Which degrees-of-freedom to plot
        dofobs = 4
        
    elif cas == 2:
        # Node coordinates - there are 3 nodes (0,1), (1,1) and (1,0)
        a = 0.4
        b = 1
        xnod = np.array([0, 1, 2, 1]) * b
        ynod = np.array([0, 1, 0, 2]) * a

        # Elements - there are two elements (1>2) and (2>3). They link nodes 1&2 and nodes 2&3.
        elemnoa = np.array([1, 2, 1, 4, 2])
        elemnob = np.array([2, 3, 4, 3, 4])

        ndof = 2 * len(xnod)   # Total number of dofs (fixed + free)
        nelem = len(elemnoa)   # Total number of elements

        # Define axial stiffness for each element
        EA = np.array([1, 1, 1, 1, 1])

        # There is one load, imposed at dof 3
        dofload = 5

        # Create vector of loading
        p = np.zeros(ndof)
        p[dofload - 1] = 1

        listfix = np.array([1, 2, 6])  # List of fixed dofs
        listfree = np.setdiff1d(np.arange(1, ndof + 1), listfix)  # List of free dofs

        # Which degrees-of-freedom to plot
        dofobs = 5

    elif cas == 3:
        # Node coordinates
        xnod = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2])
        ynod = np.array([0, 1, 2, 3, 0, 1, 2, 3, 2, 3])

        # Elements - there are 16 elements defined by node pairs
        elemnoa = np.array([1, 2, 3, 5, 6, 7, 1, 6, 2, 7, 3, 8, 7, 9, 10, 7])
        elemnob = np.array([2, 3, 4, 6, 7, 8, 6, 2, 7, 3, 8, 4, 9, 10, 8, 10])

        ndof = 2 * len(xnod)   # Total number of dofs (fixed + free)
        nelem = len(elemnoa)   # Total number of elements

        # Define axial stiffness for each element
        EA = np.ones(16)

        # There is one load, imposed at dof 3
        dofload = 18

        # Create vector of loading
        p = np.zeros(ndof)
        p[dofload - 1] = -1

        listfix = np.array([1, 2, 9, 10])  # List of fixed dofs
        listfree = np.setdiff1d(np.arange(1, ndof + 1), listfix)  # List of free dofs

        # Which degrees-of-freedom to plot
        dofobs = 18

    elif cas == 4:
        # Node coordinates
        xnod = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0.55])
        ynod = np.array([0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 4])

        # Elements
        elemnoa = np.array([1, 2, 3, 4, 6, 7, 8, 9, 6, 2, 7, 3, 8, 4, 9, 5, 10, 11])
        elemnob = np.array([2, 3, 4, 5, 7, 8, 9, 10, 2, 7, 3, 8, 4, 9, 5, 10, 11, 9])

        ndof = 2 * len(xnod)   # Total number of dofs (fixed + free)
        nelem = len(elemnoa)   # Total number of elements

        # Define axial stiffness for each element
        EA = np.ones(18)

        # There is one load, imposed at dof 3
        dofload = 22

        # Create vector of loading
        p = np.zeros(ndof)
        p[dofload - 1] = -1

        listfix = np.array([1, 2, 11, 12])  # List of fixed dofs
        
        # Which degrees-of-freedom to plot
        dofobs = [21, 22]
    
    return xnod, ynod, elemnoa, elemnob, ndof, EA, dofload, p, listfix, dofobs
