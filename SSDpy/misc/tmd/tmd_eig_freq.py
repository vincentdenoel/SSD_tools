import numpy as np
# Compute the ev and V of a plate with four springs
# Define the problem parameters:
# k = spring stiffnesses [N/m]
# m = mass surface density [kg/m2]
# b,L = plate width and length [m]

k = [1.05, 1., 1.05, 1.1]
m = 1.
b = 1.
L = 1;

# Define the equations
k11 = k[0] + k[1] + k[2] + k[3]
k22 = b**2 / 4 * k11
k33 = L**2 / 4 * k11
k12 = b / 2   * (k[0] - k[1] - k[2] + k[3])
k13 = L / 2   * (k[0] + k[1] - k[2] - k[3])
k23 = b*L / 2 * (k[0] - k[1] + k[2] - k[3])

# Define matrices K and M
K = np.array([
    [k11, k12, k13],
    [k12, k22, k23],
    [k13, k23, k33]
])

M = m * np.array([
    [1, 0, 0],
    [0, b**2 / 12, 0],
    [0, 0, L**2 / 12]
])

# Calculate eigenvalues and eigenvectors
ev, V = np.linalg.eig(np.linalg.solve(M, K))

# Sort the eigenvalues and eigenvectors, from lowest to highest
idx = np.argsort(ev)
ev = ev[idx]
V = V[:,idx]

# Diagonalize K and M
K_diag = V.T @ K @ V
M_diag = V.T @ M @ V

# Print results with 3 significant digits
print("Total mass [kg]: ", m * b * L)
for imode in range(3):
    print("Mode ", imode)
    print("  Frequency: ", np.round(np.sqrt(ev[imode]),3), "rad/s, ", np.round(np.sqrt(ev[imode]) / (2 * np.pi),3), "Hz")
    print("  Mass [kg]: ", np.round(M_diag[imode, imode],3))
    print("  V: ", V[:, imode])


