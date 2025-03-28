import numpy as np
import matplotlib.pyplot as plt

import SSDpy as ssd

import numpy as np
import matplotlib.pyplot as plt

ssd.plot_settings()

# Exmple d'utilisation

plt.figure()

L = 120 # beam length [m]
n = 500 # number of elements
k = 100 # foundation stiffness [N/m3] = gamma_water * largeur de poutre
EI = 2000 # bending stiffness [N.m2]

p0 = 1000 # charge répartie constante [N/m)]
P1 = 1e4 # charge concentrée à mi-longeur [N]
x1 = L / 2 # position
P2 = 5e3 # charge concentrée sur le bord [N]
x2 = 0.95 * L # position

p = p0 * np.ones(n + 1) # load per unit length [N/m]

p[int(x1 / L * n)] += P1 / ( L / n)
p[int(x2 / L * n)] += P2 / ( L / n)


# use this code to solve the beam on elastic foundation problem, under load p
# return the displacement v and bending moment M
x, v, M = ssd.winkler.BeamElasticFoundation(EI, k, L, p)

plt.subplot(121)
plt.plot(x, v)
plt.title('Displacement [m]')

plt.subplot(122)
plt.plot(x, M)
plt.title('Bending moment [N.m]')

plt.show()