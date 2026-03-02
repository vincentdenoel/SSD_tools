import pysindy as ps
import numpy as np

# Simple linear system
t = np.linspace(0, 10, 100)
x = np.sin(t).reshape(-1, 1)  # Ensure data is 2D
sindy = ps.SINDy()
sindy.fit(x, t=np.diff(t).mean())
sindy.print()
