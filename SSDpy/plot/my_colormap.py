import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap


# This colormap is based on the "Strawberry Fields Forever" color palette
# continuous colorbar
strbry_cont = [(0.4,0.2,0.5),(0.1,0.475,0.4),(0.3,0.65,0.1),(0.95,0.6,0.7),(1,0.85,1)]
strbry_cont_map = LinearSegmentedColormap.from_list('smap', strbry_cont, N = 256)

# diverging colorbar (for differences)
strbry_div = [(0, 0.125, 0.3),(0.1, 0.35, 0.4),(0.35,0.6,0.25),(0.7,0.8,0.2),
 (1,1,1),(1,0.6,0.7),(0.85,0.3,0.6),(0.5, 0.1, 0.6),(0.2, 0, 0.3)]
strbry_div_map = LinearSegmentedColormap.from_list('smap', strbry_div, N = 256)


strbry_seq = [ (1.00, 0.97, 1.00),  # very light pink-white
    (1.00, 0.85, 1.00),  # pale pink
    (0.95, 0.60, 0.70),  # strawberry pink
    (0.72, 0.30, 0.58),  # raspberry
    (0.40, 0.20, 0.50),  # purple
    (0.20, 0.00, 0.30),  # dark purple
]

strbry_seq_map = LinearSegmentedColormap.from_list("strawberry_sequential", strbry_seq, N=256,)

def get_cmap(map=1):
    """
    Set the colormap for the plot.
    """
    if map == 1:
        top = cm.get_cmap('Oranges_r', 128)
        bottom = cm.get_cmap('Blues', 128)
        newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                            bottom(np.linspace(0, 1, 128))))
        new_cmap = ListedColormap(newcolors, name='OrangeBlue')
    elif map == 2:
        new_cmap = strbry_cont_map
    elif map == 3:
        new_cmap = strbry_div_map
    elif map == 4:
        new_cmap = strbry_seq_map

    return new_cmap