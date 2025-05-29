from matplotlib.colors import LinearSegmentedColormap

# This colormap is based on the "Strawberry Fields Forever" color palette

# continuous colorbar
strbry_cont = [(0.4,0.2,0.5),(0.1,0.475,0.4),(0.3,0.65,0.1),(0.95,0.6,0.7),(1,0.85,1)]
strbry_cont_map = LinearSegmentedColormap.from_list('smap', strbry_cont, N = 256)

# diverging colorbar (for differences)
strbry_div = [(0, 0.125, 0.3),(0.1, 0.35, 0.4),(0.35,0.6,0.25),(0.7,0.8,0.2),
 (1,1,1),(1,0.6,0.7),(0.85,0.3,0.6),(0.5, 0.1, 0.6),(0.2, 0, 0.3)]
strbry_div_map = LinearSegmentedColormap.from_list('smap', strbry_div, N = 256)

colormaps = [strbry_cont_map, strbry_div_map]



