import matplotlib.pyplot as plt
import matplotlib as mpl


#plt.style.use('ggplot')
def plot_settings():
    # Set style options directly in the module
    mpl.rcParams.update({
        'font.size': 9,
        'figure.figsize': (10, 6),
        'lines.linewidth': 1,
        'grid.linestyle': '-'
    })
