import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from . import histv

def histv_plots(x, ax=None, bins='auto', density=True, alpha=0.5, color='C0',edgecolor=None):
    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 3))
    if ax is not None and len(ax) != 2:
        plot_both = False
    else:
        plot_both = True
    if isinstance(ax, plt.Axes):
        ax = [ax]

    # we're doing the plot manually to avoid patches outside plotting area
    # get histogram:
    nn, xx = histv(x)
    maxh = np.max(nn)
    minh = np.min(nn[nn > 0])
    nn[np.where(nn==0)] = 0.999 *minh

    points = np.column_stack((xx, nn))

    if plot_both:
        # plot a patch
        poly = Polygon(points, closed=True, facecolor=color, edgecolor=edgecolor, alpha=alpha)
        ax[0].add_patch(poly)        
        ax[0].set_title('Histogram')
        ax[0].autoscale_view()
        yl = ax[0].get_ylim()
        ax[0].set_ylim([0, yl[1]])
        

        poly2 = Polygon(points, closed=True, facecolor=color, edgecolor=edgecolor, alpha=alpha)
        ax[1].add_patch(poly2)
        ax[1].set_title('Histogram')
        ax[1].set_yscale('log')
        ax[1].autoscale_view()
        yl = ax[1].get_ylim()
        ax[1].set_ylim([minh, yl[1]])
        
        

    else:
        poly = Polygon(points, closed=True, facecolor=color, edgecolor=edgecolor, alpha=alpha)
        ax[0].add_patch(poly)        
        ax[0].set_title('Histogram')
        yl = ax[0].get_ylim()
        ax[0].set_ylim([0, yl[0]])
        ax[0].autoscale_view()


