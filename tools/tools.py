import os
import numpy as np
import matplotlib.pyplot as plt

# imports from SSDpy
from SSDpy.signal import signal

""" This is a module of the SSDpy package. It contains several tools to do intersting things for Structural and Stochastic Dynamics
Contents of this package:
- plot_tf      Plot time and frequency domain of a signal
- saveas       Save a Matplotlib figure as a PDF file with vector graphics
- find_peaks   Find peaks in a signal
"""

def plot_tf(t, x, NFFT=None, sp=None, scaling=1, output=None, freq_xlim=None, *args, **kwargs):
    # set default values for optional arguments
    if sp is None:
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
        sp = (ax1, ax2)
    if NFFT is None:
        NFFT = 512

    # Make sure x is organized as: nb of signals, nb of time steps
    if x.ndim == 1:
        x = x.reshape(1, -1) # reshape to 2D array
    else:
        if x.shape[1] < x.shape[0]:
            x = x.T

    dt = t[1] - t[0]

    XX, f, *useless = signal.vpsd(x, NFFT, dt, scaling)

    psdX = np.zeros((XX.shape[0], XX.shape[1]))
    for i in range(x.shape[0]):
        psdX[:, i] = XX[:, i, i]
        psdX[0, i] = np.nan

    sp[0].plot(t, x.T, *args, **kwargs)
    sp[0].set_xlabel('Time [s]')
    sp[0].grid()

    sp[1].semilogy(f, psdX)
    sp[1].set_xlabel('Frequency [Hz]')
    sp[1].grid()

    if freq_xlim is not None:
        sp[1].set_xlim(freq_xlim)

    plt.show(block=False)

    if output == 'psd':
        return f, psdX
    elif output == 'all':
        return f, psdX, t, x
    elif output == 'plots':
        return sp
    elif output == 'psd-plots':
        return f, psdX, sp
    else:
        return



def saveas(fig, folder_path, filename):
    """
    Saves a Matplotlib figure as a PDF file with vector graphics.

    Parameters:
        fig (matplotlib.figure.Figure): The Matplotlib figure to save.
        folder_path (str): The path to the directory where the PDF file should be saved.
        filename (str): The name of the PDF file.
    """
    # Ensure that the directory exists
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    
    # Set pdf.fonttype to 'embedded'
    plt.rcParams['pdf.fonttype'] = 42

    # Save the figure as a PDF file
    file_path = os.path.join(folder_path, filename)
    fig.savefig(file_path, format='pdf')



def find_peaks(x):
    """Find the local maxima and minima of a signal x.
    
    Returns: the maxima, their indices, the minima, and their indices."""

    x = np.array(x)
    dx = np.diff(x)
    iMax = np.where((dx[:-1] > 0) & (dx[1:] < 0))[0] + 1
    iMin = np.where((dx[:-1] < 0) & (dx[1:] > 0))[0] + 1
    M = x[iMax]
    m = x[iMin]
    return M, iMax, m, iMin