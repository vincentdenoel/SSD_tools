import os
import numpy as np
import matplotlib.pyplot as plt
# from scipy.fftpack import fft

# imports from SSDpy
from SSDpy.signal import signal

def plot_tf(t, x, NFFT=None, sp=None, scaling=1, *args, **kwargs):
    # set default values for optional arguments
    if sp is None:
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
        sp = (ax1, ax2)
    if NFFT is None:
        NFFT = 512

    # make sure time runs along columns
    if x.shape[1] > x.shape[0]:
        x = x.T
    if t.shape[1] > t.shape[0]:
        t = t.T


    dt = t[1] - t[0]

    XX, f, *useless = signal.vpsd(x, NFFT, dt, scaling)

    psdX = np.zeros((XX.shape[0], XX.shape[1]))
    for i in range(x.shape[1]):
        psdX[:, i] = XX[:, i, i]
        psdX[0, i] = np.nan

    sp[0].plot(t, x, *args, **kwargs)
    sp[0].set_xlabel('Time [s]')
    sp[0].grid()

    sp[1].semilogy(f, psdX)
    sp[1].set_xlabel('Frequency [Hz]')
    sp[1].grid()

    plt.show(block=False)

    if 'return' in kwargs:
        if kwargs['return'] == 'psd':
            return f, psdX
        elif kwargs['return'] == 'all':
            return f, psdX, t, x
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
