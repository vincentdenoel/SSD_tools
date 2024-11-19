import numpy as np
import matplotlib.pyplot as plt

# imports from SSDpy
from SSDpy.signal import vpsd


def plot_tf(t, x, NFFT=None, sp=None, scaling=1, output=None, freq_xlim=None, *args, **kwargs):
    """ 
    This function plots the time series and the power spectral density of a signal.
    This function is part of the SSDpy package
    """

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

    XX, f, *useless = vpsd(x, NFFT, dt, scaling)

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
