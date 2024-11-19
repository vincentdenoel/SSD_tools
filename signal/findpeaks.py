import numpy as np

def findpeaks(x):
    """Find the local maxima and minima of a signal x.
    This is a function of the SSDpy package.
    
    Returns: the maxima, their indices, the minima, and their indices."""

    x = np.array(x)
    dx = np.diff(x)
    iMax = np.where((dx[:-1] > 0) & (dx[1:] < 0))[0] + 1
    iMin = np.where((dx[:-1] < 0) & (dx[1:] > 0))[0] + 1
    M = x[iMax]
    m = x[iMin]
    return M, iMax, m, iMin
