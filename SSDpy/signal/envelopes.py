import numpy as np
from scipy.interpolate import interp1d

# imports from SSDpy
from SSDpy.signal.findpeaks import findpeaks

def envelopes(s):
    """ This is a function of the SSDpy package.
    It computes the lower and upper envelopes of a signal."""

    M, imax, m, imin = findpeaks(s)

    # interpolate through the locals min and max so that output has same size as input
    f_min = interp1d(imin, s[imin], kind='quadratic', fill_value="extrapolate")
    f_max = interp1d(imax, s[imax], kind='quadratic', fill_value="extrapolate")
    
    return f_min(np.arange(len(s))), f_max(np.arange(len(s)))