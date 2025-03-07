import numpy as np
from scipy.signal import hilbert, butter, filtfilt, find_peaks
from scipy.fft import fft

def phase(x):
    """
    This function is part of the SSDpy toolbox

    phase computes phase of signal(s), or relative phase, with Hilbert
    transform of mirrored signals.

    Parameters:
    x : numpy array
        Input signal. Can be a 1D array (vector) or a 2D array (n-by-2 matrix).

    Returns:
    If x is a vector: returns the unwrapped phase.
    If x is an n-by-2 matrix: returns the phase shift between the two signals.
    """
    s = np.shape(x)

    if len(s) == 1:
        return getPhase(x)
    
    elif len(s) == 2:
        if s[1] > s[0]:
            x = x.T
            s = np.shape(x)
        if s[1] > 2:
            raise Warning('Phase: input argument should be a 1-by-n vector or a 2-by-n matrix')
            return None
        sigphase1 = getPhase(x[:, 0])
        sigphase2 = getPhase(x[:, 1])
        shift = sigphase2 - sigphase1
        return shift
    
    else:
        raise Warning('myPhase: input argument should be a 1-by-n vector or a 2-by-n matrix')
        return None
    

def getPhase(x):
    """
    getPhase computes the unwrapped phase of a signal using the Hilbert transform
    of mirrored signals.

    Parameters:
    x : numpy array
        Input signal.

    Returns:
    Unwrapped phase of the signal.
    """
    n1 = len(x)

    # FFT to get the dominant frequency
    A = np.abs(fft(x))
    iA = np.argmax(A)

    # Butterworth filter
    B, A = butter(4, 2 * iA / (n1 / 2))
    x_filtered = filtfilt(B, A, x)

    # Find peaks and mirroring signal
    ip, *_ = find_peaks(x_filtered)
    n2 = ip[-1]

    tmp = np.concatenate((np.flipud(x_filtered), x_filtered[:n2], np.flipud(x_filtered[:n2-1])))
    hilbTMP = hilbert(tmp)
    sigphase = np.unwrap(np.angle(hilbTMP))
    return sigphase[n1:2*n1]

# Example usage:
# x = np.random.randn(100)  # Replace with your signal
# phase = myPhase(x)
# print(phase)
