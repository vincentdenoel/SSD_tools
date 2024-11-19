import numpy as np

def vpsd(x, NFFT, DT=None, scaling=1):
    """
    vpsd - Computes the psd matrix of signals x
    (Historically VincePSD.m)
    This is part of the SSDpy toolbox.

    Parameters:
        x : numpy.ndarray
            signals in x are divided in smaller segments with length NFFT. 
            These segments can't overlap.
            NOTE : x must be arranged as x(signal number, sampling) in such a way
            that the time is running along the columns.
        NFFT : int
            Length of the FFT window
        DT : float, optional
            time step (in seconds) of the input signal. If provided, it returns the frequency vector.
        scaling : int, optional
            Scaling type. Default is 1.

    Returns:
        PSD : numpy.ndarray
            psd matrix of signals x. An estimation of the corresponding frequencies can be obtained if the time step DT is given.
        FREQ : numpy.ndarray
            frequency vector, if the time step DT is given.
        OMEGA : numpy.ndarray
            angular frequency vector, if the time step DT is given.
    """
    
    # Make sure x is organized as: nb of signals, nb of time steps
    if x.ndim == 1:
        NS, N = 1, len(x)
        x = x.reshape(1, -1) # reshape to 2D array
    else:
        NS, N = x.shape


    if NS > N:
        x = x.T
        NS, N = x.shape
    
    # trim nan's at beginning of signal and nan's at end of signal
    for i in range(NS):
        if not np.isnan(x[i, :]).all():
            x[i, np.isnan(x[i, :])] = np.nanmean(x[i, :])

    if DT is None:
        scaling = 0
    
    Nblocs = N // NFFT

    PSD = np.zeros((NFFT // 2, NS, NS))

    # Build Hanning window
    t = np.linspace(0, 1, NFFT)
    hann = np.sin(np.pi * t) ** 2
    hann = np.hanning(NFFT)
    hann = np.ones(NFFT)
    W = np.sum(hann ** 2)

    # Loop along the blocs
    for i in range(Nblocs):
        for s1 in range(NS):
            # Extract segment for signal 1 and window it
            xx1 = x[s1, i * NFFT:(i + 1) * NFFT]
            xx1 = xx1 * hann
            XX1 = np.fft.fft(xx1 - np.mean(xx1))
            for s2 in range(NS):
                # Extract segment for signal 2 and window it
                xx2 = x[s2, i * NFFT:(i + 1) * NFFT]
                xx2 = xx2 * hann
                XX2 = np.fft.fft(xx2 - np.mean(xx2))
                periodogram = XX1 * np.conj(XX2)

                PSD[:, s1, s2] = PSD[:, s1, s2] + (periodogram[0:NFFT // 2]).real

    PSD = PSD / Nblocs / W

    if DT is not None:
        DF = 1 / (NFFT * DT)
        FREQ = np.arange(0, NFFT // 2) * DF
        OMEGA = 2*np.pi*FREQ
        if scaling == 1:
            PSD = PSD / NFFT / DF * 2  # scale PSD such that intergral of PSD over FREQ returns variance
        if scaling == 2:
            PSD = PSD / NFFT / DF * 2 / 4 / np.pi  # scale PSD such that
        return PSD, FREQ, OMEGA
    else:
        return PSD
