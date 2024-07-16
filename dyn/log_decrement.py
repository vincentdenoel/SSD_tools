import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

from SSDpy.tools import tools as ssdtools


def log_decrement(t, y, t0=None, t1=None, DoPlot=True, F=None):
    """Perform log decrement analysis on a signal.

    Inputs: t, y, t0, t1, DoPlot, F
    t : time (numpy.ndarray)
    y : signal (numpy.ndarray)
    t0 : start time (float) - optional, if not provided, the start time of the signal is used
    t1 : end time (float)   - optional, if not provided, the end time of the signal is used
    DoPlot : plot results (bool)
    F : filter frequency (float)

    Outputs: T, xi
    """
    dt = t[1] - t[0]
    fs = 1 / dt
    
    if t0 is None:
        t0 = t[0]  # start time
    if t1 is None:
        t1 = t[-1]  # end time
    if F is not None:
        B, A = butter(4, F / (fs / 2))
        y_ = filtfilt(B, A, y[::-1])
        y = y_[::-1]
    
    if t0<t[0]:
        t0 = t[0]
    if t1>t[-1]:
        t1 = t[-1] 
        
    i0 = int(np.round((t0 - t[0]) / dt))
    i1 = int(np.round((t1 - t[0]) / dt))
       
    v = y[i0:i1]
    v = v - np.mean(v)
    tw = t[i0:i1]
    
    # find peaks in v
    M, iM, *useless = ssdtools.find_peaks(v)
    m, im, *useless = ssdtools.find_peaks(-v)
    
    MM = M.copy()
    mm = m.copy()
    
    # Remove data before the first peak
    pM = np.argmax(M)
    M = M[pM:]
    MM = MM[pM:]
    iM = iM[pM:]
    
    pm = np.argmax(m)
    m = m[pm:]
    mm = mm[pm:]
    im = im[pm:]
    
    # Find first negative value in M and discard remaining data
    i = np.argmax(M < 0)
    if i > 0:
        M = M[:i]
        MM = MM[:i]
        iM = iM[:i]

    # Find first negative value in m and discard remaining data
    i = np.argmax(m < 0)
    if i > 0:
        m = m[:i]
        mm = mm[:i]
        im = im[:i]

    M = M / M[0]
    m = m / m[0]
    iiM = np.arange(len(M))
    iim = np.arange(len(m))
    
    xiM = np.zeros(len(iiM))
    for i in range(1, len(iiM)):
        P = np.polyfit(iiM[:i+1], np.log(M[:i+1]), 1)
        xiM[i] = -P[0] / (2 * np.pi)
    
    xim = np.zeros(len(iim))
    for i in range(1, len(iim)):
        P = np.polyfit(iim[:i+1], np.log(m[:i+1]), 1)
        xim[i] = -P[0] / (2 * np.pi)
    
    Method = 3
    if Method == 1:
        indx = max([1, len(xiM) - 6])
        xi1 = np.median(xiM[indx:])
        xi2 = np.median(xim[indx:])
        xi = (xi1 + xi2) / 2
    elif Method == 2:
        skip1 = min([15, len(M)])
        a = np.zeros(len(M))
        r = np.zeros(len(M))
        for j in range(skip1, len(M)):
            P = np.polyfit(iiM[:j+1], np.log(M[:j+1]), 1)
            a[j] = -P[0]
            r[j] = np.linalg.norm(np.log(np.abs(M[:j+1])) - P[1] - P[0] * iiM[:j+1]) / (j+1)
        imin1 = np.argmin(r[skip1:]) + skip1
        xi1 = a[imin1] / (2 * np.pi)
        
        skip2 = min([15, len(m)])
        a = np.zeros(len(m))
        r = np.zeros(len(m))
        for j in range(skip2, len(m)):
            P = np.polyfit(iim[:j+1], np.log(m[:j+1]), 1)
            a[j] = -P[0]
            r[j] = np.linalg.norm(np.log(np.abs(m[:j+1])) - P[1] - P[0] * iim[:j+1]) / (j+1)
        imin2 = np.argmin(r[skip2:]) + skip2
        xi2 = a[imin2] / (2 * np.pi)
        xi = (xi1 + xi2) / 2
    elif Method == 3:
        npts = 3
        sz = min([npts, len(M) // 4])
        xiM = np.full(len(M), np.nan)
        for j in range(sz+1, len(M)-sz-1):
            ii = np.arange(j-sz, j+sz+1)
            P = np.polyfit(iiM[ii], np.log(M[ii]), 1)
            xiM[j] = -P[0] / (2 * np.pi)
        xi1 = np.nanmedian(xiM)
        
        sz = min([npts, len(m) // 4])
        xim = np.full(len(m), np.nan)
        for j in range(sz+1, len(m)-sz-1):
            ii = np.arange(j-sz, j+sz+1)
            P = np.polyfit(iim[ii], np.log(m[ii]), 1)
            xim[j] = -P[0] / (2 * np.pi)
        xi2 = np.nanmedian(xim)
        xi = (xi1 + xi2) / 2
    
    n = int(np.ceil(0.05 * len(M)))
    i_ = iM[n:len(M)-n]
    T = np.mean(np.diff(i_) * dt)
    
    if DoPlot:
        logdec_plot(t, y, T, xi, M,iM,m,im,iiM,iim,i0,i1,xiM,xim)
    
    return T, xi


def logdec_plot(t, y, T, xi, M,iM,m,im,iiM,iim,i0,i1,xiM,xim):
    
    cmp = plt.rcParams['axes.prop_cycle'].by_key()['color']
    dt = t[1] - t[0]
    
    iM = iM + i0
    im = im + i0  

    # plot detailed results
    plt.figure(figsize=(10, 8))
    
    plt.subplot(2, 2, 1)
    plt.plot(t, y)
    plt.plot(t[iM], y[iM], '.-', color=cmp[1])
    plt.plot(t[im], y[im], '.-', color=cmp[2])
    plt.xlabel('time [s]')
    plt.ylabel('Max & Min')
    
    plt.subplot(2, 2, 2)
    plt.plot(iiM, np.log(M), '.', color=cmp[1])
    plt.plot(iim, np.log(m), '.', color=cmp[2])
    plt.xlabel('Index of Max/Min')
    plt.ylabel('Log of normalized Max & Min')
    
    plt.subplot(2, 2, 3)
    plt.plot(iM[1:]*dt, xiM[1:], 'x-', color=cmp[1])
    plt.plot(im[1:]*dt, xim[1:], 'x-', color=cmp[2])
    plt.xlabel('time[s]')
    plt.ylabel('Damping ratio xi []')
    plt.title(f'Natural frequency {1/T:.3f} Hz')
    plt.plot(plt.xlim(), [xi, xi], 'k:', linewidth=2)
    
    plt.subplot(2, 2, 4)
    plt.plot( y[iM[1:]], xiM[1:], color=cmp[1])
    plt.plot(-y[im[1:]], xim[1:], color=cmp[2])
    plt.xlabel('Amplitude')
    plt.ylabel('Damping ratio xi []')
    plt.plot(plt.xlim(), [xi, xi], 'k:', linewidth=2)
        
    plt.figure()
    ii = np.arange(i0, i1)
    y0 = np.mean(y)
    y = y - y0
    jj = np.arange(max(1, i0 - 200), min(len(y), i1 + 200))
    plt.plot(t[jj], y[jj], color=[0.7, 0.7, 0.7])
    plt.plot(t[ii], y[ii], color=cmp[0])
    
    pred = y[iM[0]] * np.exp(-xi*2*np.pi/T*(t[ii]-t[iM[0]]))

    plt.plot(t[ii], pred, 'k', linewidth=2)
    plt.xlabel('time [s]')
    plt.title(f'Frequ.: {1/T:.2f} Hz,   Damping: {np.mean(xi)*100:.2f}%')
    plt.show()
