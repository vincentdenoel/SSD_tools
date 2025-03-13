import numpy as np
import matplotlib.pyplot as plt

def histv(X, N=None):
    """
    Computes the probability density function as a stepped histogram.
    
    Parameters:
        X (array-like): Signal data.
        N (int, optional): Number of bins. If not provided, it is calculated using `NbBin`.
        
    Returns:
        tuple: (nn, xx)
            - nn (array): Values of the probability density function.
            - xx (array): Bin edges for the stepped histogram.
    """
    if N is None:
        N = NbBin(X)

    # Handle the case where the standard deviation is zero
    if np.std(X) == 0:
        X = X + (np.finfo(float).eps + np.mean(X)) * 1e-8 * np.random.randn(*X.shape)
        N = NbBin(X)

    # Compute histogram
    n, x = np.histogram(X, bins=N)
    dx = x[1] - x[0]

    # Normalize histogram to compute probability density
    n = n / (np.sum(n) * dx)

    # Create stepped histogram
    xx = [x[0] - dx / 2]
    nn = [0]

    for i in range(N):
        xx.append(x[i] - dx / 2)
        xx.append(x[i] + dx / 2)
        nn.append(n[i])
        nn.append(n[i])

    xx.append(x[-1] + dx / 2)
    nn.append(0)

    xx = np.array(xx)
    nn = np.array(nn)

    # If called without output, plot the histogram
    if not plt.isinteractive():
        plt.plot(xx, nn)
        plt.xlabel("Value")
        plt.ylabel("Probability Density")
        plt.show()

    return nn, xx


def NbBin(x):
    """
    Computes the optimal number of bins based on Scott's rule.
    
    Parameters:
        x (array-like): Input data.
        
    Returns:
        int: Optimal number of bins.
    """
    N = len(x)
    W = 3.5 * np.nanstd(x) / N**(1/3)
    Nbins = int(np.floor((np.max(x) - np.min(x)) / W))
    return max(Nbins, 2)
