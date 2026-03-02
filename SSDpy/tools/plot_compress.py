import numpy as np
import matplotlib.pyplot as plt

def plot_compress(t, x, block_size=200, ax=None):
    N = len(x)
    block_size = 200
    n_blocks = N // block_size  # ensure it's divisible

    # Take min and max for each block
    min_vals = np.min(x[:n_blocks*block_size].reshape(n_blocks, block_size), axis=1)
    max_vals = np.max(x[:n_blocks*block_size].reshape(n_blocks, block_size), axis=1)

    # Interleave min and max so the result is 100x smaller than original
    compressed = np.empty(2*n_blocks)
    compressed[0::2] = min_vals
    compressed[1::2] = max_vals

    t_compressed = np.linspace(t[0], t[n_blocks*block_size-1], 2*n_blocks)

    if ax is None:
        plt.figure(figsize=(10, 4))
        plt.plot(t_compressed, compressed)
        plt.title('Compressed Data')
        plt.show()
    else:
        ax.plot(t_compressed, compressed)
