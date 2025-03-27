import numpy as np
from scipy.interpolate import interp1d
from scipy.linalg import eigh
from scipy.integrate import trapz
import matplotlib.pyplot as plt


def generate_samples(hdl_s, t, f_calcul):
    # Ensure t has an even length
    if len(t) % 2:
        t = np.append(t, 2 * t[-1] - t[-2])
    n = len(t)

    dt = t[1] - t[0]
    df = 1 / t[-1]

    t = np.arange(0, n) * dt
    f = np.arange(0, n) * df

    if f_calcul[-1] < f[-1]:
        f_calcul = np.append(f_calcul, f[-1])

    # Modal decomposition of the desired power spectral density matrix
    if is_vector_input_supported(hdl_s):
        S = hdl_s(2 * np.pi * f_calcul)

        V = np.zeros((S.shape[0], S.shape[1], len(f_calcul)))
        D = np.zeros((S.shape[0], len(f_calcul)))

        for ifr, fc in enumerate(f_calcul):
            eigvals, eigvecs = eigh(S[:, :, ifr])
            order = np.argsort(eigvals)[::-1]
            eigvals = eigvals[order]
            eigvecs = eigvecs[:, order]

            if ifr > 0:
                for i in range(eigvecs.shape[1]):
                    if np.dot(eigvecs[:, i].T, V[:, i, ifr - 1]) < 0:
                        eigvecs[:, i] *= -1

            V[:, :, ifr] = eigvecs
            D[:, ifr] = eigvals
    else:
        V, D = [], []
        for ifr, fc in enumerate(f_calcul):
            S = hdl_s(2 * np.pi * fc)
            eigvals, eigvecs = eigh(S)
            order = np.argsort(eigvals)[::-1]
            eigvals = eigvals[order]
            eigvecs = eigvecs[:, order]

            V.append(eigvecs)
            D.append(eigvals)

    cov_target = trapz(4 * np.pi * f_calcul, S, axis=2)

    Nhist = S.shape[0]
    NModes = Nhist

    # Generate modal processes
    Nu = np.zeros((n, NModes), dtype=np.complex_)

    for m in range(NModes):
        Vrm = interp1d(f_calcul, D[m, :], kind="linear", fill_value="extrapolate")(f)
        theta = 2 * np.pi * np.random.rand(n // 2 - 1)
        c = 4 * np.pi / (2 * dt)

        Nu[0, m] = 0
        Nu[1:n // 2, m] = c * np.sqrt(n * dt / (2 * np.pi) * Vrm[1:n // 2]) * np.exp(1j * theta)
        Nu[n // 2, m] = c * np.sqrt(n * dt / (2 * np.pi) * Vrm[n // 2])
        Nu[n // 2 + 1:, m] = np.conj(Nu[n // 2 - 1:0:-1, m])

    # Generate nodal processes
    ZMAT = np.zeros((n, NModes), dtype=np.complex_)
    X = np.zeros((n, Nhist), dtype=np.complex_)

    for h in range(Nhist):
        for m in range(NModes):
            V_interp = interp1d(f_calcul, V[h, m, :], kind="linear", fill_value="extrapolate")
            ZMAT[1:n // 2, m] = V_interp(f[1:n // 2])
            ZMAT[n // 2, m] = np.real(V_interp(f[n // 2]))
            ZMAT[n // 2 + 1:, m] = np.conj(ZMAT[n // 2 - 1:0:-1, m])

        X[:, h] = np.sum(ZMAT * Nu, axis=1)

    samples = np.zeros((n, Nhist))
    for h in range(Nhist):
        samples[:, h] = np.fft.ifft(X[:, h]).real

    return t, samples


def is_vector_input_supported(func):
    try:
        test_vector = np.array([1, 2, 3, 4])
        func(test_vector)
        return True
    except:
        return False
