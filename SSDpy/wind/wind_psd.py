def wind_psd(f, param, spectre):
    """
    Compute spectral densities based on the provided parameters and spectrum type.

    Parameters:
    - f : array-like
        Frequencies.
    - param : list or array
        Parameters for spectrum definitions.
    - spectre : int
        Spectrum type:
          1 : Von Karman (parameters: sigma, U, L)
          2 : Davenport (parameters: sigma, U, L=1200 fixed)
          3 : Eurocode longitudinal (parameters: sigma, U, L)
          4 : Transverse Von Karman (parameters: sigma, U, L)

    Returns:
    - S : array-like
        Spectral density values.
    """
    if spectre == 1:  # Von Karman
        sig = param[0]
        U = param[1]
        L = param[2]
        S = 4 * L / U * sig**2 / (1 + 70.7 * (f * L / U)**2)**(5 / 6)

    elif spectre == 2:  # Davenport
        sig = param[0]
        U = param[1]
        L = 1200
        S = (2 / 3) * L / U * sig**2 * (f * L / U) / (1 + (f * L / U)**2)**(4 / 3)

    elif spectre == 3:  # Eurocode longitudinal component
        sig = param[0]
        U = param[1]
        L = param[2]
        S = sig**2 * 6.8 * (L / U) / (1 + 10.2 * (f * L / U))**(5 / 3)

    elif spectre == 4:  # Transverse Von Karman
        sig = param[0]
        U = param[1]
        L = param[2]
        S = 4 * L / U * sig**2 * (1 + 755.2 * (f * L / U)**2) / (1 + 283.2 * (f * L / U)**2)**(11 / 6)

    else:
        raise ValueError("Invalid spectre type. Choose between 1, 2, 3, or 4.")

    return S
