def RayleighDamping(omega1, omega2, ksi1, ksi2):
    """
    Return the coefficients of the 2-point Rayleigh damping.
    C = alpha * M + beta * K

    Parameters:
    omega1, omega2 : float
        Natural frequencies (rad/s)
    ksi1, ksi2 : float
        Damping ratios

    Returns:
    alpha, beta : float
        Rayleigh damping coefficients
    """
    beta = 2 * (omega2 * ksi2 - omega1 * ksi1) / (omega2**2 - omega1**2)
    alpha = 2 * omega1 * ksi1 - beta * omega1**2
    return alpha, beta
