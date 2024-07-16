import numpy as np

def newmark(m, xi, f, p, dt, Nstep, a=0.25, d=0.5, dep0=0, vit0=0):
    
    omega = 2 * np.pi * f
    k = m * omega ** 2
    c = 2 * m * omega * xi

    dep = np.zeros(Nstep+1)
    dep[0] = dep0
    vit = np.zeros(Nstep+1)
    vit[0] = vit0
    acc = np.zeros(Nstep+1)
    qs = np.zeros(Nstep+1)
    akf = m / a / dt ** 2 + d / a / dt * c + k

    for i in range(1, Nstep):
        dep[i] = (p[i] + m * (1 / a / dt ** 2 * dep[i-1] + 1 / a / dt * vit[i-1] + (1 / 2 / a - 1) * acc[i-1])
                + c * (d / a / dt * dep[i-1] + (d / a - 1) * vit[i-1] + dt / 2 * (d / a - 2) * acc[i-1])) / akf
        vit[i] = d / a / dt * (dep[i] - dep[i-1]) + (1 - d / a) * vit[i-1] + dt * (1 - d / 2 / a) * acc[i-1]
        acc[i] = 1 / a / dt ** 2 * (dep[i] - dep[i-1]) - 1 / a / dt * vit[i-1] - (1 / 2 / a - 1) * acc[i-1]
        qs[i] = p[i] / k

    time = np.linspace(0.,Nstep*dt,Nstep+1)
    return time, dep, vit, acc, qs


def NewmarkMDDL(m, k, c, p, dt, Nstep, a=0.25, d=0.5, depin=None, vitin=None):
    nddl = len(k)

    dep = np.zeros((nddl, Nstep))
    vit = np.zeros((nddl, Nstep))
    acc = np.zeros((nddl, Nstep))
    qs = np.zeros((nddl, Nstep))

    if depin is not None and vitin is not None:
        dep[:, 0] = depin
        vit[:, 0] = vitin

    akf = m / a / dt ** 2 + d / a / dt * c + k
    invakf = np.linalg.inv(akf)
    invk = np.linalg.inv(k)

    for i in range(1, Nstep):
        dep[:, i] = invakf @ (p[:, i] + m @ (1 / a / dt ** 2 * dep[:, i - 1] + 1 / a / dt * vit[:, i - 1] + (1 / 2 / a - 1) * acc[:, i - 1])
                             + c @ (d / a / dt * dep[:, i - 1] + (d / a - 1) * vit[:, i - 1] + dt / 2 * (d / a - 2) * acc[:, i - 1]))
        vit[:, i] = d / a / dt * (dep[:, i] - dep[:, i - 1]) + (1 - d / a) * vit[:, i - 1] + dt * (1 - d / 2 / a) * acc[:, i - 1]
        acc[:, i] = 1 / a / dt ** 2 * (dep[:, i] - dep[:, i - 1]) - 1 / a / dt * vit[:, i - 1] - (1 / 2 / a - 1) * acc[:, i - 1]
        qs[:, i] = invk @ p[:, i]

    time = np.linspace(0.,Nstep*dt,Nstep)
    return time, dep, vit, acc, qs




def diff_centrale1(m, c, k, dt, p, N, dep0=0, vit0=0):
    dep = np.zeros(N)
    vit = np.zeros(N)
    dep[0] = dep0
    vit[0] = vit0

    for i in range(N-1):
        dep[i+1] = dep[i] + dt*vit[i] + dt**2/2/m*(p[i] - c*vit[i] - k*dep[i])
        vit[i+1] = 2*(dep[i+1] - dep[i])/dt - vit[i]

    time = np.linspace(0.,N*dt,N)
    return time, dep, vit


def diff_centrale2(m, c, k, dt, p, N, dep0=0, vit0=0):
    if N <= 1:
        raise ValueError("N must be greater than 1 for this function.")

    dep = np.zeros(N+1)
    dep_m1 = dep0 - dt*vit0

    dep[0] = dep_m1
    dep[1] = dep0

    for i in range(1, N):
        dep[i+1] = (p[i] + (2*m/dt**2-k)*dep[i] - (m/dt**2-c/2/dt)*dep[i-1]) / (m/dt**2+c/2/dt)

    time = np.linspace(0.,N*dt,N)
    return time, dep[1:]

def CHNL(param):
    """ Solve nonlineat dynamic problem using Chung-Hulbert method
    param : dictionary containing the following
    - rhoInf : spectral radius at infinity
    - num : [nCp, h, aTol, iterMax] = [number of time steps, duration of simulation,
         time step, tolerance on acceleration, max number of iterations per time step]]
    - M : mass matrix
    - IC : initial conditions [u0, v0] - size = 2*nDof
    - fExt : external force function - format fExt = fExt(t, param)
    - fInt : internal force function - format: fInt, K, C = fInt(u, v, param)
    """
    # Set Chung-Hulbert parameters
    rhoInf = param['rhoInf']
    a_m = (2*rhoInf - 1) / (rhoInf + 1)
    a_f = rhoInf / (rhoInf + 1)
    beta = 0.25 * (1 - a_m + a_f)**2
    gamma = 0.5 - a_m + a_f

    # Get integration parameters
    nCp = param['num'][0]
    h = param['num'][1] / nCp
    aTol = param['num'][2]
    iterMax = param['num'][3]

    # Get model parameters
    M = param['M']
    nDof = M.shape[0]
    ind = np.arange(0, nDof)
    uOld = param['IC'][ind]
    vOld = param['IC'][ind + nDof]
    fExt = param['fExt']
    fInt = param['fInt']

    # Initialize storage arrays
    TIME = np.zeros(nCp + 1)
    STATE = np.empty((3 * nDof, nCp + 1))
    FINT = np.empty((nDof, nCp + 1))
    KE = np.empty(nCp + 1)
    WINT = np.zeros(nCp + 1)

    # Initialize local variables
    fOld, KOld,_ = fInt(uOld, vOld, param)
    aOld = np.linalg.solve(M, (fExt(0, param) - fOld))
    old = np.concatenate((uOld, vOld, aOld))
    tOld = 0
    STATE[:, 0] = old
    FINT[:, 0] = fOld
    KE[0] = 0.5 * np.dot(vOld, np.dot(M, vOld))
    incr = 2

    while incr <= nCp:
        # Increment time
        tNew = tOld + h

        # Compute predictor
        uPred = uOld + h * vOld + h**2 * (0.5 - beta) * aOld
        vPred = vOld + h * (1 - gamma) * aOld
        aNew = aOld.copy()

        # Enter iterative loop
        fNew, KNew, CNew = fInt(uPred, vPred, param)
        rNew = np.dot(M, ((1 - a_m) * aNew + a_m * aOld)) + \
               (1 - a_f) * (fNew - fExt(tNew, param)) + a_f * (fOld - fExt(tOld, param))
        iter = 1
        err = np.linalg.norm(rNew)
        if param.get('debug', 0) > 0:
            print('INCR:', incr)
            print(f'   ITER: {iter} - RES: {err:.3e}')
        while (iter <= iterMax and err > aTol) or iter == 1:
            JNew = JACOBIAN(h, KNew, CNew, M, beta, gamma, a_m, a_f)
            corr = -np.linalg.solve(JNew, rNew)
            aNew += corr
            uNew = uPred + h**2 * beta * aNew
            vNew = vPred + h * gamma * aNew
            fNew, KNew, CNew = fInt(uNew, vNew, param)
            rNew = np.dot(M, ((1 - a_m) * aNew + a_m * aOld)) + \
                   (1 - a_f) * (fNew - fExt(tNew, param)) + a_f * (fOld - fExt(tOld, param))
            err = np.linalg.norm(rNew)
            iter += 1
            if param.get('debug', 0) > 0:
                print(f'   ITER: {iter} - RES: {err:.3e}')
        assert np.linalg.norm(rNew) < aTol, 'CHNL:NR:failed to converge'
        TIME[incr] = tNew
        STATE[:, incr] = np.concatenate((uNew, vNew, aNew))
        FINT[:, incr] = fNew
        KE[incr] = 0.5 * np.dot(vNew, np.dot(M, vNew))
        WINT[incr] = WINT[incr - 1] + fIntWork(h, fOld, fNew, uOld, uNew, vOld, vNew, KOld, KNew)
        tOld = tNew
        uOld = uNew
        vOld = vNew
        aOld = aNew
        fOld = fNew
        incr += 1

    sol = {'time': TIME,
           'disp': STATE[ind, :],
           'velo': STATE[ind + nDof, :],
           'accel': STATE[ind + 2 * nDof, :],
           'fInt': FINT,
           'KE': KE,
           'WINT': WINT}

    return sol


def JACOBIAN(h, K, C, M, beta, gamma, a_m, a_f):
    J = (1 - a_m) * M + (1 - a_f) * (beta * h**2 * K + gamma * h * C)
    return J


def fIntWork(h, fOld, fNew, uOld, uNew, vOld, vNew, KOld, KNew):
    du = uNew - uOld
    dwInt = 0.5 * np.dot(du, fOld + fNew) - \
            h / 12 * (np.dot(vNew - vOld, fNew - fOld) - np.dot(du, np.dot(KNew, vNew) - np.dot(KOld, vOld)))
    return dwInt
