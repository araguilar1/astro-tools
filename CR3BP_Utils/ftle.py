"""
ftle.py



Functions to Compute Finite Time Lyapunov Exponents

Functions
---------
computeLocalLyapunovExponents : Compute LLEs

"""
import numpy as np
from .stm import computeStmParallel, computeCauchyGreenTensor

def computeLocalLyapunovExponents(ode, traj:list, time_horizon:list, nt:int=12, nthreads:int=8) -> np.ndarray:
    """
    Compute Local Lyapunov Exponents

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        ASSET ODE
    traj : list of array
        Trajectory to evaluate LLEs
    time_horizon : list
        Min and Max time horizon over which to evaluate LLEs
    nt : int
        Number of evenly space intervals between min and max time horizons
    nthreads : int
        Number of threads used in parallel integration

    Returns
    -------
    lle : np.ndarray
        Local Lyapunov Exponents
    """

    npoints = len(traj)

    lle = np.empty((npoints, npoints))

    Th = np.linspace(time_horizon[0], time_horizon[1], nt)

    for th in Th:
        tfs = [t[6] + th for t in traj]

        stms = computeStmParallel(ode, traj, tfs, nthreads)

        local_lypaunov_exp = np.zeros((npoints))

        for i in range(npoints):
            for j, stm in enumerate(stms):

                cgt = computeCauchyGreenTensor(stm)

                eig_val, _ = np.linalg.eig(cgt)

                local_lypaunov_exp[j] = np.log(np.sqrt(np.max(np.abs(eig_val))))

            lle[i,:] = local_lypaunov_exp[:]

    return lle