"""
stm.py



Tools to handle CR3BP STM

Functions
---------
computeStm                        : Compute STM
computeStmParallel                : Compute Multiple STMs in Parallel
computeCauchyGreenTensor          : Compute CGT
computeStabilityIndex             : Compute Stability Index
computeStableUnstableEigenvectors : Compute Stable/Unstable Eigenvectors

"""
import numpy as np 

def computeStm(ode, x0, tf):
    """
    Compute State Transition Matrix

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        Generic ASSET ODE
    x0 : list or np.array
        Initial State (must include initial time)
    tf : float
        Final integration time

    Returns
    -------
    stm : np.array
        State Transition Matrix
        
    """
    integrator = ode.integrator('DOPRI87',1e-4)

    _, jac = integrator.integrate_stm(x0,tf)

    stm = jac[0:6, 0:6]

    return stm

def computeStmParallel(ode, x0s, tfs, nthreads:int=8):
    """
    Compute Multiple State Transition Matrices in Parallel

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        Generic ASSET ODE
    x0s : list of list or np.array
        Initial States (must include initial time)
    tfs : list of float
        Final integration times
    nthreads : int
        Number of threads to use in parallel inegration
        
    Returns
    -------
    stms : np.array
        State Transition Matrix

    """
    integrator = ode.integrator('DOPRI87',1e-4)

    results = integrator.integrate_stm_parallel(x0s,tfs,nthreads)

    stms = []

    for _, jac in results:

        stm = jac[0:6, 0:6]
        stms.append(stm)

    return stms

def computeStmSvd(stm):
    """
    Compute STM Singular Value Decomposition
    
    Parameters
    ----------
    stm : np.array
        State Transition Matrix
    
    Returns
    -------
    u,s,v : tuple
        Singular Value Decomposition Matrices

    """

    u,s,v = np.linalg.svd(stm)

    return u,s,v

def computeCauchyGreenTensor(stm):
    """
    Compute the Cauchy Green Strain Tensor from a State Transition Matrix

    Parameters
    ----------
    stm : np.array
        State Transition Matrix

    Returns
    -------
    cgt : np.array
        Cauchy Green Strain Tensor

    """

    cgt = stm.T @ stm

    return cgt

def computeStabilityIndex(stm):
    """
    Compute a Periodic Orbit Stability Index

    Parameters
    ----------
    stm : np.array
        State Transition Matrix (Note: Must be the Monodromy matrix)

    Returns
    -------
    stab_idx : float
        Stability Index

    """

    eig_val, _ = np.linalg.eig(stm)

    stab_idx = (1.0/2.0) * (np.max(np.abs(eig_val)) + 1.0/np.max(np.abs(eig_val)))

    return stab_idx

def computeStableUnstableEigenvectors(stm):
    """
    Compute the stable and unstable eigenvectors of a STM

    Parameters
    ----------
    stm : np.array
        State Transition Matrix

    Returns
    -------
    stable_vec, unstable_vec : tuple
        Stable and Unstable Eigenvectors

    """

    eig_val, eig_vec = np.linalg.eig(stm)

    eig_vec = eig_vec.T

    idxs = list(range(0,6))

    idxs.sort(key=lambda x: np.abs(eig_val[x]))

    stable_vec   = np.real(eig_vec[idxs[0]])
    unstable_vec = np.real(eig_vec[idxs[-1]])

    return stable_vec, unstable_vec 