"""
manifolds.py



Compute CR3BP Periodic Orbit Manifolds

Functions
---------   
getManifoldsDrDv                 : returns manifolds computed by perturbing position and velocity along eigenvector
getManifoldsDr                   : returns manifolds computed by perturbing position along eigenvector
getManifoldInitialConditionsDrDv : returns manifold initial conditions computed by perturbing position and velocity along eigenvector
getManifoldInitialConditionsDr   : returns manifold initial conditions computed by perturbing position along eigenvector
"""
from .stm import computeStmParallel, computeStableUnstableEigenvectors
from ..general_utils import M_TO_KM
import numpy as np

def getManifoldsDrDv(ode, orbit_in:list, dx:float, dt:float, events:list=[], nman:int=50, stable:bool=True, ncores:int=8):
    """
    Compute Periodic Orbit Manifolds by perturbing both position and velocity of the states along the orbit

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        Generic ASSET ODE
    orbit_in : list of array
        Periodic Orbit
    dx : float
        Step size off manifolds
    dt : float
        Manifold propagation time
    events : list
        List of ASSET integrator event conditions
    nman : int
        Number of manifolds to compute
    stable : bool
        Compute stable manifolds (True) or Unstable Manifolds (False)
    ncores : int
        Number of cores to use in parallel integration (default 8)

    Returns
    -------
    manifolds : list
        Periodic Orbit Manifolds

    """
    integrator = ode.integrator('DOPRI87',1e-4)
    integrator.setAbsTol(1.0e-13)

    period = orbit_in[-1][6]
    orbit = integrator.integrate_dense(orbit_in[0],period,nman)

    tfs = [o[6]+period for o in orbit]
    stms = computeStmParallel(ode, orbit, tfs, ncores)

    ics = [] # initial conditions of perturbed points

    manifolds = []

    for i, stm in enumerate(stms):

        stable_vec,unstable_vec = computeStableUnstableEigenvectors(stm)

        if stable: vec = stable_vec
        else: vec = unstable_vec

        x_plus = np.copy(orbit[i])
        x_plus[-1] = 0.0
        x_plus[:6] += vec[:6]/np.linalg.norm(vec[:3])*dx

        x_minus = np.copy(orbit[i])
        x_minus[-1] = 0.0
        x_minus[:6] -= vec[:6]/np.linalg.norm(vec[:3])*dx

        ics.append(x_plus)
        ics.append(x_minus)

    if stable: dt = -dt

    ts = [ic[6] + dt for ic in ics]

    if events:
        results = integrator.integrate_dense_parallel(ics,ts,events,ncores)
        for result in results:
            traj, event_locs = result
            if len(event_locs[0]) > 0:
                traj.pop()
                traj.append(event_locs[0][1])

            manifolds.append(traj)

    else:
        manifolds = integrator.integrate_dense_parallel(ics,ts,ncores)

    return manifolds


def getManifoldsDr(ode, orbit_in:list, dx:float, dt:float, events:list=[], nman:int=50, stable:bool=True, ncores:int=8):
    """
    Compute Periodic Orbit Manifolds by only perturbing the position of the states along the orbit

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        Generic ASSET ODE
    orbit_in : list of array
        Periodic Orbit
    dx : float
        Step size off manifolds
    dt : float
        Manifold propagation time
    events : list
        List of ASSET integrator event conditions
    nman : int
        Number of manifolds to compute
    stable : bool
        Compute stable manifolds (True) or Unstable Manifolds (False)
    ncores : int
        Number of cores to use in parallel integration (default 8)

    Returns
    -------
    manifolds : list
        Periodic Orbit Manifolds

    """
    integrator = ode.integrator('DOPRI87',1e-4)
    integrator.setAbsTol(1.0e-13)

    period = orbit_in[-1][6]
    orbit = integrator.integrate_dense(orbit_in[0],period,nman)

    tfs = [o[6]+period for o in orbit]
    stms = computeStmParallel(ode, orbit, tfs, ncores)

    ics = [] # initial conditions of perturbed points

    manifolds = []

    for i, stm in enumerate(stms):

        stable_vec,unstable_vec = computeStableUnstableEigenvectors(stm)

        if stable: vec = stable_vec
        else: vec = unstable_vec

        x_plus = np.copy(orbit[i])
        x_plus[-1] = 0.0
        x_plus[:3] += vec[:3]/np.linalg.norm(vec[:3])*dx

        x_minus = np.copy(orbit[i])
        x_minus[-1] = 0.0
        x_minus[:3] -= vec[:3]/np.linalg.norm(vec[:3])*dx


        ics.append(x_plus)
        ics.append(x_minus)

    if stable: dt = -dt

    ts = [ic[6] + dt for ic in ics]

    if events:
        results = integrator.integrate_dense_parallel(ics,ts,events,ncores)
        for result in results:
            traj, event_locs = result
            if len(event_locs[0]) > 0:
                traj.pop()
                traj.append(event_locs[0][1])

            manifolds.append(traj)

    else:
        manifolds = integrator.integrate_dense_parallel(ics,ts,ncores)

    return manifolds


def getManifoldInitialConditionsDrDv(ode, orbit_in:list, dx:float, dt:float, nman:int=50, stable:bool=True, ncores:int=8):
    """
    Compute Periodic Orbit Manifolds Initial Conditions only by perturbing both position and velocity of the states along the orbit

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        Generic ASSET ODE
    orbit_in : list of array
        Periodic Orbit
    dx : float
        Step size off manifolds
    dt : float
        Manifold propagation time
    nman : int
        Number of manifolds to compute
    stable : bool
        Compute stable manifolds (True) or Unstable Manifolds (False)
    ncores : int
        Number of cores to use in parallel integration (default 8)

    Returns
    -------
    ics, ts : tuple of list
        Periodic Orbit Manifolds Initial Conditions and Final Integration Times

    """
    print(f'Perturbation Magnitude: {dx*ode.lstar*M_TO_KM} km | {dx*ode.vstar} m/s')

    integrator = ode.integrator('DOPRI87',1e-4)
    integrator.setAbsTol(1.0e-13)

    period = orbit_in[-1][6]
    orbit = integrator.integrate_dense(orbit_in[0],period,nman)

    tfs = [o[6]+period for o in orbit]
    stms = computeStmParallel(ode, orbit, tfs, ncores)

    ics = [] # initial conditions of perturbed points

    for i, stm in enumerate(stms):

        stable_vec,unstable_vec = computeStableUnstableEigenvectors(stm)

        if stable: vec = stable_vec
        else: vec = unstable_vec

        x_plus = np.copy(orbit[i])
        x_plus[-1] = 0.0
        x_plus[:6] += vec[:6]/np.linalg.norm(vec[:3])*dx

        x_minus = np.copy(orbit[i])
        x_minus[-1] = 0.0
        x_minus[:6] -= vec[:6]/np.linalg.norm(vec[:3])*dx

        ics.append(x_plus)
        ics.append(x_minus)

    if stable: dt = -dt

    ts = [ic[6] + dt for ic in ics]

    return ics, ts


def getManifoldInitialConditionsDr(ode, orbit_in:list, dx:float, dt:float, nman:int=50, stable:bool=True, ncores:int=8):
    """
    Compute Periodic Orbit Manifolds Initial Conditions only by only perturbing position of the states along the orbit

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        Generic ASSET ODE
    orbit_in : list of array
        Periodic Orbit
    dx : float
        Step size off manifolds
    dt : float
        Manifold propagation time
    nman : int
        Number of manifolds to compute
    stable : bool
        Compute stable manifolds (True) or Unstable Manifolds (False)
    ncores : int
        Number of cores to use in parallel integration (default 8)

    Returns
    -------
    ics, ts : tuple of list
        Periodic Orbit Manifolds Initial Conditions and Final Integration Times
    """
    print(f'Perturbation Magnitude: {dx*ode.lstar*M_TO_KM} km')

    integrator = ode.integrator('DOPRI87',1e-4)
    integrator.setAbsTol(1.0e-13)

    period = orbit_in[-1][6]
    orbit = integrator.integrate_dense(orbit_in[0],period,nman)

    tfs = [o[6]+period for o in orbit]
    stms = computeStmParallel(ode, orbit, tfs, ncores)

    ics = [] # initial conditions of perturbed points

    for i, stm in enumerate(stms):

        stable_vec,unstable_vec = computeStableUnstableEigenvectors(stm)

        if stable: vec = stable_vec
        else: vec = unstable_vec

        x_plus = np.copy(orbit[i])
        x_plus[-1] = 0.0
        x_plus[:3] += vec[:3]/np.linalg.norm(vec[:3])*dx

        x_minus = np.copy(orbit[i])
        x_minus[-1] = 0.0
        x_minus[:3] -= vec[:3]/np.linalg.norm(vec[:3])*dx

        ics.append(x_plus)
        ics.append(x_minus)

    if stable: dt = -dt

    ts = [ic[6] + dt for ic in ics]

    return ics, ts