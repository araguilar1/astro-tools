"""
propagator.py

ASSET-based orbit propagation for the CR3BP.

Functions
---------
create_ode        : Create a CR3BP ODE instance (cached)
propagate_orbit   : Propagate a single orbit using ASSET DOPRI87
propagate_family  : Propagate a list of orbits
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from Models.CR3BP import CR3BP

_ode_cache = {}


def create_ode(mu, lstar, tstar):
    """
    Create or retrieve a cached CR3BP ODE instance.

    Parameters
    ----------
    mu : float
        Mass ratio
    lstar : float
        Characteristic length (km)
    tstar : float
        Characteristic time (s)

    Returns
    -------
    CR3BP
        ASSET ODE instance
    """
    key = (mu, lstar, tstar)
    if key not in _ode_cache:
        _ode_cache[key] = CR3BP(mu, lstar, tstar)
    return _ode_cache[key]


def propagate_orbit(ode, state_6, period, n_points=300):
    """
    Propagate a single orbit using ASSET's DOPRI87 integrator.

    Parameters
    ----------
    ode : CR3BP
        ASSET ODE instance
    state_6 : list of float
        [x, y, z, vx, vy, vz]
    period : float
        Full period (non-dimensional time)
    n_points : int
        Number of output points

    Returns
    -------
    list of np.ndarray
        Each element is a 7-element state [x,y,z,vx,vy,vz,t]
    """
    integrator = ode.integrator("DOPRI87", 1e-4)
    ic_7 = list(state_6) + [0.0]
    traj = integrator.integrate_dense(ic_7, period, n_points)
    return traj


def propagate_family(ode, orbit_dicts, n_points=300):
    """
    Propagate a list of orbits from catalog rows.

    Parameters
    ----------
    ode : CR3BP
        ASSET ODE instance
    orbit_dicts : list of dict
        Catalog rows (from data_loader.filter_orbits)
    n_points : int
        Points per orbit

    Returns
    -------
    list of list
        Each inner list is a trajectory (list of 7-element states)
    """
    from services.data_loader import orbit_row_to_ic

    trajectories = []
    for row in orbit_dicts:
        state_6, period, _, _ = orbit_row_to_ic(row)
        traj = propagate_orbit(ode, state_6, period, n_points)
        trajectories.append(traj)
    return trajectories
