"""
constraints.py



Common CR3BP Inequality/Equality Constraints to be used with ASSET's PSIOPT optimizer

Functions
---------
posConTable   : Position Equality Constraint using LGLInterpTable
posConFixed   : Position Equality Constraint using Fixed State
rendConTable  : Rendezvous Equality Constraint using LGLInterpTable
p2DistCon     : P2 distance Inequality Constraint
p1DistCon     : P1 distance Inequality Constraint

"""
from asset_asrl.VectorFunctions import Arguments as Args
import asset_asrl.OptimalControl as oc
import numpy as np 

def posConTable(table):
    """
    Equality Constraint that enforces position of phase must lie on Interpolation Table

    Parameters
    ----------
    table : asset_asrl.OptimalControl.LGLInterpTable
        LGLInterpTable of State Vectors, e.g. depature orbit, target orbit, or any ASSET phase

    Returns
    -------
    con : asset_asrl.VectorFunction
        Equality Constraint
    """

    Rt = Args(4)
    R  = Rt.head(3)
    t  = Rt[3]

    func = oc.InterpFunction(table, range(0,3)).vf()

    con = R - func(t)

    return con


def posConFixed(target):
    """
    Equality Constraint that enforces position of phase must be equal to a fixed state

    Parameters
    ----------
    target : np.array or list
        Fixed State position must be equal to

    Returns
    -------
    con : asset_asrl.VectorFunction
        Equality Constraint
    """

    R = Args(3)

    con = R - target

    return con


def rendConTable(table):
    """
    Equality Constraint that enforces position and velocity of phase must lie on Interpolation Table

    Parameters
    ----------
    table : asset_asrl.OptimalControl.LGLInterpTable
        LGLInterpTable of State Vectors, e.g. depature orbit, target orbit, or any ASSET phase

    Returns
    -------
    con : asset_asrl.VectorFunction
        Equality Constraint
    """
    Xt = Args(7)
    X = Xt.head(6)
    t = Xt[6]

    func = oc.InterpFunction(table, range(0,6))

    con = X - func(t)

    return con


def p2DistCon(ode, nd_dist):
    """
    Inequality Constraint that enforces state must maintain certain distance from P2
    
    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        General ASSET CR3BP ODE
    nd_dist : float
        Non-Dimensional keep out distance

    Returns
    -------
    con : assset_asrl.VectorFunction
        Inequlity Constraint
    """
    p2_loc = np.array([1-ode.mu,0.0,0.0])
    R = Args(3)

    con = nd_dist - (R - p2_loc).norm()

    return con


def p1DistCon(ode, nd_dist):
    """
    Inequality Constraint that enforces state must maintain certain distance from P1
    
    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        General ASSET CR3BP ODE
    nd_dist : float
        Non-Dimensional keep out distance

    Returns
    -------
    con : assset_asrl.VectorFunction
        Inequlity Constraint
    """
    p1_loc = np.array([-ode.mu,0.0,0.0])
    R = Args(3)

    con = nd_dist - (R - p1_loc).norm()

    return con
