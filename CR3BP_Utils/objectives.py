"""
objectives.py



Common CR3BP Objective Functions to be used with ASSET's PSIOPT optimizer

Functions
---------
dvObjTable   : Minimize delta-V with Interpolation Table
dvObjFixed   : Minimize delta-V with Fixed State
fwdLinkDvObj : Minimize delta-V in a forward link objective between two phases

Classes
-------
massObjLt  : Low-Thrust Mass Optimal Objective
powerObjLt : Low-Thrust Power Optimal Objective
"""
import asset_asrl.VectorFunctions as vf
from asset_asrl.VectorFunctions import Arguments as Args
import asset_asrl.OptimalControl as oc


def dvObjTable(table):
    """
    Minimize delta-V with state from Interpolation Table

    Parameters
    ----------
    table : asset_asrl.OptimalControl.LGLInterpolationTable
        LGLInterpTable of State Vectors, e.g. depature orbit, target orbit, or any ASSET phase

    Returns
    -------
    obj : asset_asrl.VectorFunction
        Objective Function
        
    """
    Vt =  Args(4)
    V = Vt.head(3)
    t = Vt[3]

    func = oc.InterpFunction(table, range(3,6)).vf()

    obj = (V - func(t)).squared_norm()

    return obj


def dvObjFixed(target):
    """
    Minimize delta-V with from a fixed state

    Parameters
    ----------
    target : np.array or list
        Fixed State to mimize delta-V from

    Returns
    -------
    obj : asset_asrl.VectorFunction
        Objective Function

    """

    V = Args(3)

    obj = (V - target).squared_norm()

    return obj


def fwdLinkDvObj():
    """
    Minimize delta-V between two ASSET phases

    Parameters
    ----------
    None

    Returns
    -------
    obj : asset_asrl.VectorFunction
        Objective Function

    """
    v1, v2 = Args(6).tolist([(0,3),(3,3)])
    obj = (v2-v1).squared_norm()
    return obj


class massObjLt(vf.ScalarFunction):
    
    def __init__(self, scale):
        """
        Low-Thrust Mass Optimal Objective Function

        Parameters
        ----------
        scale : float
            Scale of control vector norm
        
        Returns
        -------
        asset_asrl.VectorFunctions.ScalarFunction

        """

        u = Args(3)
        super().__init__(u.norm() * scale)


class powerObjLt(vf.ScalarFunction):
    
    def __init__(self, scale):
        """
        Low-Thrust Power Optimal Objective Function

        Parameters
        ----------
        scale : float
            Scale of control vector norm squared
        
        Returns
        -------
        asset_asrl.VectorFunctions.ScalarFunction

        """

        u = Args(3)
        super().__init__(u.norm().squared() * scale)


class momentumIntegralObj(vf.ScalarFunction):

        def __init__(self):
            """
            Momentum Integral Objective Function

            Parameters
            ----------
            None

            Returns
            -------
            asset_asrl.VectorFunctions.ScalarFunction
        
            """
             
            RV = Args(6)

            x,y,z,vx,vy,vz = RV.tolist()

            super().__init__(x*vx + y*vy + z*vz)
        
