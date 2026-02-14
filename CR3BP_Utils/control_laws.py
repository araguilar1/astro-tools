"""
control_laws.py



Generic Control Laws used in CR3BP (typically low-thrust enabled)

Functions
---------
basicFuncVnb    : Defines VNB Basis Vectors for Thrust Directions
progradeLawVnb  : Thrusts in a VNB (V) direction
progradeLawVDir : Thrusts parallel to current velocity vector

"""
import asset_asrl.VectorFunctions as vf
from asset_asrl.VectorFunctions import Arguments as Args

def basicFuncVnb():
    """
    Generate VNB Basis Vectors at current state vector

    Parameters
    ----------
    None

    Return
    ------
    vnb_basis : asset_asrl.VectorFunction
        VNB Basis Vectors
    """

    R,V = Args(6).tolist([(0,3),(3,3)])
    Vhat = V.normalized()
    Nhat = R.cross(V).normalized()
    Bhat = V.cross(Nhat).normalized()

    vnb_basis = vf.stack(Vhat,Nhat,Bhat)

    return vnb_basis


def progradeLawVnb(throttle):
    """
    Control Law to Thrust in VNB Basis vector V

    Parameters
    ----------
    throttle : float
        Ratio of throttle to use
    
    Returns
    -------
    U : asset_asrl.VectorFunction
        Control
    """
    
    RV = Args(6).vf()
    VNBBasis = basicFuncVnb()(RV)
    U = vf.RowMatrix(VNBBasis,3,3)*RV.tail(3).normalized()*throttle
    
    return U


def progradeLawVDir(throttle):
    """
    Control Law to Thrust parallel to current Velocity Vector

    Parameters
    ----------
    throttle : float
        Ratio of throttle to use
    
    Returns
    -------
    U : asset_asrl.VectorFunction
        Control

    """
    v = Args(3)
    U = v.normalized()*throttle

    return U