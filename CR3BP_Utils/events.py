"""
events.py



Common CR3BP Integrator Events to be used with ASSET

Functions
---------
apseEvent       : event when spacecraft is at an apsis
p2DistnaceEvent : event when spacecraft is certain distance from P2
p2CullEvent     : event wehen spacecraft reaches specified distance or dips below specified altitude from P2 
xAxisCrossEvent : event when spacecraft crosses x-axis
xCrossEVent     : event when spacecraft hits specified x-position
yCrossEvent     : event when spacecraft hits specified y-position
yAxisCrossEvent : event when spacecraft crosses y-axis
zCrossEvent     : event when spacecraft hits specified z-position
xCrossP2Event   : event when spacecraft crosses P2 x position
xCrossP1Event   : event when spacecraft crosses P1 x position

"""
from asset_asrl.VectorFunctions import Arguments as Args
import numpy as np


def apseEvent(st_flag=False):
    """
    Integrator event for when P3 (spacecraft) is at periapsis or apoapsis
    1 : periapsis | -1: apoapsis
    
    Parameters
    ----------
    st_flag : bool 
        Boolean indicating whether or not the event is to be used with a Sundman Transformed ODE

    Returns
    -------
    rdotv : asset_asrl.VectorFunction 
        Dot product of position and velocity
    """
    if not st_flag: R,V = Args(7).tolist([(0,3),(3,3)]) 
    else: R,V = Args(8).tolist([(0,3),(3,3)])

    rdotv = R.dot(V)

    return rdotv


def p2DistanceEvent(ode, nd_dist=0.15, st_flag=False):
    """
    Integrator event for when P3 (spacecraft) reaches a supplied Non-Dimesnional distance from P2

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        general ASSET CR3BP ode
    
    nd_dist : float
        Non-dimensional distance of event

    st_flag : bool
        Boolean indicating whether or not the event is to be used with a Sundman Transformed ODE

    Returns
    -------
    dist : asset_asrl.VectorFunction
        Distance of P3 (spacecraft) from P2 minus nd_dist

    """
    if not st_flag: R,_ = Args(7).tolist([(0,3),(3,3)]) 
    else: R,_ = Args(8).tolist([(0,3),(3,3)])

    p2loc = np.array([1-ode.mu,0.0,0.0])
    dist  = (R-p2loc).norm() - nd_dist

    return dist


def p2CullEvent(ode, st_flag=False, alt=0.01, nd_dist=0.15):
    """
    Integrator event for when P3 (spacecraft) reaches a specified altitude or distance from P2

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        general ASSET CR3BP ode
    
    st_flag : bool
        Boolean indicating whether or not the event is to be used with a Sundman Transformed ODE
        
    alt : float
        Non-dimensional altitude of event

     nd_dist : float
        Non-dimensional distance of event

    Returns
    -------
    cull : asset_asrl.VectorFunction
        Distance and altitude of P3 (spacecraft) from P2
    """
        
    if not st_flag: X = Args(7)
    else: X = Args(8)

    altitude = (X.head3()-ode.P2).norm()-alt
    y = (X[1]-nd_dist)*(X[1]+nd_dist)
    
    cull = altitude*y

    return cull


def xAxisCrossEvent(st_flag=False):
    """
    Integrator event for when P3 (spacecraft) crosses the x-axis in the CR3BP frame

    Parameters
    ----------
    st_flag : bool 
        Boolean indicating whether or not the event is to be used with a Sundman Transformed ODE

    Returns
    -------
    y : VectorFunction.Arguments
        y value of state vector, equals 0 when spacecraft is on x-axis
    """
    if not st_flag: R,_ = Args(7).tolist([(0,3),(3,3)]) 
    else: R,_ = Args(8).tolist([(0,3),(3,3)])

    y = R[1]

    return y


def xCrossEvent(st_flag=False, val=0.0):
    """
    Integrator event for when P3 (spacecraft) state vector hits a specific x value in the CR3BP frame

    Parameters
    ----------
    st_flag : bool 
        Boolean indicating whether or not the event is to be used with a Sundman Transformed ODE
    val : float
        Event value for spacecraft state vector x component

    Returns
    -------
    x : VectorFunction.Arguments
        x value of state vector, equals 0 when spacecraft reaches given x distance
    """
    if not st_flag: R,_ = Args(7).tolist([(0,3),(3,3)]) 
    else: R,_ = Args(8).tolist([(0,3),(3,3)])

    x = R[0] - val

    return x


def yCrossEvent(st_flag=False, val=0.0):
    """
    Integrator event for when P3 (spacecraft) state vector hits a specific y value in the CR3BP frame

    Parameters
    ----------
    st_flag : bool 
        Boolean indicating whether or not the event is to be used with a Sundman Transformed ODE
    val : float
        Event value for spacecraft state vector y component

    Returns
    -------
    y : VectorFunction.Arguments
        y value of state vector, equals 0 when spacecraft reaches given y distance
    """ 
    if not st_flag: R,_ = Args(7).tolist([(0,3),(3,3)]) 
    else: R,_ = Args(8).tolist([(0,3),(3,3)])

    y = R[1] - val

    return y


def yAxisCrossEvent(st_flag=False):
    """
    Integrator event for when P3 (spacecraft) crosses the y-axis in the CR3BP frame

    Parameters
    ----------
    st_flag : bool 
        Boolean indicating whether or not the event is to be used with a Sundman Transformed ODE

    Returns
    -------
    x : VectorFunction.Arguments
        x value of state vector, equals 0 when spacecraft is on y-axis
    """
    if not st_flag: R,_ = Args(7).tolist([(0,3),(3,3)]) 
    else: R,_ = Args(8).tolist([(0,3),(3,3)])

    x = R[0]

    return x


def zCrossEvent(st_flag=False, val=0.0):
    """
    Integrator event for when P3 (spacecraft) state vector hits a specific z value in the CR3BP frame

    Parameters
    ----------
    st_flag : bool 
        Boolean indicating whether or not the event is to be used with a Sundman Transformed ODE
    val : float
        Event value for spacecraft state vector x component
        
    Returns
    -------
    z : VectorFunction.Arguments
        z value of state vector, equals 0 when spacecraft z - val = 0
    """
    if not st_flag: R,_ = Args(7).tolist([(0,3),(3,3)]) 
    else: R,_ = Args(8).tolist([(0,3),(3,3)])

    z = R[2] - val

    return z


def xCrossP2Event(ode, st_flag=False):
    """
    Integrator event for when P3 (spacecraft) crosses P2's X position

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        general ASSET CR3BP ode
    st_flag : bool 
        Boolean indicating whether or not the event is to be used with a Sundman Transformed ODE

    Returns
    -------
    x : VectorFunction.Arguments
        x value of state vector minus P2 x location
    """    
    
    if not st_flag: R,_ = Args(7).tolist([(0,3),(3,3)]) 
    else: R,_ = Args(8).tolist([(0,3),(3,3)])

    p2x = 1-ode.mu

    x = R[0] - p2x

    return x


def xCrossP1Event(ode, st_flag=False):
    """
    Integrator event for when P3 (spacecraft) crosses P1's X position

    Parameters
    ----------
    ode : asset_asrl.OptimalControl.ODEBase
        general ASSET CR3BP ode
    st_flag : bool 
        Boolean indicating whether or not the event is to be used with a Sundman Transformed ODE

    Returns
    -------
    x : VectorFunction.Arguments
        x value of state vector minus P1 x location
    """    
    
    if not st_flag: R,_ = Args(7).tolist([(0,3),(3,3)]) 
    else: R,_ = Args(8).tolist([(0,3),(3,3)])

    p1x = -ode.mu

    x = R[0] - p1x

    return x