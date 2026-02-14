""""
general_utils.py



Contains general utilities, constants, and unit conversions

Functions
---------
normalize                    : Normalize a vector
findClosestPoint6d           : Use KDTree to find closest points in configuration space (6D)
findClosestPoint3d           : Use KDTree to find closest points in configuration space (3D)
findClosestConnection6d      : Use KDTree to find closest connection in configuration space (6D)
findClosestPointBySort       : Use List sort to find closest point from a set of target states (3D)
findClosestConnectionBySort  : Use List sort to find closest connection in configuration space (3D)
computeOsculatingElements    : Compute Osculating Orbital Elements from CR3BP State Vector
queryJpl3BodyApi             : Query the JPL 3Body API to get 3-Body Initial Conditions

Classes
------- 
Jpl3BodyData                 : Holds the JPL 3Body API returned data

"""
import requests
import json
import numpy as np
from scipy.spatial import KDTree
import asset_asrl.Astro.Constants as c

from asset_asrl.Astro import cartesian_to_classic_true

M_TO_KM = 1e-3
KM_TO_M = 1e3
RAD_TO_DEG = 180.0/np.pi
DEG_TO_RAD = np.pi/180.0
ARCSEC_TO_RAD = 4.84814e-6
RAD_TO_ARCSEC = 1.0/ARCSEC_TO_RAD

EARTH_MOON_CR3BP_MU    = c.MuMoon/(c.MuEarth+c.MuMoon)
EARTH_MOON_CR3BP_MSTAR = c.MuEarth+c.MuMoon
EARTH_MOON_CR3BP_LSTAR = c.LD
EARTH_MOON_CR3BP_TSTAR = np.sqrt(EARTH_MOON_CR3BP_LSTAR ** 3 / EARTH_MOON_CR3BP_MSTAR)
EARTH_MOON_CR3BP_DAY   = EARTH_MOON_CR3BP_TSTAR / c.day
EARTH_MOON_CR3BP_HOUR  = EARTH_MOON_CR3BP_DAY / 24.0

SUN_EARTH_CR3BP_MU    = c.MuEarth/(c.MuEarth+c.MuSun)
SUN_EARTH_CR3BP_MSTAR = c.MuSun + c.MuEarth
SUN_EARTH_CR3BP_LSTAR = c.AU
SUN_EARTH_CR3BP_TSTAR = np.sqrt(SUN_EARTH_CR3BP_LSTAR ** 3 / SUN_EARTH_CR3BP_MSTAR)
SUN_EARTH_CR3BP_DAY   = SUN_EARTH_CR3BP_TSTAR / c.day
SUN_EARTH_CR3BP_HOUR  = SUN_EARTH_CR3BP_DAY / 24.0 


def normalize(x):
    """
    Normalize a vector

    Parameters
    ----------
    x : np.array or list
        vector to normalize

    Returns
    -------
    unit : np.array
        unit vector

    """
    unit = np.array(x)/np.linalg.norm(x)
    return unit


def findClosestPoint6d(x, target):
    """
    Use a KDTree to find the closest point on target from x
    considering both position and velocity

    Parameters
    ----------
    x : np.array or list
        6x1 Position/Velocity State Vector
    target : np.array or list
        Array of 6x1 Pos/Vel State Vectors

    Returns
    -------
    closest_point : np.array or list
        6x1 State Vector that is closest in configuration space to x
    
    """
    target_states = [x[:6] for x in target]

    target_array  = np.array(target_states)[:,:6]

    tree = KDTree(target_array)

    _, idx = tree.query(x)

    closest_point = target[idx]

    return closest_point


def findClosestPoint3d(x, target):
    """
    Use a KDTree to find the closest point on target from x
    only considering position

    Parameters
    ----------
    x : np.array or list
        3x1 Position State Vector
    target : np.array or list
        Array of 3x1 Pos State Vectors
    
    Returns
    -------
    closest_point : np.array or list
        3x1 State Vector that is closest in configuration space to x
    
    """
    target_states = [x[:3] for x in target]

    target_array  = np.array(target_states)[:,:3]

    tree = KDTree(target_array)

    _, idx = tree.query(x)

    closest_point = target[idx]

    return closest_point


def findClosestConnection6d(crossings_1, crossings_2, returnList=False):
    """
    Use a KDTree to Determine the closest connection on a plane given two sets of 
    plane-crossings. Considers both position and velocity when determining distance.

    Parameters
    ----------
    crossings_1 : list 
        List of Plane Crossings
    crossings_2 : list 
        List of Plane Crossings
    returnList : bool 
        Return two closests point (False), or list of closests points (True) 

    Returns
    -------
    crossings : list 
        Two closests points or list of closests points

    """
    array1 = np.array(crossings_1)[:,:6]
    array2 = np.array(crossings_2)[:,:6]

    tree = KDTree(array2)
    dist, idx = tree.query(array1)
    
    if not returnList:
        
        min_index = np.argmin(dist)

        crossings = [crossings_1[min_index], crossings_2[idx[min_index]]]

    else:
        results = [(dist[i], crossings_1[i], crossings_2[idx[i]]) for i in range(len(dist))]
        results.sort(key=lambda x: x[0])
        crossings = [(result[1], result[2]) for result in results]
    
    return crossings


def findClosestConnection3d(crossings_1, crossings_2, returnList=False):
    """
    Use a KDTree to Determine the closest connection on a plane given two sets of 
    plane-crossings. Considers only position when determining distance.

    Parameters
    ----------
    crossings_1 : list 
        List of Plane Crossings
    crossings_2 : list 
        List of Plane Crossings
    returnList : bool 
        Return two closests point (False), or list of closests points (True) 

    Returns
    -------
    crossings : list 
        Two closests points or list of closests points

    """
    array1 = np.array(crossings_1)[:,:3]
    array2 = np.array(crossings_2)[:,:3]

    tree = KDTree(array2)
    dist, idx = tree.query(array1)
    
    if not returnList:
        
        min_index = np.argmin(dist)

        crossings = [crossings_1[min_index], crossings_2[idx[min_index]]]

    else:
        results = [(dist[i], crossings_1[i], crossings_2[idx[i]]) for i in range(len(dist))]
        results.sort(key=lambda x: x[0])
        crossings = [(result[1], result[2]) for result in results]
    
    return crossings


def findClosestPointBySort(x, target):
    """
    Find the closest point on target from x only considering position using list sort

    Parameters
    ----------
    x : np.array or list
        State Vector
    target : np.array or list
        Array of State Vectors
    
    Returns
    -------
    closest_point : np.array or list
        State Vector from target that is closest in configuration space to x
    
    """

    disti =  []
    for i, s in enumerate(target):
        dist = np.linalg.norm(x[0:3]-s[0:3])
        disti.append([dist,i])

    disti.sort(key=lambda x: x[0])
    closest_point = target[disti[0][1]]

    return closest_point


def findClosestConnectionBySort(x_1, x_2):
    """
    Find closest connection from the endpoints of two sets of orbits or manifolds

    Parameters
    ----------
    x_1 : list or np.array
        First data set of orbits or manifolds
    x_2 : list or np.array
        Second data set of orbits or manifolds

    Returns 
    -------
    
    closest_connection : list
        Two closests points of each data set
    """
    distij = []
    for i in range(len(x_1)):
        for j in range(len(x_2)):
            dist = np.linalg.norm(x_1[i][-1][0:6] - x_2[i][-1][0:6])
            distij.append([dist,i,j])

    distij.sort(key=lambda x: x[0])

    closest_connection = [x_1[distij[0][1]], x_2[distij[0][2]]]

    return closest_connection


def computeOsculatingElements(ndx, ode, mu, center='P2'):
    """
    Compute instantenous Osculating Orbital Elements with ASSET
    built-in cartesian_to_classic_true() function. Assumes 
    ODE class has P1, P2, and Libration Point Locations.

    Parameters
    ----------
    ndx : np.ndarray 
        Non-dimensional CR3BP state vector
    ode : asset_asrl.OptimalControl.ODEBase
        CR3BP ODE
    mu : float
        Gravitational Parameter
    center : str
        Point at which to center the state vector on, e.g., P2=Moon

    Raises
    ------
    ValueError : If center location not found in ODE class

    Returns
    -------
    elems : np.ndarray
        Osculating Orbit Elements (angles in radians) in following order
            - a (sma)
            - e (eccentricity)
            - i (inclination)
            - Omega (RAAN)
            - w (argument of pergiee)
            - v (true anomaly)
    """

    # Non-dim CR3BP State
    x = np.copy(ndx)

    # center about specific point
    match center:
        case 'P2':
            x[0:3] -= ode.P2
        case 'P1': 
            x[0:3] -= ode.P1
        case 'L1':
            x[0:3] -= ode.L1
        case 'L2':
            x[0:3] -= ode.L2
        case 'L3':
            x[0:3] -= ode.L3
        case 'L4':
            x[0:3] -= ode.L4
        case 'L5':
            x[0:3] -= ode.L5
        case _:
            raise ValueError("Invalid Center Location")
    
    # Define CR3BP Z-Axis
    z_axis = np.array([0.0, 0.0, 1.0])
    
    # Compute Angle (gamma) between State Vector and Z-axis
    rvec = x[0:3]
    rvec_norm = np.linalg.norm(rvec)
    rdotz = np.dot(rvec,z_axis)
    gamma = np.arccos(rdotz/rvec_norm)

    # Formulate rotation matrix to rotate CR3BP state to arbitrary inertial frame
    rotMat = np.array([
                    [np.cos(gamma), -np.sin(gamma), 0.0, 0.0, 0.0, 0.0],
                    [np.sin(gamma), np.cos(gamma), 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                    [-np.sin(gamma), -np.cos(gamma), 0.0, np.cos(gamma), -np.sin(gamma), 0.0],
                    [np.cos(gamma), -np.sin(gamma), 0.0, np.sin(gamma), np.cos(gamma), 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
    
    # Compute the inertial State
    inert_state = rotMat @ x[0:6]

    # Dimensionalize the inertial state
    inert_state[0:3] *= ode.lstar
    inert_state[3:6] *= ode.vstar

    elems = cartesian_to_classic_true(inert_state, mu)

    return elems


def generateEvenSphericalPoints(n=50):
    """
    Generate evenly spaced points on a sphere using the Fibonacci Lattice

    Parameters
    ----------
    n : int
        number of points

    Returns
    -------
    x, y, z : tuple
        x, y, z coordinates of evenly distributed points

    """

    golden_ratio = (1+5**0.5)/2
    i = np.arange(0, n)

    theta = 2 * np.pi * i / golden_ratio
    phi = np.arccos(1 - 2*(i+0.5)/n)

    x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)

    return x, y, z


def queryJpl3BodyApi(sys="earth-moon", family="halo", libr="2", branch="N", periodmin="", periodmax="", periodunits="", jacobimin="", jacobimax="", stabmin="", stabmax=""):
    """
    Get information from the JPL Horizons periodic 3-body orbits API.

    Parameters
    ----------
    sys : str 
        three-body system defined in lower-case as “primary-secondary,” e.g. earth-moon, mars-phobos, sun-earth.
    family : str 
        name of the orbit family: halo,vertical,axial,lyapunov,longp,short,butterfly,dragonfly,resonant,dro,dpo,lpo
    libr : str 
        libration point. Required for lyapunov,halo (1,2,3), longp, short (4,5), and axial,vertical (1,2,3,4,5).
    branch : str 
        branch of orbits within the family: N/S for halo,dragonfly,butterfly, E/W for LPO, and pq integer sequence for resonant (e.g., 12 for 1:2).
    periodmin : str 
        minimum period (inclusive). Units defined by periodunits.
    periodmax : str 
        maximum period (inclusive). Units defined by periodunits.
    periodunits : str 
        units of pmin and pmax: s for seconds, h for hours, d for days, TU for nondimensional.
    jacobimin : str 
        minimum Jacobi constant (inclusive). Nondimensional units.
    jacobimax : str 
        maximum Jacobi constant (inclusive). Nondimensional units.
    stabmin : str 
        minimum stability index (inclusive).
    stabmax : str 
        maximum stability index (inclusive).

    Returns
    -------
    data : json 
        JSON object containing the requested data.
    
    """
    # Inputs
    args = locals()

    # Define the base URL for the JPL Horizons API.
    baseUrl = "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?"
    for key in args:
        if args[key]:
            baseUrl += "{}={}&".format(key, args[key])

    baseUrl = baseUrl[:-1]

    # Get the data from the JPL Horizons API.
    r = requests.get(baseUrl)
    data = json.loads(r.text)
    # Check response status
    if r.status_code != 200:
        print("Error: {}".format(r.status_code))
        return None
    else:
        return data


class Jpl3BodyData:
    """
    Class to store the information from the JPL Horizons periodic 3-body orbits API return.

    Parameters
    ----------
    inData : json
        JSON object from JPL API query

    Methods
    -------
    sys()       : returns system string name, e.g. 'earth-moon'
    primary()   : returns primary string name, e.g. 'Earth'
    secondary() : returns secondary string name, e.g. 'Moon'
    family()    : returns family string name, e.g. 'halo'
    libpoint()  : retruns libration point string name, e.g. 'L2'
    mu()        : returns mass parameter float, e.g. 0.0121 
    lpoints()   : returns list of libration point locations
    ics()       : returns list of orbit family initial conditions
    """

    def __init__(self, inData):
        self.data = inData

    def sys(self):
        return self.data['system']['name']

    def primary(self):
        return self.data['system']['name'].split('-')[0]

    def secondary(self):
        return self.data['system']['name'].split('-')[1]

    def family(self):
        return self.data['family']

    def libpoint(self):
        return self.data['libration_point']

    def mu(self):
        return float(self.data['system']['mass_ratio'])

    def lpoints(self):
        lpoints = [self.data['system']['L{}'.format(i)] for i in range(1,6)]
        return [list(map(float, pos)) for pos in lpoints]

    def ics(self):
        initial_conditions = self.data['data']
        return [list(map(float, ic)) for ic in initial_conditions]
    

if __name__ == "__main__":

    import plotly.graph_objects as go

    n_ranges = [50,100,250,500,1000]

    for n in n_ranges:
        x,y,z = generateEvenSphericalPoints(n=n)

        fig = go.Figure()

        fig.add_trace(go.Scatter3d(x=x,y=y,z=z,mode='markers',marker=dict(size=5),showlegend=False))

        fig.show()