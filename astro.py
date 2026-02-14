"""
Start astro.py: Script containing functions for astrodynamic applications

Author: Abram Aguilar

Class Objects and Method Objects have the following naming convention:
    Class Objects: UpperCamelCase
    Method Objects: lower_case_with_underscores
    ex. class CelestialObject, method axial_rotation_period()
    
Functions have the following naming convention:
    Functions: camelCase
    ex. function characteristicQuantities()

12/29/2019: Creation
"""

import numpy as np
from scipy.constants import G
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import newton


class CelestialObject:
    """
    Class for extracting celestial body data. 
    
    Attitributes:
    ------------
        body: str 
            Name of celestial body
        ID: int 
            ID number of celestial body
        data: list 
            List of data for celestial body

    Methods:
    --------
        axial_rotation_period(): 
            Returns axial rotation period of celestial body
        equatorial_radius(): 
            Returns equatorial radius of celestial body
        mu(): 
            Returns gravitational parameter of celestial body
        mass(): 
            Returns mass of celestial body
        orbit_semi_major_axis(): 
            Returns semi-major axis of celestial body's orbit
        orbital_period(): 
            Returns orbital period of celestial body
        orbit_eccentricity(): 
            Returns eccentricity of celestial body's orbit
        orbit_inclination(): 
            Returns inclination of celestial body's orbit
        name(): 
            Returns name of celestial body


    The data table below contains the following information:

        data[ID][0] = Axial Rotation Period (Rev/Day)
        data[ID][1] = Equatorial Radius (km)
        data[ID][2] = Gravitational Parameter (mu = G*m (km^3/s^2) )
        data[ID][3] = Semi-major Axis of Orbit (km)
        data[ID][4] = Orbital Period (days)
        data[ID][5] = Orbital Eccentricity
        data[ID][6] = Inclination of Orbit to Ecliptic (deg)

        Where ID = Body ID # specified by an integer
    """
    raw_data = [
                [0.039401, 695700.0, 132712440041.94, '--', '--', '--', '--', "Sun"],  # Sun
                [0.0366, 1738.0, 4902.8, 384400.0, 27.32, 0.0549, 5.145, "Moon"],  # Moon 
                [0.017051, 2440.0, 22031.78, 4902.8, 87.97, 0.205647, 7.0039, "Mercury"],  # Mercury
                [0.004115, 6051.89, 324858.59, 108208441.28, 224.7, 0.006794, 3.39449, "Venus"],  # Venus
                [1.002737, 6378.14, 398600.44, 149657721.28, 365.26, 0.016192, 0.0045, "Earth"],  # Earth
                [0.974700, 3394.0, 42828.38, 227937168.37, 686.98, 0.093343, 1.84814, "Mars"],  # Mars
                [2.418111, 71492.0, 126686534.91, 778350260.22, 4332.59, 0.048708, 1.30394, "Jupiter"],  # Jupiter
                [2.252205, 60268.0, 37931207.80, 1433203470.67, 10755.70, 0.050663, 2.48560, "Saturn"],  # Saturn
                [-1.3921089, 25559.0, 5793951.32, 2863429878.70, 30685.40, 0.048551, 0.77151, "Uranus"],  # Uranus
                [1.489754, 24766.0, 6835099.50, 4501859020.15, 60189.0, 0.007183, 1.77740, "Neptune"],  # Neptune
                [0.156563, 1188.30, 869.34, 6018076570.89, 91101.50, 0.258313, 17.22524, "Pluto"],  # Pluto
                [0.15625, 605.0, 102.27, 19596.84, 6.39, 0.00005, 112.89596, "Charon"],  # Charon (about Pluto)
                [0.5291, 43.33, 0.0003, 48690.0, 24.85, 0.238214, 112.88839, "Nix"],  # Nix (about Pluto)
                [2.328, 65.0, 0.000320, 64738.0, 38.2, 0.0058652, 0.24200, "Hydra"],  # Hydra (about Pluto)
                ['Synch', 2634.0, 9891.0, 1070042.8, 7.15, 0.0006, 0.186, "Ganymede"],  # Ganymede (about Jupiter)
                ['Synch', 2575.5, 8978.13, 1221870.0, 15.95, 0.0288, 0.28, "Titan"],  # Titan (about Jupiter)
                ['Synch', 788.9, 235.4, 435800.0, 8.71, 0.0022, 0.1, "Titania"],  # Titania (about Uranus)
                [2.644860, 469.7, 62.63, 413968739.37, 1680.22, 0.076103, 10.6007, "Ceres"],  # Ceres 
                ['Synch', 252.30, 1.21135, 238040.0, 1.370218, 0.0047, 0.009, "Enceladus"],  # Enceladus (about Saturn)
                ['Synch', 13.10, 0.000721, 9377.20, 0.32, 0.0151, 1.082, "Phobos"],  # Phobos (about Mars)
                ['Synch', 1352.60, 1432.93, 354760.0, 5.876854, 0.000016, 156.834, "Triton"],  # Triton (about Neptune)
                [0.05992, 2403.0, 7181.32, 1883000.0, 16.69, 0.007, 0.281, "Callisto"],  # Callisto (about Jupiter)
                [4.2625, 151.959, 3203.56, 671100.0, 3.552, 0.0094, 7.483, "Europa"], # Europa (about Jupiter)
                ]

    def __init__(self, name):
        self.body = name
        if self.body.upper() == "SUN":
            self.ID = 0
        elif self.body.upper() == "MOON":
            self.ID = 1
        elif self.body.upper() == "MERCURY":
            self.ID = 2
        elif self.body.upper() == "VENUS":
            self.ID = 3
        elif self.body.upper() == "EARTH":
            self.ID = 4
        elif self.body.upper() == "MARS":
            self.ID = 5
        elif self.body.upper() == "JUPITER":
            self.ID = 6
        elif self.body.upper() == "SATURN":
            self.ID = 7
        elif self.body.upper() == "URANUS":
            self.ID = 8
        elif self.body.upper() == "NEPTUNE":
            self.ID = 9
        elif self.body.upper() == "PLUTO":
            self.ID = 10
        elif self.body.upper() == "CHARON":
            self.ID = 11
        elif self.body.upper() == "NIX":
            self.ID = 12
        elif self.body.upper() == "HYDRA":
            self.ID = 13
        elif self.body.upper() == "GANYMEDE":
            self.ID = 14
        elif self.body.upper() == "TITAN":
            self.ID = 15
        elif self.body.upper() == "TITANIA":
            self.ID = 16
        elif self.body.upper() == "CERES":
            self.ID = 17
        elif self.body.upper() == "ENCELADUS":
            self.ID = 18
        elif self.body.upper() == "PHOBOS":
            self.ID = 19
        elif self.body.upper() == "TRITON":
            self.ID = 20
        elif self.body.upper() == "CALLISTO":
            self.ID = 21
        elif self.body.upper() == "EUROPA":
            self.ID = 22
        else:
            print("Error: Celesital Object not found")
            self.ID = 999999

        self.data = self.raw_data[self.ID]

    def axial_rotation_period(self):
        return self.data[0]

    def equatorial_radius(self):
        return self.data[1]

    def mu(self):
        return self.data[2]

    def mass(self):
        return self.data[2] / (G / 1000**3)

    def orbit_semi_major_axis(self):
        return self.data[3]

    def orbital_period(self,unit="day"):
        if unit == "day":
            return self.data[4]
        elif unit == "year":
            return self.data[4]/365.25
        else:
            print("Warning: Unit not found. Defaulting to days.")
            return self.data[4]

    def orbit_eccentricity(self):
        return self.data[5]

    def orbit_inclination(self):
        return self.data[6]

    def name(self):
        return self.data[7]


def characteristicQuantities(P1, P2):
    """
    Returns array of characteristic quantities of the CR3BP for a primary-secondary system:
    1. System Mass Ratio (mu) 
    2. Characteristic Length
    3. Characteristic Time
    """
    m_star = P1.mass() + P2.mass()
    m2 = P2.mass()
    mu = m2 / m_star
    lstar = P2.orbit_semi_major_axis()
    tstar = np.sqrt(lstar ** 3 / (G/1000**3) / m_star)
    return [mu, lstar, tstar]


def cr3bpOde(ndx, t, mu):
    """
    CR3BP Equations of Motion

    :param ndx: State Vector
    :param t: Time
    :param mu: Mass Ratio

    :return: State Vector Derivative

    """
    # Position
    x = ndx[0]
    y = ndx[1]
    z = ndx[2]

    # Velocity
    vx = ndx[3]
    vy = ndx[4]
    vz = ndx[5]

    # Scalar distances from P1 to P3 and P2 to P3 respectively
    d = np.sqrt((x + mu) ** 2 + y ** 2 + z ** 2)
    r = np.sqrt((x + mu - 1) ** 2 + y ** 2 + z ** 2)

    # Acceleration
    ax = -(1 - mu) * (x + mu) / (d ** 3) - mu * (x - 1 + mu) / (r ** 3) + 2 * vy + x
    ay = -(1 - mu) * y / (d ** 3) - mu * y / (r ** 3) - 2 * vx + y
    az = -(1 - mu) * z / (d ** 3) - mu * z / (r ** 3)

    return np.array([vx, vy, vz, ax, ay, az])


def cr3bpOdeStm(ndx, t, mu):
    
    # Position
    x = ndx[0]
    y = ndx[1]
    z = ndx[2]

    # Velocity
    vx = ndx[3]
    vy = ndx[4]
    vz = ndx[5]

    # STM
    phi = np.reshape(ndx[6:], (6, 6))

    # Scalar Distances from P1 to P3, and P2 to P3 respectively
    d = np.sqrt((x + mu) ** 2 + y ** 2 + z ** 2)
    r = np.sqrt((x + mu - 1) ** 2 + y ** 2 + z ** 2)

    d3, d5 = d**3, d**5
    r3, r5 = r**3, r**5

    sigXX = 1 - (1-mu)/d3 - mu/r3 + 3*(1-mu)*(x+mu)**2/d5 + 3*mu*(x-1+mu)**2/r5
    sigYY = 1 - (1-mu)/d3 - mu/r3 + 3*(1-mu)*y**2/d5 + 3*mu*y**2/r5
    sigZZ = -(1-mu)/d3 - mu/r3 + 3*(1-mu)*z**2/d5 + 3*mu*z**2/r5
    sigXY = 3*(1-mu)*(x+mu)*y/d5 + 3*mu*(x-1+mu)*y/r5
    sigXZ = 3*(1-mu)*(x+mu)*z/d5 + 3*mu*(x-1+mu)*z/r5
    sigYZ = 3*(1-mu)*y*z/d5 + 3*mu*y*z/r5

    # Jacboian
    A = np.array([ [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                 [sigXX, sigXY, sigXZ, 0.0, 2.0, 0.0],
                 [sigXY, sigYY, sigYZ, -2.0, 0.0, 0.0],
                 [sigXZ, sigYZ, sigZZ, 0.0, 0.0, 0.0] ])

    # Acceleration 
    ax = -(1 - mu) * (x + mu) / d3 - mu * (x - 1 + mu) / r3 + 2 * vy + x
    ay = -(1 - mu) * y / d3 - mu * y / r3 - 2 * vx + y
    az = -(1 - mu) * z / d3 - mu * z / r3

    yd = np.array([vx, vy, vz, ax, ay, az])
    phiDotVec = np.dot(A, phi).reshape(36) 
    yd = np.hstack((yd, phiDotVec))
    return yd


def uStar(ndx, mu):
    """ Returns the pseudo-potential of the CR3BP """
    x = ndx[0]
    y = ndx[1]
    z = ndx[2]
    d = np.sqrt((x + mu) ** 2 + y ** 2 + z ** 2)
    r = np.sqrt((x + mu - 1) ** 2 + y ** 2 + z ** 2)
    return (1 - mu) / d + mu / r + (x ** 2 + y ** 2) / 2


def jacobi(ndx, mu):
    """Returns the Jacobi constant for some coordinates and mass ratio"""
    ustar = uStar(ndx, mu)
    C = 2 * ustar - np.dot(ndx[3:], ndx[3:])
    # C = 2 - uStar - np.linalg.norm(ndx[3:]) ** 2 (Also works, just to show different methods)
    return C


def lagrangePoints(mu):
    """Returns the x,y coordinates of the 5 Lagrange points for any system"""
    tol = 1e-12

    def colin(x):
        d = abs(x + mu)
        r = abs(x - 1 + mu)
        a = x
        b = (1 - mu) * (x + mu) / (d ** 3)
        c = mu * (x - 1 + mu) / (r ** 3)
        return a - b - c

    def fd_colin(x):
        xm = x - tol
        xp = x + tol
        ym = colin(xm)
        y = colin(x)
        yp = colin(xp)
        dm = (y - ym) / (x - xm)
        dp = (yp - y) / (xp - x)
        return (dm + dp) / 2

    # L1
    def f1(gam):
        x = 1 - mu - gam
        return colin(x)

    def f1p(gam):
        x = 1 - mu - gam
        return -fd_colin(x)

    gam1 = mu
    delta = 100
    while (abs(delta) > tol):
        delta = f1(gam1) / f1p(gam1)
        gam1 -= delta
    x1 = 1 - mu - gam1

    # L2
    def f2(gam):
        x = 1 - mu + gam
        return colin(x)

    def f2p(gam):
        x = 1 - mu + gam
        return fd_colin(x)

    gam2 = mu
    delta = 100
    while (abs(delta) > tol):
        delta = f2(gam2) / f2p(gam2)
        gam2 -= delta
    x2 = 1 - mu + gam2

    # L3
    def f3(gam):
        x = -mu - gam
        return colin(x)

    def f3p(gam):
        x = -mu - gam
        return -fd_colin(x)

    gam3 = 1
    delta = 100
    while (abs(delta) > tol):
        delta = f3(gam3) / f3p(gam3)
        gam3 -= delta
    x3 = -mu - gam3

    # L4/L5
    x45 = 0.5 - mu
    y4 = np.sqrt(3) / 2
    y5 = -np.sqrt(3) / 2

    return np.array([[x1, 0], [x2, 0], [x3, 0], [x45, y4], [x45, y5]])


def lagrangePoints2(mu):

    def collinear(x, mu):
        d = np.abs(x + mu)
        r = np.abs(x - 1 + mu)
        a = x
        b = (1 - mu) / d ** 3 * (x + mu)
        c = mu / r ** 3 * (x - 1 + mu)
        return a - b - c
    
    x1 = newton(func=collinear, x0=0, args=(mu,))
    x2 = newton(func=collinear, x0=1, args=(mu,))
    x3 = newton(func=collinear, x0=-1, args=(mu,))

    x45 = 0.5 - mu
    y4 = np.sqrt(3) / 2
    y5 = -np.sqrt(3) / 2

    return np.array([[x1, 0], [x2, 0], [x3, 0], [x45, y4], [x45, y5]])


def nonDim(x, lstar, tstar, mu, shift=True):
    """Takes a dimensional state (km, km/s) and non-dimensionalizes it using the CR3BP
        Characteristic Quantities
        
        if shift = True, origin shifted to system barycenter"""

    ndx      = np.copy(x)
    ndx[0:3] = ndx[0:3] / lstar
    ndx[3:6] = ndx[3:6] * tstar / lstar
    if shift: ndx[0] -= mu
    return ndx


def dim(ndx, lstar, tstar, mu, shift=True):
    """Converts a non-dimensional CR3BP state to a dimensional state (km, km/s)
    
    if shift = True, origin shifted to primary"""

    x = np.copy(ndx)
    if shift: x[0] += mu
    x[0:3] = x[0:3] * lstar
    x[3:6] = x[3:6] * lstar / tstar
    return x


def cr3bpToJ2k(rp2, vp2, direction=0):
    """" 
    Returns the Rotation Matrix for the following transformation:   
    CR3BP State -> J2000 State
    States must be dimensionalized first.
    Rotation matrix is epoch dependent, need spice and info from DExxx file.

    if direction = 0, then the transformation is from CR3BP to J2000, i.e. returns R
    if direction = 1, then the transformation is from J2000 to CR3BP, i.e. returns R^-1
    Based on algorithm in section 2.3.2 of Thomas Pavlak's PhD disseration. 
    'TRAJECTORY DESIGN AND ORBIT MAINTENANCE STRATEGIES IN MULTI-BODY DYNAMICAL REGIMES'
    """
    lstar = np.linalg.norm(rp2)
    h_ = np.cross(rp2, vp2)
    h = h_/lstar**2
    xhat = rp2/lstar 
    zhat = h/np.linalg.norm(h)
    yhat = np.cross(zhat, xhat)
    R0 = np.array([xhat[0], yhat[0], zhat[0], xhat[1], yhat[1], zhat[1], xhat[2], yhat[2], zhat[2]]).reshape(3,3)
    R1 = np.zeros((3,3))
    hn = np.linalg.norm(h)
    R2 = np.array([hn*yhat[0], -hn*xhat[0], 0, hn*yhat[1], -hn*xhat[1], 0, hn*yhat[2], -hn*xhat[2], 0]).reshape(3,3)
    Rtop = np.hstack((R0, R1))
    Rbot = np.hstack((R2, R0))
    R = np.vstack((Rtop, Rbot))
    
    if direction == 0:
        return R
    elif direction == 1:
      return np.linalg.inv(R)


def periodicOrbitDf():
    """
    Function returns a dataframe of period Lagrange Orbit information created by Dan Grebow.
        No input is required.
        NOTE: Assumes that the csv file is in the same directory as this script"""
    data = r'periodicLagrangeOrbits.csv'
    df = pd.read_csv(data)
    return df


def plotZvc(mu, C, fill=0):
    """ Plots the ZVC surface for a given mu and Jacobi constant. NOTE: uses matplotlib by default"""
    ax = plt.gca()

    npoints = 1000
    x = np.linspace(-1.5, 1.5, npoints)
    y = np.linspace(0.0, 1.5, npoints)
    ZZ = np.empty((npoints, npoints))
    ZZ[:] = np.NaN
    zzPos = np.zeros((npoints, npoints))
    eps = 0.1
    for i in range(npoints):
        for j in range(npoints):
            ndx = np.array([x[i], y[j], 0.0])
            temp = jacobi(ndx, mu)
            if (C >= temp) and (temp > (C * (1 - eps))):
                ZZ[j, i] = temp
                zzPos[j, i] = 1
    
    if fill == 1:
        ax.contourf(x, y, ZZ)
        ax.contourf(x, -y, ZZ)
    else:
        ax.contour(x, y, ZZ)
        ax.contour(x, -y, ZZ)

def twoBodyOde(dx, t, mu):
    r = dx[0:3]
    r3 = np.linalg.norm(r)**3
    v = dx[3:6]
    a = (-mu / r3) * r
    return np.hstack((v, a))


def lagrangePointEigendecomp(ndx, mu):
    x = ndx[0]
    y = ndx[1]
    z = ndx[2]

    d = np.sqrt((x + mu) ** 2 + y ** 2 + z ** 2)
    r = np.sqrt((x + mu - 1) ** 2 + y ** 2 + z ** 2)

    d3, d5 = d ** 3, d ** 5
    r3, r5 = r ** 3, r ** 5

    Uxx = 1 - (1 - mu) / d3 - mu / r3 + 3 * (1 - mu) * (x + mu) ** 2 / d5 + 3 * mu * (x - 1 + mu) ** 2 / r5
    Uyy = 1 - (1 - mu) / d3 - mu / r3 + 3 * (1 - mu) * y ** 2 / d5 + 3 * mu * y ** 2 / r5
    Uzz = -(1 - mu) / d3 - mu / r3 + 3 * (1 - mu) * z ** 2 / d5 + 3 * mu * z ** 2 / r5
    Uxy = 3 * (1 - mu) * (x + mu) * y / d5 + 3 * mu * (x - 1 + mu) * y / r5
    Uxz = 3 * (1 - mu) * (x + mu) * z / d5 + 3 * mu * (x - 1 + mu) * z / r5
    Uyz = 3 * (1 - mu) * y * z / d5 + 3 * mu * y * z / r5

    A = np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                  [Uxx, Uxy, Uxz, 0.0, 2.0, 0.0],
                  [Uxy, Uyy, Uyz, -2.0, 0.0, 0.0],
                  [Uxz, Uyz, Uzz, 0.0, 0.0, 0.0]])
    eigVal, eigVec = np.linalg.eig(A)

    return eigVal, eigVec


def lyapunovConstraints(ndx):
    x0 = np.copy(ndx)
    # y = xDot = 0.0
    return np.array([x0[1], x0[3]])


def xCrossing(t, ndx, mu):
    return ndx[1]
xCrossing.terminal = True


def solveLyapunov(ndx, t, mu=0.0121, tol=1e-12, direction=1):
    xCrossing.direction = direction
    xF = solve_ivp(cr3bpOdeStm, [0, t], ndx, args=[mu], method='DOP853', atol=tol, rtol=tol, events=xCrossing, dense_output=True)
    t0 = np.copy(xF.t[-1])
    x0 = xF.y[:,0]
    F = lyapunovConstraints(xF.y[:, -1])
    Fnorm = np.linalg.norm(F)
    while Fnorm>tol:
        x = xF.y[0, -1]
        y = xF.y[1, -1]
        yDot = xF.y[4, -1]
        d = np.sqrt((x+mu)**2 + y**2)
        r = np.sqrt((x-1+mu)**2 + y**2)
        xDotDot = -(1 - mu) * (x + mu) / (d ** 3) - mu * (x - 1 + mu) / (r ** 3) + 2 * yDot + x
        phi = xF.y[6:,-1]
        phi = np.reshape(phi, (6,6))


        F = lyapunovConstraints(xF.y[:,-1])
        Fnorm = np.linalg.norm(F)
        DF = np.array([[phi[1,4], yDot],
                             [phi[3,4], xDotDot]])
        dx = np.linalg.inv(DF).dot(F)

        x0[4] -= dx[0]
        t0    -= dx[1]
        xF = solve_ivp(cr3bpOdeStm, [0, t0], x0, args=[mu], method='DOP853', atol=tol, rtol=tol, dense_output=True)
    # Return full trajectory and period at the end
    xF = solve_ivp(cr3bpOdeStm, [0, 2*t0], x0, args=[mu], method='DOP853', atol=tol, rtol=tol, dense_output=True)
    return xF.y, xF.t[-1]


def lyapunovContinuation(ndx, t, mu=0.0121, dx=0.001, lim=0.9, direction=1):
    trajList = []
    pList = []
    xCrossing.direction = direction
    x0 = np.copy(ndx)
    t0 = np.copy(t)
    # first orbit
    xF1, IP = solveLyapunov(x0, t0, mu=mu, direction=direction)
    trajList.append(xF1)
    pList.append(IP)
    sign = np.sign(lim - trajList[-1][0][0])
    signLast = sign
    while sign == signLast:
       # initial guess is result of last orbit
        IG = np.copy(trajList[-1][0:,0])
        tIG = np.copy(IP/2)

        IG[0] += dx
        try:
            traj, IP = solveLyapunov(IG, tIG, mu=mu, direction=direction)

            trajList.append(np.copy(traj))
            signLast = sign
            sign = np.sign(lim - trajList[-1][0][0])
        except:
            continue
    return trajList


def haloConstraints(ndx):
    # y = xDot = zDot = 0.0
    return np.array([ndx[1], ndx[3], ndx[5]])


def solveHalo(ndx, t,  mu=0.0121, direction=-1):
    xCrossing.direction = direction
    # Initial guess
    xF = solve_ivp(cr3bpOdeStm, [0, t], ndx, args=[mu], method='DOP853', rtol=1e-12, atol=1e-12, events=xCrossing, dense_output=True)
    t0 = xF.t[-1]
    x0 = xF.y[:,0]
    F = haloConstraints(xF.y[:, -1])
    Fnorm = np.linalg.norm(F)
    while Fnorm>1e-12:
        x = xF.y[0, -1]
        y = xF.y[1, -1]
        z = xF.y[2, -1]
        yDot = xF.y[4, -1]
        d = np.sqrt((x+mu)**2 + y**2 + z**2)
        r = np.sqrt((x-1+mu)**2 + y**2 + z**2)
        xDotDot = -(1 - mu) * (x + mu) / (d ** 3) - mu * (x - 1 + mu) / (r ** 3) + 2 * yDot + x
        zDotDot = -(1 - mu) * z / (d ** 3) - mu * z / (r ** 3)
        phi = xF.y[6:,-1]
        phi = np.reshape(phi, (6,6))

        F = haloConstraints(xF.y[:, -1])
        Fnorm = np.linalg.norm(F)
        DF = np.array([
            [phi[1,0], phi[1,4], yDot],
            [phi[3,0], phi[3,4], xDotDot],
            [phi[5,0], phi[5,4], zDotDot]])

        dx = np.linalg.inv(DF) @ F

        x0[0] -= dx[0]
        x0[4] -= dx[1]
        t0    -= dx[2]

        xF = solve_ivp(cr3bpOdeStm, [0, t0], x0, args=[mu], method='DOP853', rtol=1e-12, atol=1e-12, dense_output=True)
    # Return full traj
    xF = solve_ivp(cr3bpOdeStm, [0, 2*t0], x0, args=[mu], method='DOP853', rtol=1e-12, atol=1e-12, dense_output=True)
    return xF.y, xF.t[-1]


def haloContinuation(ndx, t, mu=0.0121, dx=-0.001, lim=0.04, direction=-1):
    trajList = []
    traj1, IP = solveHalo(ndx, t, direction=direction)
    trajList.append(traj1)
    sign = np.sign(lim - trajList[-1][2][0])
    signLast = sign
    while sign == signLast:
        IG = np.copy(trajList[-1][0:,0])
        tIG = IP

        IG[2] += dx # adjust along Z dir
        try:
            traj, IP = solveHalo(IG, tIG, mu=mu, direction=direction)
            trajList.append(traj)

            signLast = sign
            sign = np.sign(lim - trajList[-1][2][0])
        except:
            continue
    return trajList

""" 
End astro.py
"""
