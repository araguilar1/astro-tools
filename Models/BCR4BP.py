"""
BCR4BP.py

Implementation of the Bicircular Restricted Four-Body Problem (BCR4BP) in ASSET Format

Classes
-------
    BCR4BPFrame : Generic Frame Implementation in ASSET format
    BCR4BP      : Generic ODE implementation of the BCR4BPFrame in ASSET format (Earth-Moon Barycenter)
    BCR4BPSB1   : Generic ODE implementation of the BCR4BPFrame in ASSET format (Sun-B1 Barycenter)
"""
import asset_asrl as ast
import numpy as np
from asset_asrl.OptimalControl import ODEBase

vf = ast.VectorFunctions
oc = ast.OptimalControl
Args = vf.Arguments

class BCR4BPFrame:
    """
    Basic Frame Implementation

    Parameters
    ----------
    mu1 : float
        Earth Gravitational Parameter
    mu2 : float
        Moon Gravitational Parameter
    mu4 : float
        P4 (Sun) Gravitational Parameter
    lstarEM : float
        Earth-Moon System Characteristic Length
    lstarSB1 : float:
        Sun-Earth-Moon Barycenter System Characteristic Length
    theta0 : float
        Initial Sun Angle (radians)
    eta : float
        Ratio of amount of perturbing force the Sun enforces on P3 (Spacecraft) in the Earth-Moon rotating frame (0=CR3BP,1=Full BCR4BP)
    """

    def __init__(self, mu1, mu2, mu4, lstarEM, lstarSB1, theta0, eta):

        self.P1mu = mu1 # Earth Graviational Parameter
        self.P2mu = mu2 # Moon Graviational Parameter
        self.P4mu = mu4 # Sun Graviational Parameter
        self.lstarEM = lstarEM # Mean distance from Earth to Moon
        self.lstarSB1 = lstarSB1 # Mean distance from Sun to Earth-Moon Barycenter

        self.mu = mu2 / (mu1 + mu2) # Earth-Moon mass ratio
        self.muprime = (mu1 + mu2) / (mu1 + mu2 + mu4) # Sun-Earth-Moon mass ratio

        self.tstarEM = np.sqrt((lstarEM ** 3) / (mu1 + mu2)) # Earth-Moon characteristic time
        self.tstarSB1 = np.sqrt((lstarSB1 ** 3) / (mu1 + mu2 + mu4)) # Sun-Earth-Moon Barycenter characteristic time

        self.m4 = mu4 / (mu1 + mu2) # Non-dimensional mass of Sun
        self.a4 = self.lstarSB1 / self.lstarEM # Non-dimensional distance from Sun to Earth-Moon Barycenter

        self.theta0 = theta0 # Initial Angle between Earth-Moon X-Axis and the Sun's Position (radians), also used as the initial angle between the Sun-B1 X-axis and the Earth-Moon X-Axis (radians)
        self.thetadot = np.sqrt((1.0+self.m4)/self.a4**3) - 1.0 # Mean Motion of Sun - 1.0
        self.eta = eta # Ratio of the amount of perturbiring force the Sun enforces on the s/c in the E-M frame

        self.CalcPoints()


    def CalcPoints(self):
        """
        
        Calculates location of system points in Earth-Moon frame. Called
        when BCR4BP is initialized.

        This function calculates the location of the points in the Earth-Moon
        rotating frame. Since the Sun will rotate about this system, only the initial position wrt the 
        Earth-Moon Barycenter is calculated.

        """

        mu = self.mu
        mu_prime = self.muprime
        theta0 = self.theta0

        # Earth-Moon Rotating Frame
        self.P1EM  = np.array([-mu, 0, 0]) # Constant
        self.P2EM  = np.array([1.0 - mu, 0, 0]) # Constant
        self.P4EM0 = self.a4*np.array([np.cos(theta0), np.sin(theta0), 0]) # Variable

        # Sun-B1 Rotating Frame
        self.P4SB1  = np.array([-mu_prime, 0, 0]) # Constant
        self.PB1SB1 = np.array([1.0 - mu_prime, 0, 0]) # Constant
        self.P1SB1_0  = np.array([1.0 - mu_prime - (1.0/self.a4)*mu*np.cos(theta0), -(1.0/self.a4)*mu*np.sin(theta0), 0]) # Variable    
        self.P2SB1_0  = np.array([1.0 - mu_prime + (1.0/self.a4)*(1.0-mu)*np.cos(theta0), (1.0/self.a4)*(1.0-mu)*np.sin(theta0), 0]) # Variable


    def HamiltonianEM(self, ndx, theta):
        """
        
        Hamiltonian of the Earth-Moon system

        """

        a4 = self.a4
        m4 = self.m4
        mu = self.mu

        x = ndx[0]
        y = ndx[1]
        z = ndx[2]

        # args = Args(6)

        # r = args.head3()
        # v = args.tail3()


        # P1EM = self.P1EM
        # P2EM = self.P2EM
        # P4EM = vf.stack

        # gt1 = (r-P1EM).inverse_norm()*(1.0-mu)
        # gt2 = (r-p2loc).inverse_norm()*(mu)
        # gt4 = (r-)

        gamma = (1.0-mu)/np.sqrt((x+mu)**2 + y**2 + z**2) + mu/np.sqrt((x-1.0+mu)**2 + y**2 + z**2) + (x**2 + y**2)/2.0 + m4/np.sqrt((x-a4*np.cos(theta))**2 + (y-a4*np.sin(theta))**2 + z**2) - m4*(x*np.cos(theta) + y*np.sin(theta))/a4**2
        v = np.linalg.norm(ndx[3:6])

        return 2*gamma - v**2


    def BCR4BPEOMS(self, r, v, theta, eta=1, otherAccs=[]):

        a4 = self.a4
        m4 = self.m4
        mu = self.mu
        P1EM = self.P1EM
        P2EM = self.P2EM

        x = r[0]
        y = r[1]
        xdot = v[0]
        ydot = v[1]
        

        P4EM = vf.stack([a4*vf.cos(theta), a4*vf.sin(theta)]).padded_lower(1)

        rterms = vf.stack([2.0*ydot+x, -2.0*xdot+y]).padded_lower(1)
        g1 = (r-P1EM).normalized_power3()*(mu-1.0)
        g2 = (r-P2EM).normalized_power3()*(-mu)
        g4 = (r-P4EM).normalized_power3()*(-m4)*eta
        g5 = vf.stack([-m4*vf.cos(theta)/(a4**2), -m4*vf.sin(theta)/(a4**2)]).padded_lower(1) * eta

        thetadot = np.sqrt((1.0+m4)/a4**3) - 1.0

        accterms = [g1, g2, g4, g5, rterms] + otherAccs
        acc = vf.sum(accterms)

        return vf.stack(v, acc, thetadot)


    def BCR4BPSB1EOMS(self, r, v, theta, otherAccs=[]):
        
        mu = self.mu
        muprime = self.muprime
        a4 = self.a4
        m4 = self.m4
        P4mu = self.P4mu
        P4SB1 = self.P4SB1
        
        x = r[0]
        y = r[1]
        xdot = v[0]
        ydot = v[1]

        P1SB1 = vf.stack([1.0-muprime-(1.0/a4*mu*vf.cos(theta)), -(1.0/a4*mu*vf.sin(theta))]).padded_lower(1)
        P2SB1 = vf.stack([1.0-muprime+(1.0/a4*(1.0-mu)*vf.cos(theta)), (1.0/a4)*(1.0-mu)*vf.sin(theta)]).padded_lower(1)
    
        rterms = vf.stack([2.0*ydot+x, -2.0*xdot+y]).padded_lower(1)
        g4 = (r-P4SB1).normalized_power3()*(muprime-1.0)
        g1 = (r-P1SB1).normalized_power3()*(muprime)*(mu-1.0)
        g2 = (r-P2SB1).normalized_power3()*(muprime*-mu)

        thetadot = 1.0 - np.sqrt((1.0+m4)/a4**3)
        # thetadot = np.sqrt((1.0+m4)/a4**3) - 1.0
        # thetadot = a4**3/np.sqrt(1.0+P4mu)

        accterms = [g1, g2, g4, rterms] + otherAccs
        acc = vf.sum(accterms)
        return vf.stack(v, acc, thetadot)


class BCR4BP(ODEBase, BCR4BPFrame):
    def __init__(self, P1mu, P2mu, P3mu, lstarEM, lstarSB1, theta0, eta):
        BCR4BPFrame.__init__(self, P1mu, P2mu, P3mu, lstarEM, lstarSB1, theta0, eta)

        args = oc.ODEArguments(7,0)

        r, v, theta = args.XVec().tolist([(0,3),(3,3),(6,1)])

        ode = self.BCR4BPEOMS(r, v, theta, eta=eta, otherAccs=[])

        ODEBase.__init__(self, ode, 7)


class BCR4BPSB1(ODEBase, BCR4BPFrame):
    def __init__(self, P1mu, P2mu, P3mu, lstarEM, lstarSB1, theta0, eta=1):
        BCR4BPFrame.__init__(self, P1mu, P2mu, P3mu, lstarEM, lstarSB1, theta0, eta)

        args = oc.ODEArguments(7,0)

        r, v, theta = args.XVec().tolist([(0,3),(3,3),(6,1)])

        ode = self.BCR4BPSB1EOMS(r, v, theta, otherAccs=[])

        ODEBase.__init__(self, ode, 7)