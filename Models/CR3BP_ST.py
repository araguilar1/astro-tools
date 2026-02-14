"""
CR3BP_ST.py

Copy of the base ASSET Circular Restricted Three-Body Problem (CR3BP) with a Sundman Transform Implementation.

Classes
-------
    CR3BP_ST_FRAME : Generic Frame Implementation in ASSET format
    CR3BP_ST       : Generic ODE implementation of the CR3BP_ST_FRAME in ASSET format
"""
import numpy as np
import asset_asrl as ast
from asset_asrl.OptimalControl import ODEBase

vf = ast.VectorFunctions
oc = ast.OptimalControl
Args = vf.Arguments


class CR3BP_ST_FRAME:
    """"
    Basic Frame Implementation

    Parameters
    ----------
    mu1 : float
        P1 Gravitational Parameter
    mu2 : float
        P2 Gravitational Parameter
    lstar : float
        Characteristic Length of System
    """

    def __init__(self,mu1,mu2,lstar):

        self.P1mu  = mu1
        self.P2mu  = mu2
        self.lstar = lstar
        self.mstar = mu1+mu2
        self.mu    = mu2/self.mstar
        self.tstar = np.sqrt(lstar**3/self.mstar)
        self.vstar = lstar/self.tstar
        self.astar = self.mstar/lstar**2
        
        self.CalcPoints()

    def CalcPoints(self):
        """Calculates location of system Lagrange Points. Called
        when CR3BPFrame is initialized.


        Parameters
        ----------

        Returns
        -------

        """
        mu = self.mu
        self.P1 = np.array([-mu,0,0])
        self.P2 = np.array([1.0-mu,0,0])
        self.L4 = np.array([.5 - mu, np.sqrt(3.0) / 2.0, 0.0])
        self.L5 = np.array([.5 - mu, -np.sqrt(3.0) / 2.0, 0.0])
        
        ###L1
        gamma0 = pow((mu*(1.0 - mu)) / 3.0, 1.0 / 3.0)
        guess = gamma0 + 1.0
        while (abs(guess - gamma0) > 10e-15):
            gamma0 = guess
            guess = pow((mu*(gamma0 - 1)*(gamma0 - 1.0)) / (3.0 - 2.0 * mu - gamma0 * (3.0 - mu - gamma0)), 1.0 / 3.0);
        self.L1 = np.array([1 - mu - guess, 0, 0]);

		###L2
        gamma0 = pow((mu*(1.0 - mu)) / 3.0, 1.0 / 3.0);
        guess = gamma0 + 1.0;
        while (abs(guess - gamma0) > 10e-15):
            gamma0 = guess;
            guess = pow((mu*(gamma0 + 1)*(gamma0 + 1.0)) / (3.0 - 2.0 * mu + gamma0 * (3.0 - mu + gamma0)), 1.0 / 3.0);
        self.L2 =  np.array([1 - mu + guess, 0, 0]);

		#### L3
        gamma0 = pow((mu*(1.0 - mu)) / 3.0, 1.0 / 3.0);
        guess = gamma0 + 1.0;
        while (abs(guess - gamma0) > 10e-15): 
            gamma0 = guess;
            guess = pow(((1.0 - mu)*(gamma0 + 1)*(gamma0 + 1.0)) / (1.0 + 2.0 * mu + gamma0 * (2.0 + mu + gamma0)), 1.0 / 3.0);  
        self.L3 =  np.array([-mu - guess, 0, 0])

    def CR3BP_ST_EOMS(self,r,v):

        # Sundmann Transform term (distance from P2)
        ST = vf.norm(r-np.array([1-self.mu,0.0,0.0]))

        x    = r[0]
        y    = r[1]
        xdot = v[0]
        ydot = v[1]
        
        rterms = vf.stack([2*ydot+x,-2*xdot+y]).padded_lower(1)
        g1 = (r-self.P1).normalized_power3()*(self.mu-1.0)
        g2 = (r-self.P2).normalized_power3()*(-self.mu)

        accterms   = [g1,g2,rterms]
        acc = vf.sum(accterms)
        terms = [v*ST,acc*ST,ST] # add ST term here
        return vf.stack(terms)


class CR3BP_ST(ODEBase, CR3BP_ST_FRAME):
    def __init__(self,P1mu,P2mu,lstar):
        CR3BP_ST_FRAME.__init__(self,P1mu,P2mu,lstar)
        args = oc.ODEArguments(7,0)
        r,v = args.tolist([(0,3),(3,3)])
        ode = self.CR3BP_ST_EOMS(r,v)
        ODEBase.__init__(self,ode,7)