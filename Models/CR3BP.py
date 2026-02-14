"""
CR3BP.py



Copy of the base ASSET Circular Restricted Three-Body Problem (CR3BP), but allows users to specify mu, l*, and t*

Classes
-------
    CR3BP_FRAME : Generic Frame Implementation in ASSET format
    CR3BP       : Generic ODE implementation of the CR3BP_FRAME in ASSET format
"""

import asset_asrl as ast
from asset_asrl.OptimalControl import ODEBase
import numpy as np

vf = ast.VectorFunctions
Args = vf.Arguments
oc = ast.OptimalControl

class CR3BPFrame:
    """
    Basic Frame Implementation
    
    Attributes
    ----------
    

    """

    def __init__(self, mu, lstar, tstar):
        
        """
        
        Parameters
        ----------
        mu : float
            System Mass Parameter
        lstar : float
            Characteristic distance between primaries in kilometers (km)
        tstar : float
            Characteristic time of system in seconds (s)
        """
        self.mu = mu
        self.tstar = tstar
        self.lstar = lstar
        self.vstar = lstar/self.tstar

        
        self.CalcPoints()
    
    def CalcPoints(self):
        """
        Calculates location of system Lagrange Points. Called
        when CR3BPFrame is initialized.


        Parameters
        ----------
        None

        Returns
        -------
        None

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
        
    def GenLissajousImpl(self,func,X,xnd,znd,phideg,psideg,nplanrev,npo,t0 = 0):
        """
        Lissajou orbit implementation function. Usually not directly used,
        instead see GenL1Lissajous or GenL2Lissajous.


        Parameters
        ----------
        func: :obj:`function`
            Function used to evaluate Dynamics
        X : :obj:`list` of float
            State to evaluate Jacobian of fun
        xnd : float
            Size of x dimension of Lissajou orbit (non-dim)
        znd : float
            height of z dimension of Lissajou orbit (non-dim)
        phideg : float
            Phasing?
        psideg :float
            Phasing
        nplanrev : float
            Number of revolutions of Lissajou orbit
        npo : float
            ???
        t0 : float
            Initial time of Lissajou orbit (non-dim)
            
        Returns
        -------
        traj : :obj:`list` of :obj:`list`
            Lissajou orbit states

        """
        J = func.jacobian(X)
        Oxx = J[3,0]
        Oyy = J[4,1]
        Ozz = J[5,2]
        
        pi = np.pi;
        b1 = 2.0 - (Oxx + Oyy) / 2.0;
        b2sq = -Oxx * Oyy;
        s = np.sqrt(b1 + np.sqrt(b1*b1 + b2sq));
        b3 = (s*s + Oxx) / (2.0*s);
        pp = 2.0*pi / s;
        nu = np.sqrt(abs(Ozz));
        
        traj = []
        dtr = np.pi / 180.0;
        phi = phideg * dtr;
        psi = psideg * dtr;
        
        ynd = xnd*b3
        tt  = nplanrev * pp;
        dt  = tt / float(npo - 1);
        for i in range(0,npo): 
            st = np.zeros((7));
            ti  = float(i)*dt;
            st[0] = -(ynd / b3)*np.cos(s*ti + phi);
            st[1] = ynd * np.sin(s*ti + phi);
            st[2] = znd * np.sin(nu*ti + psi);
            st[3] = (ynd / b3)*s*np.sin(s*ti + phi);
            st[4] = ynd * s*np.cos(s*ti + phi);
            st[5] = znd * nu*np.cos(nu*ti + psi);
            st[6] =  t0 + ti;
            st[0:3] += X[0:3];
            traj.append(st)
        return traj;
    
    def GenL1Lissajous(self,xnd,znd,phideg,psideg,nplanrev,npo,t0 = 0):
        """Creates a  L1 Lissajou orbit


        Parameters
        ----------
        xnd : float
            Size of x dimension of Lissajou orbit (non-dim)
        znd : float
            height of z dimension of Lissajou orbit (non-dim)
        phideg : float
            Phasing?
        psideg :float
            Phasing
        nplanrev : float
            Number of revolutions of Lissajou orbit
        npo : float
            ???
        t0 : float
            Initial time of Lissajou orbit (non-dim)
            
        Returns
        -------
        traj : :obj:`list` of :obj:`list`
            Lissajou orbit states

        """
        args = Args(6)
        func = self.CR3BPEOMs(args.head(3),args.tail(3))
        X = np.zeros((6))
        X[0:3] = self.L1
        return self.GenLissajousImpl(func,X,xnd,znd,phideg,psideg,nplanrev,npo,t0)
    
    def GenL2Lissajous(self,xnd,znd,phideg,psideg,nplanrev,npo,t0 = 0):
        """Creates a L2 Lissajou orbit


        Parameters
        ----------
        xnd : float
            Size of x dimension of Lissajou orbit (non-dim)
        znd : float
            height of z dimension of Lissajou orbit (non-dim)
        phideg : float
            Phasing?
        psideg :float
            Phasing
        nplanrev : float
            Number of revolutions of Lissajou orbit
        npo : float
            ???
        t0 : float
            Initial time of Lissajou orbit (non-dim)
            
        Returns
        -------
        traj : :obj:`list` of :obj:`list`
            Lissajou orbit states

        """
        args = Args(6)
        func = self.CR3BPEOMs(args.head(3),args.tail(3))
        X = np.zeros((6))
        X[0:3] = self.L2
        return self.GenLissajousImpl(func,X,xnd,znd,phideg,psideg,nplanrev,npo,t0)
        
    def CalcSubPoint(self,func,IG):
        """Compute Sub Lagrange points, such as with a solar sail


        Parameters
        ----------
        func : :obj:`function`
            ASSET VectorFunction modeling dynamics
        IG : :obj:`list` of float
            Initial guess of location of sub Lagrange points
            
        Returns
        -------
        X : :obj:`list` of float 
            Position of sub Lagrange points

        """
        X = np.zeros((func.IRows()))
        X[0:3]=IG
        while(max(abs(func.compute(X)))> 1.0e-14):
            F = func.compute(X)[3:5]
            Jt = func.jacobian(X)
            J = np.array([[Jt[3,0],Jt[3,1]],
                          [Jt[4,0],Jt[4,1]]])
            X[0:2] = X[0:2] - np.dot(np.linalg.inv(J),F)
        return np.array(X[0:3])
    
    def CalcSubPoints(self,func):
        """Compute Sub Lagrange points, such as with a solar sail


        Parameters
        ----------
        func : :obj:`function`
            ASSET VectorFunction modeling dynamics
            
        Returns
        -------

        """
        self.SubL1 = self.CalcSubPoint(func,self.L1)
        self.SubL2 = self.CalcSubPoint(func,self.L2)
        self.SubL3 = self.CalcSubPoint(func,self.L3)
        self.SubL4 = self.CalcSubPoint(func,self.L4)
        self.SubL5 = self.CalcSubPoint(func,self.L5)
        
    def CR3BPEOMs(self,r,v,otherAccs = [],otherEOMs = []):
        """
        CR3BP Equations of motion


        Parameters
        ----------
        r : `list` of `float`
            Position components of state in the CR3BP
        v : `list` of `float`
            Position components of state in the CR3BP
        otherAccs : `list` of `float`, optional
            Other accelerations to be included in Equations of motion,
            such as low-thrust accelerations, perturbations etc.
        otherEOMs : `list` of `float`, optional
            Other equations of motion to be included in the dynamics
        
            
        Returns
        -------
        X : `asset_asrl.VectorFunction`
            ASSET VectorFunction output vector containing derivatives of states
            in the CR3BP

        """
        x    = r[0]
        y    = r[1]
        xdot = v[0]
        ydot = v[1]
        
       
        rterms = vf.stack([2*ydot+x,-2*xdot+y]).padded_lower(1)
        g1 = (r-self.P1).normalized_power3()*(self.mu-1.0)
        g2 = (r-self.P2).normalized_power3()*(-self.mu)

        accterms   = [g1,g2,rterms] + otherAccs
        acc = vf.sum(accterms)
        terms = [v,acc] + otherEOMs
        X = vf.stack(terms)
        return X


class CR3BP(ODEBase,CR3BPFrame):
    
    def __init__(self,mu,lstar,tstar):
        CR3BPFrame.__init__(self,mu,lstar,tstar)
        ###################################
        args = oc.ODEArguments(6,0)
        r = args.XVec().head3()
        v = args.XVec().tail3()
        ode = self.CR3BPEOMs(r,v,
                             otherAccs=[],
                             otherEOMs=[])
        ODEBase.__init__(self,ode,6)
        ###################################

