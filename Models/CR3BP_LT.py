""""
CR3BP_LT.py

Generic Low-Thrust ODEs

Classes
-------
CR3BP_LT    : Low-Thrust enabled CR3BP ODE
CR3BP_LT_ST : Low-Thrust enabled CR3BP ODE with Sundman Transform Implementation

"""
import asset_asrl.OptimalControl as oc
import asset_asrl.VectorFunctions as vf
import numpy as np

class CR3BP_LT(oc.ODEBase):

    def __init__(self, mu, tLt, tLTstar, lstar, tstar, isp, g0):
        """
        Low-Thrust augmented CR3BP Equations of Motion

        Parameters
        ----------
        mu : float
            System mass parameter
        tLt: float
            Dimensional low-thrust (N)
        tLTstar : float
            Low-thrust non-dimensionalized in length and time units
        lstar : float
            Characteristic length (m)
        tstar : float
            Characteristic time (s)
        isp : float
            engine specific impulse (s)
        g0 : float
            Acceleration due to gravity (9.81 m/s**2)

        Returns
        -------
        ODEBase Vector Function
        """

        self.mu = mu
        self.lstar = lstar
        self.tstar = tstar

        # ode is x,t,u
        x_vars = 7 # pos, vel, mass
        u_vars = 3 # control

        args = oc.ODEArguments(x_vars,u_vars)

        r = args.head3()
        v = args.segment3(3)
        u = args.tail3()

        x = args[0]
        y = args[1]
        xdot = args[3]
        ydot = args[4]
        m = args[6]

        rterms = vf.stack([2 * ydot + x,
                        -2.0 * xdot + y]).padded_lower(1)
        p1loc = np.array([-mu, 0, 0])
        p2loc = np.array([1.0 - mu, 0, 0])

        g1 = (r - p1loc).normalized_power3() * (mu - 1.0)
        g2 = (r - p2loc).normalized_power3() * (-mu)

        acc_lt = u*tLTstar/m # acceleration from low-thrust, non-dimensional 
        
        mdot   = - (u.norm() * tLt / (isp*g0)) * tstar
        
        acc = vf.sum([rterms,g1,g2,acc_lt])
        ode = vf.stack([v,acc,mdot])
        super().__init__(ode,x_vars,u_vars)


class CR3BP_LT_ST(oc.ODEBase): 
    
    def __init__(self, mu, tLt, tLTstar, lstar, tstar, isp, g0):
        """
        Low-Thrust augmented CR3BP Equations of Motion with Sundman Transform

        Parameters
        ----------
        mu : float
            System mass parameter
        tLt: float
            Dimensional low-thrust (N)
        tLTstar : float
            Low-thrust non-dimensionalized in length and time units
        lstar : float
            Characteristic length (m)
        tstar : float
            Characteristic time (s)
        isp : float
            engine specific impulse (s)
        g0 : float
            Acceleration due to gravity (9.81 m/s**2)

        Returns
        -------
        ODEBase Vector Function
        """
        
        self.mu = mu
        self.lstar = lstar
        self.tstar = tstar

        # ode is x,t,u
        x_vars = 8 # pos, vel, mass, st
        u_vars = 3 # control

        XtU = oc.ODEArguments(x_vars,u_vars)

        r = XtU.head3()
        v = XtU.segment3(3)
        u = XtU.tail3()

        x    = XtU[0]
        y    = XtU[1]
        xdot = XtU[3]
        ydot = XtU[4]
        m    = XtU[6]

        rterms = vf.stack([2 * ydot + x,
                        -2.0 * xdot + y]).padded_lower(1)
        p1loc = np.array([-mu, 0, 0])
        p2loc = np.array([1.0 - mu, 0, 0])

        g1 = (r - p1loc).normalized_power3() * (mu - 1.0)
        g2 = (r - p2loc).normalized_power3() * (-mu)

        acc_lt = u*tLTstar/m # acceleration from low-thrust, non-dimensional 
        
        mdot   = - (u.norm() * tLt / (isp*g0)) * tstar

        ST = vf.norm(r-np.array([1-mu,0.0,0.0])) # Sundman Transform term

        acc = vf.sum([rterms,g1,g2,acc_lt])
        ode = vf.stack([v*ST,acc*ST,mdot*ST,ST])
        super().__init__(ode,x_vars,u_vars)