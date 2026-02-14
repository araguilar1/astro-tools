"""
LVLH_CR3BP.py



Implementation of the Local Vertical Local Horizontal (LVLH) Model in the Circular Restricted Three-Body Problem (CR3BP) 
in ASSET format. This frame is leveraged to study relative motion of a chaser/loiterer spacecraft about a chief/target 
spacecraft in the CR3BP.

Classes:
-------
    LVLHFrame : Generic Frame Implementation in ASSET format
    LVLHCR3BP : Generic ODE implementation of the LVLHFrame in ASSET format

"""
import asset_asrl as ast
import asset_asrl.Astro.Constants as c
from asset_asrl.OptimalControl import ODEBase
import numpy as np

vf = ast.VectorFunctions
oc = ast.OptimalControl

class LVLHFrame:
    '''
    Basic Frame Implementation

    Parameters
    ----------

    p1mu : float
        P1 Gravitational Parameter
    p2mu : float
        P2 Gravitational Parameter
    lstar : float
        Characteristic Length of System

    '''

    def __init__(self, p1mu, p2mu, lstar):

        self.p1mu = p1mu
        self.p2mu = p2mu
        self.lstar = lstar
        self.mu = p2mu / (p1mu + p2mu)
    
    def dd(self, vec, vecdot):
        mag = vf.norm(vec)
        return (vecdot/(mag*mag*mag)+vec*(-3/(mag*mag*mag*mag*mag)*(vec.dot(vecdot))))
    
    def dd_np(self, vec, vecdot):
        mag = np.linalg.norm(vec)
        return (vecdot/(mag*mag*mag)+vec*(-3/(mag*mag*mag*mag*mag)*(vec.dot(vecdot))))
    
    def RPO_CR3BP_EOMS(self, r, v, pr, pv):
        """
        LVLH Frame Equations of Motion for the Loitering Spacecraft
        
        Parameters
        ----------
            r : np.array
                Loiterer 1x3 position state
            v : np.array
                Loiterer 1x3 velocity state
            pr : np.array
                Target 1x3 position state
            pv : np.array
                Target 1x3 velocity state

        Returns
        -------
            acc : asset_asrl.VectorFunction
                ASSET Vector Function ODE describing equation of motion
        """

        mu = self.mu

        h = r.cross(v)
        r_mag = vf.norm(r)
        h_mag = vf.norm(h) 
        # v_mag = vf.norm(v)

        # assume CR3BP for now, angular velocities
        w_MI = np.array([0, 0, 1])
        wd_MI = np.array([0, 0, 0])
        wdd_MI = np.array([0, 0, 0])

        # target motion relative to moon
        r_em = np.array([-1, 0, 0]) # check this
        r_em_mag = np.linalg.norm(r_em)

        r_acc = (vf.cross(-2*w_MI,v) - vf.cross(wd_MI,r) - vf.cross(w_MI,vf.cross(w_MI,r)) - mu*(r.normalized_power3())
                -(1-mu)*((r+r_em).normalized_power3()-r_em/(r_em_mag**3)))


        # define LVLH frame in Moon frame coordinates
        k_hat = -(r.normalized()) # MOON
        j_hat = -(h.normalized()) # MOON 
        i_hat = j_hat.cross(k_hat) # MOON

        hd = r.cross(r_acc) 
        hd_mag = (-hd).dot(j_hat)

        v_em = np.array([0,0,0]) 
        v_mag = (-v).dot(k_hat)

        r_j1 = vf.cross(-2*w_MI,r_acc)+vf.cross(-3*wd_MI,v)-vf.cross(wdd_MI,r)-vf.cross(wd_MI,vf.cross(w_MI,r))-vf.cross(w_MI,vf.cross(wd_MI,r))-vf.cross(w_MI,vf.cross(w_MI,v))
        r_j2 = -mu*(self.dd(r, v))-(1-mu)*(self.dd(r+r_em, v+v_em))+(1-mu)*self.dd_np(r_em,v_em) 

        r_jerk = (r_j1+r_j2)
        w_LM = vf.stack([-h_mag/(r_mag*r_mag), ((-r_mag/(h_mag*h_mag)*h).dot(r_acc))]).padded_upper(1) # local

        wd_LM = vf.stack([-1/r_mag*(hd_mag/r_mag+2*v_mag*w_LM[1]), (v_mag/r_mag-2*hd_mag/h_mag)*w_LM[2]-((r_mag/(h_mag*h_mag)*h).dot(r_jerk))]).padded_upper(1) # local
        wd_MI = np.array([0, 0, 0])


        w_MI_L = vf.stack([vf.dot(w_MI,i_hat), vf.dot(w_MI,j_hat), vf.dot(w_MI,k_hat)])
        w_LI = w_LM+w_MI_L
        wd_LI = wd_LM+wd_MI-vf.cross(w_LM,w_MI_L) # local

        # # need, r and r_em in LVLH frame
        r_L = vf.stack([-r_mag]).padded_upper(2) # since in -k direction
        r_em_L = vf.stack([vf.dot(r_em,i_hat), vf.dot(r_em,j_hat), vf.dot(r_em,k_hat)])

        p_acc = (vf.cross(-2*w_LI,pv)-vf.cross(wd_LI,pr)-vf.cross(w_LI,vf.cross(w_LI,pr))+mu*(r_L.normalized_power3()-(r_L+pr).normalized_power3())
                +(1-mu)*((r_L+r_em_L).normalized_power3()-(r_L+pr+r_em_L).normalized_power3()))

        acc = vf.stack([v,r_acc,pv,p_acc])

        return acc


class LVLHCR3BP(ODEBase, LVLHFrame):
    """
    Basic ASSET ODE Implementation

    Parameters
    ----------
        ODEBase (asset_asrl.OptimalControl.ODEBase): Generic ASSET ODEBase Class
        LVLHFrame (class LVLHFrame): LVLHFrame Class including equations of motion

    Returns
    -------
        ODEBase.__init__ (asset_asrl.OptimalControl.ODEBase): Initialized Generic ASSET ODEBase Class
    """

    def __init__(self, p1mu, p2mu, lstar):
        LVLHFrame.__init__(self, p1mu, p2mu, lstar)

        args = oc.ODEArguments(12,0)

        r, v, pr, pv = args.XVec().tolist([(0,3),(3,3),(6,3),(9,3)])

        ode = self.RPO_CR3BP_EOMS(r, v, pr, pv)

        ODEBase.__init__(self, ode, 12)