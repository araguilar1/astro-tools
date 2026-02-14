"""
solve_periodic_orbit.py



Solve a Periodic Orbit

Functions
---------
solvePeriodic : Solve Periodic Orbit (feasible solution)

"""
import numpy as np

def solvePeriodic(ig: np.array, tf: float, ode, print_level=0, fix_init=[0,1,2,3,6], fix_end=[1,3,5])->list:
    """
    Solve a Periodic Orbit in the CR3BP assuming a perpendicular x-axis crossing

    Parameters
    ----------
    ig : np.array
        Initial Guess
    tf : float
        Time of Flight initial guess
    ode : asset_asrl.OptimalControl.ODEBase
        Generic ASSET ODE
    print_level : int
        phase optimizer print level (0=full, >2=None)
    fix_init : list
        Index of variables to fix at start of orbit, e.g., 0 = fix X position

    Returns
    -------
    traj : list
        Solved Periodic Orbit
    """
    # define integrator
    integrator = ode.integrator('DOPRI87',1e-4)

    # integrate initial guess
    steps = 10000
    traj_guess = integrator.integrate_dense(ig, tf, steps)

    # create phase
    phase = ode.phase("LGL5")
    phase.Threads = 8
    
    nseg = 500

    phase.setTraj(traj_guess,nseg)
    for idx in fix_init:
        phase.addBoundaryValue("Front",[idx],[0.0])
    
    for idx in fix_end:
        phase.addBoundaryValue("Back",[idx],[0.0])

    # solve
    tol = 1e-12
    phase.optimizer.set_EContol(tol)
    phase.optimizer.PrintLevel = print_level
    phase.optimize()

    # return solution
    traj  = phase.returnTraj()
    cflag = phase.optimizer.get_ConvergenceFlag()

    return traj, cflag
