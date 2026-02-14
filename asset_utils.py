""""
asset_utils.py



Utilities built to easily interact with the ASSET Python Package

Functions
---------
runJetJobConvergedOnly : run a jet job and only return converged solutions
runJetJobAllSolutions  : run a jet job and return all solutions regardless of convergence status
"""

import asset_asrl as ast
vf      = ast.VectorFunctions
oc      = ast.OptimalControl
solvers = ast.Solvers
Args    = vf.Arguments
cflags  = solvers.ConvergenceFlags
tmodes  = oc.TranscriptionModes


def runJetJobConvergedOnly(ocps:list, ncores:int, printOutput:bool=True)->list:
    """
    Run a Jet job and return converged solutions

    Parameters
    ----------
    ocps : list
        list of ASSET oc.OptimalControlProblem
    ncores : int
        number of CPU cores to utilize
    printOutput : bool
        print Jet output

    Returns
    -------
    converged_ocps : list
        ocps that successfully converged
    """
    ocps = solvers.Jet.map(ocps, ncores, printOutput)

    converged_ocps = []
    for ocp in ocps:
        if ocp.optimizer.get_ConvergenceFlag() == cflags.CONVERGED:
            converged_ocps.append(ocp)

    return converged_ocps


def runJetJobAllSolutions(ocps:list, ncores:int, printOutput:bool=True)->list:
    """
    Run a Jet job and return all solutions regardless of convergence status
    
    Parameters
    ----------
    ocps : list
        list of ASSET oc.OptimalControlProblem
    ncores : int
        number of CPU cores to utilize
    printOutput : bool
        print Jet output

    Returns
    -------
    ocps : list
        list of ASSET oc.OptimalControlProblem
    """
    ocps = solvers.Jet.map(ocps, ncores, printOutput)

    return ocps