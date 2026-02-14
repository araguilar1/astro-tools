"""
Solve_NRHOs.py



Solve and save the 9:2 and 4:1 L2 Southern NRHOs

Initial Conditions from pages 245-246 of Dr. Emily Spreen's PhD Disseration:
TRAJECTORY DESIGN AND TARGETING FOR APPLICATIONS TO THE EXPLORATION PROGRAM IN CISLUNAR SPACE

"""

import numpy as np
from asset_asrl.Astro.AstroModels import CR3BP
import asset_asrl.Astro.Constants as c
from prettytable import PrettyTable

# Add utils to path
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parents[1].absolute()))

from CR3BP_Utils.solve_periodic_orbit import solvePeriodic
from data_utils import writeDataToPkl
from astro import jacobi
from general_utils import EARTH_MOON_CR3BP_DAY

########### Define CR3BP ODE ###########
ode = CR3BP(c.MuEarth, c.MuMoon, c.LD)

########### Create Initial Guesses ###########
# 9:2 L2 Southern NRHO
x0_92_L2s_NRHO    = np.zeros((7))
x0_92_L2s_NRHO[0] = 1.02134
x0_92_L2s_NRHO[2] = -0.18162
x0_92_L2s_NRHO[4] = -0.10176
x0_92_L2s_NRHO[5] = 9.76561e-07
tf_92_L2s_NRHO    = 1.50206

# 4:1 L2 Southern NRHO
x0_41_L2s_NRHO    = np.zeros((7))
x0_41_L2s_NRHO[0] = 1.03545
x0_41_L2s_NRHO[2] = -0.19003
x0_41_L2s_NRHO[4] = -0.13071
x0_41_L2s_NRHO[5] = 5.62991e-07
tf_41_L2s_NRHO    = 1.68981

########### Solve Each Orbit ###########
traj_92_L2s_NRHO, converged_92_L2s_NRHO = solvePeriodic(x0_92_L2s_NRHO, tf_92_L2s_NRHO, ode, fix_init=[1,5])
traj_41_L2s_NRHO, converged_41_L2s_NRHO = solvePeriodic(x0_41_L2s_NRHO, tf_41_L2s_NRHO, ode, fix_init=[1,5])

########### Validate Orbit Period and Jacobi Constant ###########
traj_92_jc = jacobi(traj_92_L2s_NRHO[0], ode.mu)
traj_41_jc = jacobi(traj_41_L2s_NRHO[0], ode.mu)

traj_92_period_days = traj_92_L2s_NRHO[-1][-1]*EARTH_MOON_CR3BP_DAY
traj_41_period_days = traj_41_L2s_NRHO[-1][-1]*EARTH_MOON_CR3BP_DAY

table = PrettyTable()
table.field_names = ["Orbit", "Computed JC", "Disseration JC", "Computed Period", "Disseration Period"]
table.add_row(['9:2 L2 Southern NRHO',traj_92_jc,3.04719,traj_92_period_days,6.5556])
table.add_row(['4:1 L2 Southern NRHO',traj_41_jc,3.03476,traj_41_period_days,7.375])
print(table)

########### Save Data ###########
traj_92_path = './data/Orbits/traj_92_L2s_NRHO_CR3BP.pkl'
writeDataToPkl(traj_92_L2s_NRHO, traj_92_path)

traj_41_path = './data/Orbits/traj_41_L2s_NRHO_CR3BP.pkl'
writeDataToPkl(traj_41_L2s_NRHO, traj_41_path)