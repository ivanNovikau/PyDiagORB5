import Mix as mix
import read_data as rd
import zf_gam as zf
import equil_profiles as equ
import ControlPlot as cpr
import transport
import matplotlib.pyplot as mpl
import numpy as np
from scipy import constants


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(zf)
    mix.reload_module(equ)
    mix.reload_module(cpr)
    mix.reload_module(transport)


# ----------------------------------
# --- TCV case ---
# path_Xavier = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/test-Krook/ref/'
path_Xavier = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/test-Krook/krook/'

dd = {
    'path': path_Xavier,
    'a0': 0.13,
    'R0': 1.3,
    'B0': 1.9,
    'mass_pf': 2 * constants.proton_mass
}

# ******************************************************
# *** Initialization of the simulation's data ***
reload()
rd.init(dd)

# ******************************************************