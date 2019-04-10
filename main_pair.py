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
# --- AUG20787 case ---
# path_pair = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'pair-plasma/NL/debye-001/krook-1e4/'
# path_pair = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'pair-plasma/NL/debye-001/krook-3e4/'
# path_pair = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'pair-plasma/NL/debye-001/krook-3e5/'
path_pair = 'd:/Work-Projects/MyProgs/ORB_data/' \
           'pair-plasma/NL/debye-300/krook-3e4/'

dd = {
    'path': path_pair,
    'a0': 0.612,
    'R0': 1.7,
    'B0': 2.0,
    'mass_pf': constants.electron_mass
}

# ******************************************************
# *** Initialization of the simulation's data ***
reload()
rd.init(dd)

# ******************************************************
