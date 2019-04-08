import Mix as mix
import read_data as rd
import zf_gam as zf
import equil_profiles as equ
import ControlPlot as cpr
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


# ----------------------------------
# --- AUG20787 case ---
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/AUG20787/adiab/Krook/s81/'
path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
           'turbulence/AUG20787/adiab/Krook/s51-ns1024/'
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/AUG20787/adiab/Krook/s751-ns512/'
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/AUG20787/adiab/Krook/flux4-s51-ns1024/'

dd = {
    'path': path_AUG,
    'a0': 0.5,
    'R0': 1.65,
    'B0': 2.0,
    'mass_pf': 2 * constants.proton_mass
}

# ******************************************************
# *** Initialization of the simulation's data ***
reload()
rd.init(dd)

# ******************************************************