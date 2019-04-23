import Mix as mix
import read_data as rd
import zf_gam as zf
import equil_profiles as equ
import ControlPlot as cpr
import transport
import work_profiles
import ITG_gamma as itg
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
    mix.reload_module(work_profiles)
    mix.reload_module(itg)


## --- NL SIMULATIONS ---
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/ref-results-43516/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/n128/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/n80/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/n40/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/n60/'

## --- LINEAR GAM ---
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-GAM/n0/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-GAM/n0-Krook/'

## --- LINEAR ITG ---
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n58'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n60'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n62'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n80'


## --- LINEAR EGAM ---
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/flat-v4'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/flat-v8'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-v4-f01'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-v8-f01'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#             'TCV/linear/scan-EGAM/ge-v65-f02'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/gc-v4-f02'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-v2-f02'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-v4-f02'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-vf4-T2'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-vf4-T4'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-v6-f02'  # <-- !!!
path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
            'TCV/linear/scan-EGAM/gc-v6-f02'
path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
            'TCV/linear/scan-EGAM/ge-v55-f02'



dd = {
    'path': path_TCV,
    'a0': 0.2475,
    'R0': 0.88,
    'B0': 1.4348,
    'mass_pf': 2 * constants.proton_mass
}

# ******************************************************
# *** Initialization of the simulation's data ***
reload()
rd.init(dd)

# ******************************************************

