import Mix as mix
import read_data as rd
import zf_gam as zf
import equil_profiles as equ
import ControlPlot as cpr
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
    mix.reload_module(work_profiles)
    mix.reload_module(itg)


# ----------------------------------
# --- AUG20787 case ---
    
# --- LAPTOP ---
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/AUG20787/adiab/Krook/s81/'
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/AUG20787/adiab/Krook/s51-ns1024/'  # <-- !!!
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/AUG20787/adiab/Krook/s751-ns512/'
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/AUG20787/adiab/Krook/flux4-s51-ns1024/'

# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#         'turbulence/AUG20787/adiab/linear/n0-k2-cos/'
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#         'turbulence/AUG20787/adiab/linear/n0-k2-sin/'


# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#         'turbulence/AUG20787/adiab/linear/n10/'
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#         'turbulence/AUG20787/adiab/linear/n30/'
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
#         'turbulence/AUG20787/adiab/linear/n60/'
path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/' \
        'turbulence/AUG20787/adiab/linear/n90/'

# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'AUG20787/adiab/Krook/scan-rotation/pos/'
# path_AUG = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'AUG20787/adiab/Krook/scan-rotation/neg/'
    
# --- DRACO ---
#path_AUG = '/ptmp/ivannovi/turbulence/AUG20787/adiab/linear/n80'
#path_AUG = '/ptmp/ivannovi/turbulence/AUG20787/adiab/linear/n70'
#path_AUG = '/ptmp/ivannovi/turbulence/AUG20787/adiab/linear/n50'
#path_AUG = '/ptmp/ivannovi/turbulence/AUG20787/adiab/linear/n40'
#path_AUG = '/ptmp/ivannovi/turbulence/AUG20787/adiab/linear/n20'
#path_AUG = '/ptmp/ivannovi/turbulence/AUG20787/adiab/linear/n1'
path_AUG = '/ptmp/ivannovi/turbulence/AUG20787/adiab/linear/n90'

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
# Main principles:
# -> i can keep a local folder (on a given machine) with some global variables:
# paths to created projects, their names, list of variables, description of a project and
# its variables, project configuration
# -> Using main script without a project, one can save data to temporal files
# (to avoid overloading of the operational memory)
# -> With a project one can define its configuration, save data to permanent files to load
# these data afterwards

# Questions:
# How to set project configuration, that would be specific for a particular project without
# making complex the whole project?