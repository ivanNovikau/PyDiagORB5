import Mix as mix
import read_data as rd
import zf_gam as zf
import numpy as np
from scipy import constants


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(zf)

dd = {}
path_TCV = '/ptmp/ivannovi/turbulence/TCV/ref-results-43516/'

# ----------------------------------
# --- TCV case ---
dd['path'] = path_TCV
dd['a0'] = 0.2475
dd['R0'] = 0.88
dd['B0'] = 1.4348
dd['mass_pf'] = 2 * constants.proton_mass

# ----------------------------------


# ******************************************************
# *** Initialization of the simultion's data ***
rd.init(dd)
# ******************************************************

