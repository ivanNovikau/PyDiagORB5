import Mix as mix
import read_data as rd
import zf_gam as zf
import equil_profiles as equ
import ControlPlot as cpr
import transport
import work_profiles
import ITG_gamma as itg
import general
import common
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
    mix.reload_module(general)
    mix.reload_module(common)


# region --- NL SIMULATIONS ---

root_path = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/'
str_comp = {
    'path_ITPA_1': '/muT12-mun5/n80/',
    'project_name_1': 'ORIG',
    # 'path_ITPA_2': 'muT12-mun5/EGAM/n80-v55-s8',
    'path_ITPA_2': 'muT12-mun5/EGAM/f001-n80-v55-s8',
    'project_name_2': 'EGAM-f001',
}

# endregion

# region --- Project structure initialization ---

dd_init = {
    'a0': 0.2475,
    'R0': 0.88,
    'B0': 1.4348,
    'mass_pf': 2 * constants.proton_mass
}

# endregion

# region --- NON-LINEAR STRUCTURES ---

dd_orig = dict(dd_init)
dd_orig.update({
    'path': root_path + str_comp['path_ITPA_1'],
    'project_name': str_comp['project_name_1'],
})
dd_egam = dict(dd_init)
dd_egam.update({
    'path': root_path + str_comp['path_ITPA_2'],
    'project_name': str_comp['project_name_2'],
})
reload()
rd.init(dd_orig)
rd.init(dd_egam)

# endregion