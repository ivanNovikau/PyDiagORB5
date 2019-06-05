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
import MPR as mpr
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


# region --- LINEAR SIMULATIONS ---
root_path = 'd:/Work-Projects/MyProgs/ORB_data/MPR/nled/'
str_comp = {
    # 'path_ITPA_1': 'kin/linear/REF-WORK-CPS2019/',
    'path_ITPA_1': 'adiab/vf8/',
    'project_name_1': 'NLED-adiab',
    'path_ITPA_2': 'kin/linear/ref-f0-short/',
    'project_name_2': 'NLED-adiab',
}
# endregion


# region --- Project structure initialization ---
dd_init = {
    'a0': 0.482,
    'R0': 1.67,
    'B0': 2.2,
    'mass_pf': 2 * constants.proton_mass
}
# endregion


# region --- LINEAR STRUCTURES ---

reload()
dd = dict(dd_init)
dd.update({
    'path': root_path + str_comp['path_ITPA_1'],
    'project_name': str_comp['project_name_1'],
})
rd.init(dd)

dd_f0 = dict(dd_init)
dd_f0.update({
    'path': root_path + str_comp['path_ITPA_2'],
    'project_name': str_comp['project_name_2'],
})
rd.init(dd_f0)

# endregion
