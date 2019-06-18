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
root_path = 'd:/Work-Projects/MyProgs/ORB_data/MPR/'
str_comp = {
    'path_ITPA_1': 'nled/kin/linear/REF-WORK-CPS2019/',
    'project_name_1': 'NLED',
    # 'path_ITPA_1': 'nled/adiab/vf8/',
    # 'project_name_1': 'NLED-adiab',
    # 'path_ITPA_2': 'nled/kin/linear/ref-f0-short/',
    # 'project_name_2': 'NLED-short',
    'path_ITPA_2': 'GAMs/adiab/q15/',
    'project_name_2': 'GAM,\ q = 1.5',
}
# endregion


# region --- Project structure initialization ---
dd_init_egam = {
    'a0': 0.482,
    'R0': 1.62,
    'B0': 2.2,
    'mass_pf': 2 * constants.proton_mass
}
dd_init_gam = {
    'a0': 0.5,
    'R0': 1.65,
    'B0': 2,
    'mass_pf': 2 * constants.proton_mass
}
# endregion


# region --- LINEAR STRUCTURES ---

reload()

# dd = dict(dd_init_egam)
# dd.update({
#     'path': root_path + str_comp['path_ITPA_1'],
#     'project_name': str_comp['project_name_1'],
# })
# rd.init(dd)
#
# dd_f0 = dict(dd_init_egam)
# dd_f0.update({
#     'path': root_path + str_comp['path_ITPA_2'],
#     'project_name': str_comp['project_name_2'],
# })
# rd.init(dd_f0)

dd_egam = dict(dd_init_egam)
dd_egam.update({
    'path': root_path + str_comp['path_ITPA_1'],
    'project_name': str_comp['project_name_1'],
})
rd.init(dd_egam)

dd_gam = dict(dd_init_gam)
dd_gam.update({
    'path': root_path + str_comp['path_ITPA_2'],
    'project_name': str_comp['project_name_2'],
})
rd.init(dd_gam)

# endregion