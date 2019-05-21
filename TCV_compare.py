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


# region --- PROJECTS TO COMPARE ---
root_path = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/'
# str_comp = {
#     'path_TCV1': '/NL/muT12-mun5/n80/',
#     'project_name1': 'n=80',
#     'path_TCV2': '/NL/muT12-mun5/n80-modT/',
#     'project_name2': 'n=80,\ modified\ T',
# }
str_comp = {
    'path_TCV1': '/NL/muT12-mun5/n80/',
    'project_name1': 'n=80',
    'path_TCV2': '/NL/muT12-mun5/EGAM/n80-v55-s8/',
    'project_name2': 'n=80,\ EGAM',
}
# str_comp = {
#     'path_TCV1': '/NL/muT12-mun5/scan-dt-n80/dt15/',
#     'project_name1': 'dt = 15',
#     'path_TCV2': '/NL/muT12-mun5/n80/',
#     'project_name2': 'dt = 20',
#     'path_TCV3': '/NL/muT12-mun5/scan-dt-n80/dt25/',
#     'project_name3': 'dt = 25',
# }
# str_comp = {
#     'path_TCV1': '/NL/muT12-mun5/scan-nptot-n80/n1e8/',
#     'project_name1': 'N_i = 1e8',
#     'path_TCV2': '/NL/muT12-mun5/scan-nptot-n80/n4e8/',
#     'project_name2': 'N_i = 4e8',
#     'path_TCV3': '/NL/muT12-mun5/n80/',
#     'project_name3': 'N_i = 6e8',
# }

# endregion

# region --- Project structure initialization ---
dd_init = {
    'a0': 0.2475,
    'R0': 0.88,
    'B0': 1.4348,
    'mass_pf': 2 * constants.proton_mass
}

dd1 = dict(dd_init)
dd1.update({
    'path': root_path + str_comp['path_TCV1'],
    'project_name': str_comp['project_name1'],
})

dd2 = dict(dd_init)
dd2.update({
    'path': root_path + str_comp['path_TCV2'],
    'project_name': str_comp['project_name2'],
})

# dd3 = dict(dd_init)
# dd3.update({
#     'path': root_path + str_comp['path_TCV3'],
#     'project_name': str_comp['project_name3'],
# })
# endregion

# region --- Initialization of the simulation's data ---
reload()
rd.init(dd1)
rd.init(dd2)
# rd.init(dd3)
# endregion
