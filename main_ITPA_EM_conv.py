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

root_path = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/ITPA/'
str_comp = {
    'path_ITPA_1': '/ES/',
    'project_name_1': 'ES',

    'path_ITPA_2': '/EM-woFast-canMaxw',
    'project_name_2': 'EM-m200',

    'path_ITPA_3': '/EM-m800-ni1-ne4/',  # exploded
    'project_name_3': 'EM-m800-ni1-ne4',

    'path_ITPA_4': '/EM-m800-ni3-ne10/',  # exploded
    'project_name_4': 'EM-m800-ni3-ne10',

    'path_ITPA_5': '/EM-m800',  # exploded
    'project_name_5': 'EM-m800',

    'path_ITPA_6': '/EM-m800-dt05',  # exploded
    'project_name_6': 'EM-m800-dt05',

    'path_ITPA_7': '/EM-m800-ns2-dt05',
    'project_name_7': 'EM-m800-ns2-dt05',
}

# endregion

# region --- Project structure initialization ---

dd_init = {
    'a0': 1.0,
    'R0': 10.0,
    'B0': 3.0,
    'mass_pf': 2 * constants.proton_mass
}

# endregion

# region --- NON-LINEAR STRUCTURES ---

reload()

dd_es = dict(dd_init)
dd_es.update({
    'path': root_path + str_comp['path_ITPA_1'],
    'project_name': str_comp['project_name_1'],
})
rd.init(dd_es)

dd_m200 = dict(dd_init)
dd_m200.update({
    'path': root_path + str_comp['path_ITPA_2'],
    'project_name': str_comp['project_name_2'],
})
rd.init(dd_m200)

# dd_m800_ni1_ne4 = dict(dd_init)
# dd_m800_ni1_ne4.update({
#     'path': root_path + str_comp['path_ITPA_3'],
#     'project_name': str_comp['project_name_3'],
# })
# rd.init(dd_m800_ni1_ne4)

# dd_m800_ni3_ne10 = dict(dd_init)
# dd_m800_ni3_ne10.update({
#     'path': root_path + str_comp['path_ITPA_4'],
#     'project_name': str_comp['project_name_4'],
# })
# rd.init(dd_m800_ni3_ne10)
#
# dd_m800 = dict(dd_init)
# dd_m800.update({
#     'path': root_path + str_comp['path_ITPA_5'],
#     'project_name': str_comp['project_name_5'],
# })
# rd.init(dd_m800)
#
# dd_m800_dt05 = dict(dd_init)
# dd_m800_dt05.update({
#     'path': root_path + str_comp['path_ITPA_6'],
#     'project_name': str_comp['project_name_6'],
# })
# rd.init(dd_m800_dt05)

dd_m800_ns2 = dict(dd_init)
dd_m800_ns2.update({
    'path': root_path + str_comp['path_ITPA_7'],
    'project_name': str_comp['project_name_7'],
})
rd.init(dd_m800_ns2)

# dd_es_flux = dict(dd_init)
# dd_es_flux.update({
#     'path': root_path + str_comp['path_ITPA_2'],
#     'project_name': str_comp['project_name_2'],
# })
# rd.init(dd_es_flux)

# dd_kt08 = dict(dd_init)
# dd_kt08.update({
#     'path': root_path + str_comp['path_ITPA_2'],
#     'project_name': str_comp['project_name_2'],
# })
# rd.init(dd_kt08)

# dd_kn04 = dict(dd_init)
# dd_kn04.update({
#     'path': root_path + str_comp['path_ITPA_2'],
#     'project_name': str_comp['project_name_2'],
# })
# rd.init(dd_kn04)

# endregion
