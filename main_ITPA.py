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
    # 'path_ITPA_1': '/ES-lin-n26/',
    # 'project_name_1': 'ES-LINER-n26',
    # 'path_ITPA_2': '/EM/',
    # 'project_name_2': 'EM',
    # 'path_ITPA_2': '/EM-fast/',
    # 'project_name_2': 'EM',

    'path_ITPA_2': '/EM-woFast-canMaxw/',
    'project_name_2': 'EM',

    # 'path_ITPA_2': '/EM-m3670-woFast/',
    # 'project_name_2': 'EM-m3670',

    # 'path_ITPA_2': '/ES-flux-n25/',
    # 'project_name_2': 'ES-n25',

    # 'path_ITPA_2': '/ES-kT08/',
    # 'project_name_2': 'ES-kT08',

    # 'path_ITPA_2': '/ES-kn04/',
    # 'project_name_2': 'ES-kn04',
}

# endregion

# region --- LINEAR SIMULATIONS ---

# root_path = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/ITPA/'
# str_comp = {
#     # 'path_ITPA_1': '/linear/n30-m200/',
#     # 'path_ITPA_1': '/linear/n30-m3670-dt3/',
#
#     # 'path_ITPA_1': '/linear/ES/kt08-kn03/n30',
#     'path_ITPA_1': '/linear/ES/kn04-kt10/n10',
#
#     'project_name_1': 'LIN-n10',
# }

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

dd_em = dict(dd_init)
dd_em.update({
    'path': root_path + str_comp['path_ITPA_2'],
    'project_name': str_comp['project_name_2'],
})
rd.init(dd_em)

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

# region --- LINEAR STRUCTURES ---

# reload()
# dd = dict(dd_init)
# dd.update({
#     'path': root_path + str_comp['path_ITPA_1'],
#     'project_name': str_comp['project_name_1'],
# })
# rd.init(dd)

# endregion