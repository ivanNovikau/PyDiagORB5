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
    'path_ITPA_2': '/EM/',
    'project_name_2': 'EM',
}
# endregion

# region --- Project structure initialization ---
dd_init = {
    'a0': 1.0,
    'R0': 10.0,
    'B0': 3.0,
    'mass_pf': 2 * constants.proton_mass
}

dd_es = dict(dd_init)
dd_es.update({
    'path': root_path + str_comp['path_ITPA_1'],
    'project_name': str_comp['project_name_1'],
})

dd_em = dict(dd_init)
dd_em.update({
    'path': root_path + str_comp['path_ITPA_2'],
    'project_name': str_comp['project_name_2'],
})
# endregion

# region --- Initialization of the simulation's data ---
reload()
rd.init(dd_es)
rd.init(dd_em)

# endregion