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


def set_dd(dd_loc, root, path_loc, name_loc):
    proj_loc = {
        'path': root + path_loc,
        'name': name_loc,
    }
    dd_loc.update({'path': proj_loc['path'], 'project_name': proj_loc['name']})
    rd.init(dd_loc)
    return dd_loc


reload()
root_path = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/ITPA/'

dd_init = {
    'a0': 1.0,
    'R0': 10.0,
    'B0': 3.0,
    'mass_pf': 2 * constants.proton_mass
}

# --- Alessandro's simulation --
dd_noise = set_dd(dict(dd_init), root_path,
                  '/s256-chi384-phi192-t050-part8-krook-0-initnoise',
                  'ORB5')

