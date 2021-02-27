import Mix as mix
import read_data as rd
import zf_gam as zf
import equil_profiles as equ
import ControlPlot as cpr
import transport
import work_profiles
import ITG_gamma as itg
import common
import Global_variables as GLO
import Geom as GEO
import MPR as mpr
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
    mix.reload_module(common)
    mix.reload_module(GLO)
    mix.reload_module(GEO)
    mix.reload_module(mpr)

def set_dd(dd_loc, root, path_loc, name_loc):
    proj_loc = {
        'path': root + path_loc,
        'name': name_loc,
    }
    dd_loc.update({'path': proj_loc['path'], 'project_name': proj_loc['name']})
    rd.init(dd_loc)
    return dd_loc

reload()


# size of the machine, magnetic field at the axis and the mass of the main species:
dd_init = {
    'a0': 0.482,
    'R0': 1.62,
    'B0': 2.2,
    'mass_pf': 2 * constants.proton_mass
}

# indicate path to the folder, where your ORB5-data are placed
project_root_path       = './TEMPLATES/'
project_relative_path   = 'f001-mpr/'
project_name            = 'LINEAR\ ES\, ADIABATIC\ ELECTRONS'

# create a structure with project data
dd  = set_dd(dict(dd_init), project_root_path, project_relative_path, project_name) 