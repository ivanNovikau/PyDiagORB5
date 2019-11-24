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


def set_dd(dd_loc, root, path_loc, name_loc):
    proj_loc = {
        'path': root + path_loc,
        'name': name_loc,
    }
    dd_loc.update({'path': proj_loc['path'], 'project_name': proj_loc['name']})
    rd.init(dd_loc)
    return dd_loc


reload()
root_ref = 'd:/Work-Projects/MyProgs/ORB_data/MPR/'
root_lin = 'd:/Work-Projects/MyProgs/ORB_data/NLED/LIN/'
root_nl = 'd:/Work-Projects/MyProgs/ORB_data/NLED/NL/'
dd_init = {
    'a0': 0.482,
    'R0': 1.62,
    'B0': 2.2,
    'mass_pf': 2 * constants.proton_mass
}

# --- STANDARD PROJECT ---
project_geometry        = dict(dd_init)
project_root_path       = root_lin
project_relative_path   = '/n0-egam/025b/f001-mpr'
project_name            = 'LINEAR\ ES\, ADIABATIC\ ELECTRONS'
b025_f001  = set_dd(project_geometry, project_root_path, project_relative_path, project_name)

# --- PROJECT WITH A LINEAR EGAM IN NLED-AUG WITH DRIFT-KINETIC ELECTRONS ---
b025k_f001 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-fpart/f001',
                                   'LINEAR\ EM\, DRIFT-KINETIC\ ELECTRONS')

# --- SOME NONLINEAR EGAM PROJECTS FROM NLED-AUG
nb025_f001 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f001', 'NL:\ ES\ EGAMb:\ f = 0.01,\ s_f = 0.50')
nb025_f009_orig = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f009-orig', 'NL:\ ES\ EGAMb:\ f = 0.0949')
nb025k_f001_dt5 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/KIN/f001', 'NL\ KIN:\ n_{EP}/n_e = 0.01')
