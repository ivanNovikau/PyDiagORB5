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
root_LA = 'd:/Work-Projects/MyProgs/ORB_data/MPR/nled/LIN/adiab/'
root_LK = 'd:/Work-Projects/MyProgs/ORB_data/MPR/nled/LIN/kin/'

root_new = 'd:/Work-Projects/MyProgs/ORB_data/NLED/LIN/'

dd_init = {
    'a0': 0.482,
    'R0': 1.62,
    'B0': 2.2,
    'mass_pf': 2 * constants.proton_mass
}

# --- REFERENCE EM LINEAR EGAM SIMULATION ---
dd_egam_ref = set_dd(dict(dd_init), root_ref,
                 'nled/kin/linear/REF-WORK-CPS2019/', 'REF')

# --- GAM --
n0a_speak00  = set_dd(dict(dd_init), root_LA,  '/scan-n/speak00/n0/', 'ES\ GAM:\ speak = 0.0')
n0a      = set_dd(dict(dd_init), root_new, '/scan-n/n0a',         'ES\ GAM:\ speak = 0.90')
n0k      = set_dd(dict(dd_init), root_new, '/scan-n/n0k',     'EM\ GAM:\ speak = 0.90')

# --- EGAM: corrected density profiles ---
# b025_dt20 = set_dd(dict(dd_init), root_new, '/n0-egam/025b-dt20', 'ES\ EGAMb: dt = 20')
# b025_dt10 = set_dd(dict(dd_init), root_new, '/n0-egam/025b-dt10', 'ES\ EGAMb: dt = 10')
# b025_dt5 = set_dd(dict(dd_init), root_new,  '/n0-egam/025b-dt5',  'ES\ EGAMb: dt = 5')
# m025_dt20 = set_dd(dict(dd_init), root_new, '/n0-egam/025m-dt20', 'ES\ EGAMm: dt = 20')
# m025_dt10 = set_dd(dict(dd_init), root_new, '/n0-egam/025m-dt10', 'ES\ EGAMm: dt = 10')
# m025_dt5 = set_dd(dict(dd_init), root_new,  '/n0-egam/025m-dt5',  'ES\ EGAMm: dt = 5')

b025_f001_dt10 = set_dd(dict(dd_init), root_new, '/n0-egam/025b-f001-dt10', 'ES\ EGAMb: f = 0.01, dt = 10')
b025_f001_N37  = set_dd(dict(dd_init), root_new, '/n0-egam/025b-f001-N37', 'ES\ EGAMb: f = 0.01, N = 3e7')

b025_f009 = set_dd(dict(dd_init), root_new, '/n0-egam/025b-dt20', 'ES\ EGAMb: f = 0.0949')
b025_f007 = set_dd(dict(dd_init), root_new, '/n0-egam/025b-f007', 'ES\ EGAMb: f = 0.07')
b025_f005 = set_dd(dict(dd_init), root_new, '/n0-egam/025b-f005', 'ES\ EGAMb: f = 0.05')
b025_f002 = set_dd(dict(dd_init), root_new, '/n0-egam/025b-f002', 'ES\ EGAMb: f = 0.02')
b025_f001 = set_dd(dict(dd_init), root_new, '/n0-egam/025b-f001', 'ES\ EGAMb: f = 0.01')

m025_f009 = set_dd(dict(dd_init), root_new, '/n0-egam/025m-dt20', 'ES\ EGAMm: f = 0.0949')
m025_f007 = set_dd(dict(dd_init), root_new, '/n0-egam/025m-f007', 'ES\ EGAMm: f = 0.07')
m025_f005 = set_dd(dict(dd_init), root_new, '/n0-egam/025m-f005', 'ES\ EGAMm: f = 0.05')
m025_f002 = set_dd(dict(dd_init), root_new, '/n0-egam/025m-f002', 'ES\ EGAMm: f = 0.02')
m025_f001 = set_dd(dict(dd_init), root_new, '/n0-egam/025m-f001', 'ES\ EGAMm: f = 0.01')


