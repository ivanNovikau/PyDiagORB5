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
root_ref = 'd:/Work-Projects/MyProgs/ORB_data/MPR/GAMs/ES-EM/'

dd_init = {
    'a0': 0.5,
    'R0': 1.65,
    'B0': 2.0,
    'mass_pf': 2 * constants.proton_mass
}

# --- GAM --
q15 = set_dd(dict(dd_init), root_ref,  '/q15-mpr',  'ES:\ q = 1.5')
q3 = set_dd(dict(dd_init), root_ref,   '/q3-mpr',   'ES:\ q = 3.0')
q5 = set_dd(dict(dd_init), root_ref,   '/q5-mpr',   'ES:\ q = 5.0')
q15k = set_dd(dict(dd_init), root_ref, '/q15k-mpr', 'EM:\ q = 1.5')
q3k = set_dd(dict(dd_init), root_ref,  '/q3k-mpr',  'EM:\ q = 3.0')
q5k = set_dd(dict(dd_init), root_ref,  '/q5k-mpr',  'EM:\ q = 5.0')

# q15_dt20_ns64_N17  = set_dd(dict(dd_init), root_ref,  '/conv/q15-dt20-ns64-N1e7',
#                             'dt = 20, ns = 64, N=1e7')
# q15_dt10  = set_dd(dict(dd_init), root_ref,  '/conv/q15-dt10',
#                             'dt = 10')
# q15_dt15  = set_dd(dict(dd_init), root_ref,  '/conv/q15-dt15',
#                             'dt = 15')
# q15_N2e7  = set_dd(dict(dd_init), root_ref,  '/conv/q15-N2e7',
#                             'N = 2e7')
# q15_ns100_N3e7  = set_dd(dict(dd_init), root_ref,  '/conv/q15-ns100',
#                             'ns = 100, N = 3e7')
# q15_ns120_N3e7  = set_dd(dict(dd_init), root_ref,  '/conv/q15-ns120',
#                             'ns = 120, N = 3e7')
# q15_ns150_N45e7  = set_dd(dict(dd_init), root_ref,  '/conv/q15-ns150',
#                             'ns = 150, N = 45e7')
# q15_ns200_N45e7  = set_dd(dict(dd_init), root_ref,  '/conv/q15-ns200',
#                             'ns = 200, N = 6e7')

