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
import fields3d
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
    mix.reload_module(fields3d)


def set_dd(dd_loc, root, path_loc, name_loc):
    proj_loc = {
        'path': root + path_loc,
        'name': name_loc,
    }
    dd_loc.update({'path': proj_loc['path'], 'project_name': proj_loc['name']})
    rd.init(dd_loc)
    return dd_loc


reload()
root_lin  = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/linear'
root_mu12 = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/'
root_ant  = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/antenna/'

dd_init = {
    'a0': 0.2475,
    'R0': 0.88,
    'B0': 1.4348,
    'mass_pf': 2 * constants.proton_mass
}

# --- LINEAR --
n80_lin  = set_dd(dict(dd_init), root_lin,  '/n80/', 'Linear:\ n80')

# --- NONLINEAR ---
n80           = set_dd(dict(dd_init), root_mu12, '/n80/', 'n80')
n80_potsc     = set_dd(dict(dd_init), root_mu12, '/n80-potsc/', 'n80')
n80_circular  = set_dd(dict(dd_init), root_mu12, '/n80-circular-chease/', 'n80-circular')

# --- ANTENNA ---
ant_test  = set_dd(dict(dd_init), root_ant,  '/TCV-test/',       'Antenna:\ test')
ant_n80   = set_dd(dict(dd_init),  root_ant, '/TCV-n80-start1/', 'Antenna:\ n=80')
ant_itpa_n20 = set_dd(dict(dd_init),  root_ant, '/ITPA-lin-n20-start/', 'Antenna:\ ITPA:\ n=20')
ant_itpa_n20_test = set_dd(dict(dd_init),  root_ant, '/ant-ITPA-n20/', 'Antenna:\ ITPA:\ n=20: test')
ant_itpa_n20_from_beginning = set_dd(dict(dd_init),  root_ant,
                                     '/launch-ant-ITPA-lin20/from-beginning',
                                     'Antenna:\ ITPA:\ n=20:\ from\ beginning')
ant_itpa_n20_from_beginning_cobra = set_dd(dict(dd_init),  root_ant,
                                     '/launch-ant-ITPA-lin20/from-beginning-cobra',
                                     'Antenna:\ ITPA:\ n=20:\ from\ beginning\ (cobra)')


ant_itpa_n20_start_1 = set_dd(dict(dd_init),  root_ant,
                            '/launch-ant-ITPA-lin20/start_1',
                                        'Antenna:\ ITPA:\ n=20:\ start-1')
ant_itpa_n20_start_2 = set_dd(dict(dd_init),  root_ant,
                            '/launch-ant-ITPA-lin20/start_2',
                                        'Antenna:\ ITPA:\ n=20:\ start-2')
ant_itpa_n20_start_2_ant = set_dd(dict(dd_init),  root_ant,
                            '/launch-ant-ITPA-lin20/start_2_ant',
                                        'Antenna:\ ITPA:\ n=20:\ start-2-ant')
ant_itpa_n20_start_2_nonfilt_ant = set_dd(dict(dd_init),  root_ant,
                            '/launch-ant-ITPA-lin20/start_2_nonfilt_ant',
                                        'Antenna:\ ITPA:\ n=20:\ start-2-nonfilt-ant')
ant_n20_n24_1 = set_dd(dict(dd_init),  root_ant,
                            '/two-modes/n20-n24-start-1',
                            'Antenna:\ n = 20, 24')


