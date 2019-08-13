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
root_ref = 'd:/Work-Projects/MyProgs/ORB_data/MPR/'
root_lin = 'd:/Work-Projects/MyProgs/ORB_data/NLED/LIN/'
root_nl = 'd:/Work-Projects/MyProgs/ORB_data/NLED/NL/'

dd_init = {
    'a0': 0.482,
    'R0': 1.62,
    'B0': 2.2,
    'mass_pf': 2 * constants.proton_mass
}

# # --- REFERENCE EM LINEAR EGAM SIMULATION ---
# dd_egam_ref = set_dd(dict(dd_init), root_ref,
#                  'nled/kin/linear/REF-WORK-CPS2019/', 'REF')

# # --- GAM --
# n0a      = set_dd(dict(dd_init), root_lin, '/scan-n/n0a',         'ES\ GAM:\ speak = 0.90')
# n0k      = set_dd(dict(dd_init), root_lin, '/scan-n/n0k',     'EM\ GAM:\ speak = 0.90')

# --- EGAM: LINEAR ---
b025_f001_mpr = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/f001-mpr',
                            'ES\ EGAMb:\ f = 0.01:\ v = 8.0')

s09_b025_f001_mpr = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/s09-f001',
                           'ES\ EGAMb:\ f = 0.01')

k_s09_b025_f001_dt5_ns256 = set_dd(dict(dd_init), root_lin,
                                   '/n0-egam/025b/conv-k/smax-90/dt5-ns256',
                                   'EM\ EGAMb:\ f = 0.01')

# # negative shifted Maxwellian
# mm025_f001 = set_dd(dict(dd_init), root_lin, '/n0-egam/025m/mvpar-f001-mpr',
#                     'ES\ EGAMmm:\ f = 0.01')

# scan on v_parallel
# b025_f001_v10 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/v10-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 1.0')  # no EGAMs
# b025_f001_v20 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/v20-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 2.0')  # no EGAMs
# b025_f001_v30 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/v30-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 3.0')  # no EGAMs
# b025_f001_v35 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/v35-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 3.5')  # no EGAMs
# b025_f001_v40 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/v40-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 4.0')  # no EGAMs
# b025_f001_v50 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/v50-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 5.0')  # no EGAMs
#
# b025_f009_v35 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/v35-f009', 'ES\ EGAMb:\ f = 0.0949,\ v = 3.5') # no EGAMs
# b025_f009_v40 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/v40-f009', 'ES\ EGAMb:\ f = 0.0949,\ v = 4.0') # no EGAMs

b025_f001_v70 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/v70-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 7.0')
b025_f001_v75 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/v75-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 7.5')
b025_f001_v78 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/v78-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 7.8')
b025_f001_v82 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/v82-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 8.2')
b025_f001_v85 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/v85-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 8.5')
b025_f001_v90 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/v90-f001', 'ES\ EGAMb:\ f = 0.01,\ v = 9.0')

# scan on concentration
b025_f0002 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f0002-mpr', 'ES\ EGAMb:\ f = 0.002')  # no EGAMs
b025_f0004 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f0004-mpr', 'ES\ EGAMb:\ f = 0.004')
b025_f0006 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f0006-mpr', 'ES\ EGAMb:\ f = 0.006')
b025_f0008 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f0008-mpr', 'ES\ EGAMb:\ f = 0.008')
b025_f002 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f002', 'ES\ EGAMb:\ f = 0.02')
b025_f005 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f005', 'ES\ EGAMb:\ f = 0.05')
b025_f007 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f007', 'ES\ EGAMb:\ f = 0.07')
b025_f009 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f009', 'ES\ EGAMb:\ f = 0.0949')

# scan on EP localisation:
b035_f001 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/sf035-f001-mpr', 'ES\ EGAMb:\ \\rho_f = 0.35')
# b045_f001 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/sf045-f001-mpr', 'ES\ EGAMb:\ \\rho_f = 0.45')  # no EGAM
# b050_f001 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/sf050-f001-mpr', 'ES\ EGAMb:\ \\rho_f = 0.50')  # no EGAM
# b055_f001 = set_dd(dict(dd_init), root_lin,
#                        '/n0-egam/025b/sf055-f001-mpr', 'ES\ EGAMb:\ \\rho_f = 0.55')  # no EGAM

b045_f005_v9  = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/sf045-f005-vp9',
                      'ES\ EGAMb:\ \\rho_f = 0.45,\ v_{\parallel} = 9.0,\ f = 0.05')
b045_f005_v10 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/sf045-f005-vp10',
                      'ES\ EGAMb:\ \\rho_f = 0.45,\ v_{\parallel} = 10.0,\ f = 0.05')


# --- EGAM: NON-LINEAR ---
nb025_f001 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f001-mpr',
                        'NL:\ ES\ EGAMb:\ f = 0.01,\ s_f = 0.50'
                        )
nb025_f001_wp = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f001-mpr-wp',
                        'NL:\ ES\ EGAMb:\ f = 0.01,\ s_f = 0.50,\ WPN'
                        )

nb025_f0004 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f0004',
                        'NL:\ ES\ EGAMb:\ f = 0.004,\ s_f = 0.50'
                        )
nb025_f005 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f005',
                        'NL:\ ES\ EGAMb:\ f = 0.05,\ s_f = 0.50'
                        )
nb025_f001_v7 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/v70-f001',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 7.0,\ f = 0.01,\ s_f = 0.50'
                        )
nb025_f001_v9 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/v90-f001',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 9.0,\ f = 0.01,\ s_f = 0.50'
                        )



