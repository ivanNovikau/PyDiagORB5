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
# n0a      = set_dd(dict(dd_init), root_lin, '/scan-n/n0a', 'ES\ GAM:\ speak = 0.90')
# n0k      = set_dd(dict(dd_init), root_lin, '/scan-n/n0k', 'EM\ GAM:\ speak = 0.90')

# --- EGAM: LINEAR ---
b025_f001_mpr = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/f001-mpr',
                            'ES\ EGAMb:\ f = 0.01:\ v = 8.0')

s09_b025_f001_mpr = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/s09-f001',
                           'ES\ EGAMb:\ f = 0.01')

# kinetic electrons
b025k_f001 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/f001',
                                   'EM\ EGAMb:\ f = 0.01,\ m_i/m_e = 3676,\ \\beta_e = 2.7e-4')

b025k_f001_mie500  = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-mime/mie500',
                                   'EM\ EGAMb:\ f = 0.01,\ m_i/m_e = 500')
b025k_f001_mie1000 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-mime/mie1000',
                                   'EM\ EGAMb:\ f = 0.01,\ m_i/m_e = 1000')
b025k_f001_mie2000 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-mime/mie2000',
                                   'EM\ EGAMb:\ f = 0.01,\ m_i/m_e = 2000')

# scan on v_parallel
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
b025_f0004 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f0004-mpr', 'f = 0.004,\ v = 8.0,\ T = 1.0')
b025_f0006 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f0006-mpr', 'ES\ EGAMb:\ f = 0.006')
b025_f0008 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f0008-mpr', 'ES\ EGAMb:\ f = 0.008')
b025_f002 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f002', 'ES\ EGAMb:\ f = 0.02')
b025_f003 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f003', 'ES\ EGAMb:\ f = 0.03')
b025_f004 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f004', 'ES\ EGAMb:\ f = 0.04')
b025_f005 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f005', 'ES\ EGAMb:\ f = 0.05')
b025_f007 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f007', 'ES\ EGAMb:\ f = 0.07')
b025_f009 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f009', 'ES\ EGAMb:\ f = 0.0949')

# scan on EP localisation:
b035_f001 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/sf035-f001-mpr', 'ES\ EGAMb:\ \\rho_f = 0.35')

b045_f005_v9  = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/sf045-f005-vp9',
                      'ES\ EGAMb:\ \\rho_f = 0.45,\ v_{\parallel} = 9.0,\ f = 0.05')
b045_f005_v9_dt10  = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/sf045-f005-vp9-dt10',
                      'ES\ EGAMb:\ \\rho_f = 0.45,\ v_{\parallel} = 9.0,\ f = 0.05,\ dt = 10')
b045_f005_v10 = set_dd(dict(dd_init), root_lin,
                    '/n0-egam/025b/sf045-f005-vp10',
                    'ES\ EGAMb:\ \\rho_f = 0.45,\ v_{\parallel} = 10.0,\ f = 0.05')
b045_f005_v10_dt10 = set_dd(dict(dd_init), root_lin,
                    '/n0-egam/025b/sf045-f005-vp10-dt10',
                    'ES\ EGAMb:\ \\rho_f = 0.45,\ v_{\parallel} = 10.0,\ f = 0.05,\ dt = 10')

# scan on EP temperature:
b025_v6_T1 = set_dd(dict(dd_init),  root_lin,
                       '/n0-egam/025b/v60-f001-T1',  'ES\ EGAMb:\ v = 6.0,\ T = 1.0')
b025_v6_T08 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/v60-f001-T08', 'ES\ EGAMb:\ v = 6.0,\ T = 0.8')
b025_v6_T06 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/v60-f001-T06', 'ES\ EGAMb:\ v = 6.0,\ T = 0.6')
b025_v6_T04 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/v60-f001-T04', 'ES\ EGAMb:\ v = 6.0,\ T = 0.4')

b025_v3_T01 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-v30-T010',
                     'f = 0.0949,\ v = 3.0,\ T = 0.10')
b025_v3_T015 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-v30-T015',
                     'f = 0.0949,\ v = 3.0,\ T = 0.15')
b025_v3_T02 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-v30-T02',
                     'f = 0.0949,\ v = 3.0,\ T = 0.20')

b025_T02_v28 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-T02-v28',
                     'f = 0.0949,\ T = 0.20,\ v = 2.8')
b025_T02_v32 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-T02-v32',
                     'f = 0.0949,\ T = 0.20,\ v = 3.2')

b025_v35_T015 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-v35-T015',
                     'f = 0.0949,\ v = 3.5,\ T = 0.15')
b025_v35_T020 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-T02-v35',
                     'f = 0.0949,\ v = 3.5,\ T = 0.20')
b025_v35_T022 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-v35-T022',
                     'f = 0.0949,\ v = 3.5,\ T = 0.22')
b025_v35_T025 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-v35-T025',
                     'f = 0.0949,\ v = 3.5,\ T = 0.25')

b025_v80_T04 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f001-v80-T04',
                     'f = 0.01,\ v = 8.0,\ T = 0.4')
b025_v80_T06 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f001-v80-T06',
                     'f = 0.01,\ v = 8.0,\ T = 0.6')
b025_v80_T08 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f001-v80-T08',
                     'f = 0.01,\ v = 8.0,\ T = 0.8')

# --- LINEAR TWO BEAMS ---
b025_b2_T15 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/TWO-BEAMS/v80-v35-T015',
                     'B2:\ v = 3.5,\ T = 0.15')
b025_b2_T25 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/TWO-BEAMS/v80-v35-T025',
                     'B2:\ v = 3.5,\ T = 0.25')


# --- EGAM: NON-LINEAR ---
nb025_f0004 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f0004',
                        'NL:\ ES\ EGAMb:\ f = 0.004,\ s_f = 0.50'
                        )
nb025_f001 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f001-mpr',
                        'NL:\ ES\ EGAMb:\ f = 0.01,\ s_f = 0.50'
                        )
nb025_f001_wp = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f001-mpr-wp',
                        'NL:\ ES\ EGAMb:\ f = 0.01,\ s_f = 0.50,\ WPN'
                        )
nb025_f002 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f002',
                        'NL:\ ES\ EGAMb:\ f = 0.02,\ s_f = 0.50'
                        )
nb025_f003 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f003',
                        'NL:\ ES\ EGAMb:\ f = 0.03,\ s_f = 0.50'
                        )
nb025_f004 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f004',
                        'NL:\ ES\ EGAMb:\ f = 0.04,\ s_f = 0.50'
                        )
nb025_f005 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f005',
                        'NL:\ ES\ EGAMb:\ f = 0.05,\ s_f = 0.50'
                        )
nb025_f009 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f009',
                        'NL:\ ES\ EGAMb:\ f = 0.0949,\ s_f = 0.50'
                        )
nb025_f001_v7 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/v70-f001',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 7.0,\ f = 0.01,\ s_f = 0.50'
                        )
nb025_f001_v75 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/v75-f001',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 7.5,\ f = 0.01,\ s_f = 0.50'
                        )
nb025_f001_v9 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/v90-f001',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 9.0,\ f = 0.01,\ s_f = 0.50'
                        )
nb025_f009_v3 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f009-v30-T02',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 3.0,\ f = 0.095,\ T = 0.20'
                        )
nb025_f009_v35 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f009-v35-T02',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 3.5,\ f = 0.095,\ T = 0.20'
                        )

nb025_v60_T04 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v60/T04',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 6.0,\ f = 0.01,\ T = 0.40'
                        )
nb025_v60_T06 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v60/T06',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 6.0,\ f = 0.01,\ T = 0.60'
                        )
nb025_v60_T08 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v60/T08',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 6.0,\ f = 0.01,\ T = 0.80'
                        )
nb025_v60_T10 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v60/T10',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 6.0,\ f = 0.01,\ T = 1.0'
                        )

nb025_b2_T15 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/TWO-BEAMS/v80-v35-T015',
                        'NL:\ two\ beams:\ T = 0.15'
                        )
nb025_b2_T25 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/TWO-BEAMS/v80-v35-T025',
                        'NL:\ two\ beams:\ T = 0.25'
                        )



