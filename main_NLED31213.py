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

# --- EGAM: LINEAR ---
# KIN. ELE. (beta_e = 2.7e-4): Realistic electrons: scan on f_part
b025k_f0006 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-fpart/f0006',
                                   'EM\ EGAMb:\ n_{EP}/n_e = 0.006,\ T_{EP} = 1.0, v_{\parallel} = 8.0')
b025k_f001 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-fpart/f001',
                                   'EM\ EGAMb:\ n_{EP}/n_e = 0.01,\ T_{EP} = 1.0, v_{\parallel} = 8.0')
b025k_f002 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-fpart/f002',
                                   'EM\ EGAMb:\ n_{EP}/n_e = 0.02,\ T_{EP} = 1.0, v_{\parallel} = 8.0')
b025k_f003 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-fpart/f003',
                                   'EM\ EGAMb:\ n_{EP}/n_e = 0.03,\ T_{EP} = 1.0, v_{\parallel} = 8.0')
b025k_f005 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-fpart/f005',
                                   'EM\ EGAMb:\ n_{EP}/n_e = 0.05,\ T_{EP} = 1.0, v_{\parallel} = 8.0')
b025k_f007 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-fpart/f007',
                                   'EM\ EGAMb:\ n_{EP}/n_e = 0.07,\ T_{EP} = 1.0, v_{\parallel} = 8.0')
b025k_f009 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-fpart/f009',
                                   'EM\ EGAMb:\ n_{EP}/n_e = 0.09,\ T_{EP} = 1.0, v_{\parallel} = 8.0')
b025k_f009_dt3 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-fpart/f009-dt3',
                                   'EM\ EGAMb:\ n_{EP}/n_e = 0.09,\ dt = 3')
b025k_f009_ne24 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/scan-fpart/f009-ne24',
                                   'EM\ EGAMb:\ n_{EP}/n_e = 0.09,\ Ne = 2.4e8')

# KIN. ELE. (beta_e = 2.7e-4): Realistic electrons: different v
b025k_v35_T015 = set_dd(dict(dd_init), root_lin, '/n0-egam/025b/KIN-smax90/f009-v35-T015/',
                                   'EM\ EGAMb:\ n_{EP}/n_e = 0.09,\ T_{EP} = 0.15, v_{\parallel} = 3.5')

# scan on concentration: v = 8.0, T = 1.0
b025_f0004 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f0004-mpr', 'f = 0.004,\ v = 8.0,\ T = 1.0')
b025_f0006 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f0006-mpr', 'ES\ EGAMb:\ f = 0.006')
b025_f0008 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f0008-mpr', 'ES\ EGAMb:\ f = 0.008')
b025_f001 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/f001-mpr', 'ES\ EGAMb:\ f = 0.01')
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

# scan on Gaussian width:
b025_f001_w005 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/scan-sigma-f001/sigma-005',
                        'n_{EP}/n_e = 0.01,\ \sigma = 0.05')
b025_f001_w015 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/scan-sigma-f001/sigma-015',
                        'n_{EP}/n_e = 0.01,\ \sigma = 0.15')
b025_f001_w020 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/scan-sigma-f001/sigma-020',
                        'n_{EP}/n_e = 0.01,\ \sigma = 0.20')

b025_f009_w005 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/scan-sigma-f009/sigma-005',
                        'n_{EP}/n_e = 0.09,\ \sigma = 0.05')
b025_f009_w015 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/scan-sigma-f009/sigma-015',
                        'n_{EP}/n_e = 0.09,\ \sigma = 0.15')
b025_f009_w020 = set_dd(dict(dd_init), root_lin,
                       '/n0-egam/025b/scan-sigma-f009/sigma-020',
                        'n_{EP}/n_e = 0.09,\ \sigma = 0.20')

# # scan on EP temperature:
b025_v35_T015 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-v35-T015',
                     'f = 0.0949,\ v = 3.5,\ T = 0.15')
b025_v35_T020 = set_dd(dict(dd_init), root_lin,
                     '/n0-egam/025b/f009-v35-T020',
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

# --- EGAM: NON-LINEAR ---
nb025_f0004 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f0004',
                        'NL:\ ES\ EGAMb:\ f = 0.004,\ s_f = 0.50'
                        )
nb025_f001 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f001',
                        'NL:\ ES\ EGAMb:\ f = 0.01,\ s_f = 0.50'
                        )
nb025_f001_x2 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f001-ns2-N4-dtd2', 'NL\ ES:\ n_{EP}/n_e = 0.01, x2')
# nb025_f001_wp = set_dd(dict(dd_init), root_nl,
#                         '/n0-egam/025b/f001-fwp',
#                         'NL:\ ES\ EGAMb:\ f = 0.01,\ s_f = 0.50,\ WPN'
#                         )
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
# nb025_f009 = set_dd(dict(dd_init), root_nl,
#                         '/n0-egam/025b/f009',
#                         'NL:\ ES\ EGAMb:\ f = 0.0949,\ s_f = 0.50'
#                         )
# nb025_f009_x2 = set_dd(dict(dd_init), root_nl,
#                         '/n0-egam/025b/f009-x2',
#                         'NL:\ ES\ EGAMb:\ f = 0.0949,\ s_f = 0.50, x2'
#                         )
nb025_f009_orig = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/f009-orig',
                        'NL:\ ES\ EGAMb:\ f = 0.0949'
                        )
# nb025_f009_wp = set_dd(dict(dd_init), root_nl,
#                         '/n0-egam/025b/f009-wp',
#                         'NL:\ ES\ EGAMb:\ f = 0.0949,\ s_f = 0.50,\ WPN'
#                         )
# nb025_f001_v7 = set_dd(dict(dd_init), root_nl,
#                         '/n0-egam/025b/v70-f001',
#                         'NL:\ ES\ EGAMb:\ v_{\parallel} = 7.0,\ f = 0.01,\ s_f = 0.50'
#                         )
# nb025_f001_v75 = set_dd(dict(dd_init), root_nl,
#                         '/n0-egam/025b/v75-f001',
#                         'NL:\ ES\ EGAMb:\ v_{\parallel} = 7.5,\ f = 0.01,\ s_f = 0.50'
#                         )
# nb025_f001_v9 = set_dd(dict(dd_init), root_nl,
#                         '/n0-egam/025b/v90-f001',
#                         'NL:\ ES\ EGAMb:\ v_{\parallel} = 9.0,\ f = 0.01,\ s_f = 0.50'
#                         )

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
nb025_v60_T04_x2 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v60/T04-x2',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 6.0,\ f = 0.01,\ T = 0.40,\ x2'
                        )

nb025_v35_T015 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v35/T015',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 3.5,\ f = 0.095,\ T = 0.15'
                        )
nb025_v35_T020 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v35/T020',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 3.5,\ f = 0.095,\ T = 0.20'
                        )
nb025_v35_T022 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v35/T022',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 3.5,\ f = 0.095,\ T = 0.22'
                        )
nb025_v35_T025 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v35/T025',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 3.5,\ f = 0.095,\ T = 0.25'
                        )
nb025_v35_T015_x2 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v35/T015-x2',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 3.5,\ f = 0.095,\ T = 0.15,\ x2'
                        )

nb025_v35_T025_bwn = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v35/T025-bulk-wpn',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 3.5,\ f = 0.095,\ T = 0.25:\ wave-Bulk\ NL'
                        )
nb025_v35_T025_fwn = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v35/T025-fast-wpn',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 3.5,\ f = 0.095,\ T = 0.25:\ wave-Fast\ NL'
                        )
nb025_v35_T015_fwn = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/scan-T-v35/T015-fast-wpn',
                        'NL:\ ES\ EGAMb:\ v_{\parallel} = 3.5,\ f = 0.095,\ T = 0.15:\ EGAM-EP\ NL'
                        )

# KINETIC ELECTRONS
nb025k_f001_dt5 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/KIN/f001', 'NL\ KIN:\ n_{EP}/n_e = 0.01')
nb025k_f001_dt3 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/KIN/f001-dt3', 'NL\ KIN:\ n_{EP}/n_e = 0.01, dt = 3')
nb025k_f001_dt1 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/KIN/f001-dt1', 'NL\ KIN:\ n_{EP}/n_e = 0.01, dt = 1')
nb025k_f009_mie500_dt3 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/KIN/f009-mie500-dt3', 'NL\ KIN:\ n_{EP}/n_e = 0.095, mie = 500, dt = 3')
nb025k_f001_x2 = set_dd(dict(dd_init), root_nl,
                        '/n0-egam/025b/KIN/f001-ns2-N4', 'NL\ KIN:\ n_{EP}/n_e = 0.01, x2')



