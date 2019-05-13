import Mix as mix
import read_data as rd
import zf_gam as zf
import equil_profiles as equ
import ControlPlot as cpr
import transport
import work_profiles
import ITG_gamma as itg
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


# region --- NL SIMULATIONS ---
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/ref-results-43516/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/n40/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/n60/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/n80/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/n80-potsc/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/n128/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/adhoc-n80/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT6-mun3/n80'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/EGAM/n80-v6-s35'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/EGAM/n80-v6-s45'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/n80-modT/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT6-mun5/n80/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/scan-nptot-n80/n1e8/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/scan-nptot-n80/n3e8/'
# endregion

# region --- LINEAR GAM ---
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-GAM/n0/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-GAM/n0-Krook/'
# endregion

# region --- LINEAR ITG ---
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n40'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n45'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n50'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n55'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n58'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n60'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n62'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n65'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n75'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n78'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n80'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n82'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n95'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n98'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n100'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/more-precise/n102'

# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT12-mun3/n40'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT12-mun3/n50'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT12-mun3/n60'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT12-mun3/n70'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT12-mun3/n75'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT12-mun3/n80'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT12-mun3/n90'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT12-mun3/n95'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT12-mun3/n100'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT12-mun3/n110'

# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT6-mun5/n10'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT6-mun5/n20'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT6-mun5/n30'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT6-mun5/n40'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT6-mun5/n50'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT6-mun5/n60'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT6-mun5/n75'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT6-mun5/n80'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#            'turbulence/TCV/linear/scan-ITG/muT6-mun5/n95'

# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT6-mun3/n50'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT6-mun3/n60'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT6-mun3/n75'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT6-mun3/n80'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT6-mun3/n95'

# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT9-mun5/n50'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT9-mun5/n60'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT9-mun5/n75'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT9-mun5/n80'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/' \
#             'turbulence/TCV/linear/scan-ITG/muT9-mun5/n95'

# endregion

# region --- LINEAR EGAM ---
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/flat-v4'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/flat-v8'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-v4-f01'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-v8-f01'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#             'TCV/linear/scan-EGAM/ge-v65-f02'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/gc-v4-f02'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-v2-f02'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-v4-f02'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-vf4-T2'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-vf4-T4'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#            'TCV/linear/scan-EGAM/ge-v6-f02'                    # <-- !!!  s_gauss ~ 0.7, s_egam ~ 0.85,
#                                                                # w[cs/a] ~ 0.64, g[cs/a] ~ 5.6e-3
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#             'TCV/linear/scan-EGAM/gc-v6-f02'                     # <-- !!! s_gauss ~ 0.54, s_egam ~ 0.68,
#                                                                  # w[cs/a] ~ 0.93, g[cs/a] ~ 1.8e-2
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#             'TCV/linear/scan-EGAM/ge-v55-f02'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/scan-position-f02-v55/s4'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#             'TCV/linear/scan-EGAM/scan-position-f02-v55/s5'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#             'TCV/linear/scan-EGAM/scan-position-f02-v55/s6'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#             'TCV/linear/scan-EGAM/scan-position-f02-v55/s7'  # !!! s_egam ~ 0.86,
#                                                              # w[cs/a] ~ 0.66, g[cs/a] ~ 6.8e-3
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#             'TCV/linear/scan-EGAM/scan-position-f02-v55/s8'  # !!! s_egam ~ 0.95,
#                                                              # w[cs/a] ~ 0.5, g[cs/a] ~ 1.0e-2
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#             'TCV/linear/scan-EGAM/scan-position-f02-v55/s9'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#             'TCV/linear/scan-EGAM/scan-position-f02-v55/s8-ds2'  # !!! s_egam ~ 0.93,
#                                                                  # w[cs/a] ~ 0.45, g[cs/a] ~ 1.3e-2
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/scan-position-f02-v6/s4'  # very spread s_egam_1 ~ 0.48, s_egam_2 ~ [0.6, 1.0]
#                                                              # w[cs/a] ~ 1.0, g[cs/a] ~ 8e-3
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/scan-position-f02-v6/s6'  # s_egam ~ 0.68, [0.74, 1.0]
#                                                              # w[cs/a] ~ 0.87, g[cs/a] ~ 8.7e-3
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/scan-position-f02-v6/s7'  # s_egam ~ 0.77, [0.83, 1.0]
#                                                              # w[cs/a] ~ 0.72, g[cs/a] ~ 6e-3
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/scan-position-f02-v6/s8'  # s_egam ~ 0.86, 0.95
#                                                              # w[cs/a] ~ 0.54, g[cs/a] ~ 1.5e-3
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s7/f005/'  # s_egam ~ [0.7, 0.8]
#                                                                      # w[cs/a] ~ 0.78, g[cs/a] ~ 9.3e-3
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s7/f01/'  # s_egam ~ [0.7, 0.8], [0.8, 1.0]
#                                                                     # w[cs/a] ~ 0.72, g[cs/a] ~ 1.8e-2
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s7/f015/'  # s_egam ~ [0.7, 0.8], [0.8, 1.0]
#                                                                      # w[cs/a] ~ 0.69, g[cs/a] ~ 1.3e-2
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s8/f005/'   # s_egam ~ [0.8, 1.0]
#                                                                       # w[cs/a] ~ 0.59, g[cs/a] ~ 1.9e-2
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s8/f01/'  # s_egam ~ [0.8, 1.0]
#                                                                     # w[cs/a] ~ 0.54, g[cs/a] ~ 2.4e-2
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s8/f015/'  # s_egam ~ [0.8, 1.0]
#                                                                      # w[cs/a] ~ 0.52, g[cs/a] ~ 1.8e-2
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/scan-position-f02-v55/s75-ds25/'  # very close to the marginality
#                                                                      # (damping ~ growth)
#                                                                      # s_egam ~ [0.85, 1.0]
#                                                                      # w[cs/a] ~ 0.48, g[cs/a] ~ [2.4, 4.4]e-3
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/scan-position-f02-v55/s7-ds3/'    # gam-like oscillations periodically appear
#                                                                      # w[cs/a] ~ 0.55
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/scan-position-f02-v6/s4-ds6/'   # no EGAMs
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#               'TCV/linear/scan-EGAM/scan-v-s4/v58/'                  # s_egam_1 ~ [0.4, 1.0]
#                                                                      # w[cs/a] ~ 1.0, g[cs/a] ~ 3.5e-3
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#               'TCV/linear/scan-EGAM/scan-v-s4/v62/'                 # s_egam_1 ~ [0.4, 1.0]
#                                                                     # w[cs/a] ~ 1.1, g[cs/a] ~ 1.25e-2
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#               'TCV/linear/scan-EGAM/scan-position-f02-v6/s35/'      # s_egam_1 ~ [0.2, 1.0]
#                                                                     # w[cs/a] ~ 1.18, g[cs/a] ~ 2.5e-2
#                                                                     # at s = 0.95, there several frequencies in FFT
#                                                                     # including ~0.35, ~0.44 ... ~1.2
#                                                                     # very small ZF contribution
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#               'TCV/linear/scan-EGAM/scan-position-f02-v6/s45/'      # s_egam_1 ~ [0.45, 1.0]
#                                                                     # w[cs/a] ~ 1.048, g[cs/a] ~ 9.0e-3
#                                                                     # high ZF contribution (specially, at the center)
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v6-s35/f001/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v6-s35/f005/'

# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s8/f001/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s8/f003/'       # w[wci] ~ 4.7e-03, g[wci] ~ 4e-05
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s8/f005-5e7/' # -> the same as with nptot = 1e8
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s8/f002/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s8/f0024/'
# path_TCV = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/' \
#              'TCV/linear/scan-EGAM/concentration-scan-v55-s8/f0026/'    # very close to marginality from GAM side

# endregion

# region --- Project structure initialization ---
dd = {
    'path': path_TCV,
    'a0': 0.2475,
    'R0': 0.88,
    'B0': 1.4348,
    'mass_pf': 2 * constants.proton_mass
}
# endregion

# region --- Initialization of the simulation's data ---
reload()
rd.init(dd)

# endregion

