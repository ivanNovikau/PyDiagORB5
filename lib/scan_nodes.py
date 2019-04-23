import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np
from scipy import interpolate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)


def scan_tcv_n40():
    n_case = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])

    nphi = np.array([160, 160, 160,  192, 256, 160,  160, 192, 160,
                        192, 160, 160,  192, 256])
    nchi = np.array([320, 320, 320,  384, 512, 320,  320, 384, 320,
                        384, 320, 320,  384, 512])

    nclones = np.array([1, 1, 1,  1, 3, 1,  3, 1, 1,  1, 1, 3,  1, 3])
    ncpn = np.array([32, 32, 32,  48, 48, 32,  48, 48, 32,
                      48, 32, 48,  48, 48])
    nnodes = np.array([5, 10, 15,  16, 16, 20,  20, 24, 25,
                      28, 30, 30,  32, 32])

    ntime = np.array([1.9882360E+02, 1.0559034E+02, 7.2066126E+01,
                      6.4477990E+01, 6.3376985E+01, 5.7538552E+01,
                      4.6553262E+01, 4.2268782E+01, 4.5905534E+01,
                      3.8858627E+01, 3.9604043E+01, 3.2844944E+01,
                      3.3203459E+01, 3.5985387E+01])

    ncpu = ncpn * nnodes
    ncpu_time = ncpu * ntime
    effectiveness = 1 / ncpu_time

    curves = crv.Curves().xlab('case').ylab('norm.\ values').tit('_')\
        .xt(n_case)
    curves.flag_norm = True
    curves.new('y1').XS(n_case).YS(effectiveness)\
        .leg('effectiveness: 1/(ncpu\cdot time)').sty('o-')
    curves.new('y2').XS(n_case).YS(ntime).leg('time').sty('o:')
    cpr.plot_curves(curves)