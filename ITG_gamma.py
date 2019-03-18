import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)


def read_data(dd, path_to_folder):
    dd = rd.potsc(dd, path_to_folder)
    return dd


def plot_t(dd, s1, s2, chi1, oo={}):
    # coordinate axes
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    _, ids_s = mix.get_array(s, s1, s2)
    t, ids_t = mix.get_array_oo(oo, t, 't')
    id_chi, _ = mix.find(chi, chi1)

    # non-zonal (at phi = 0) Phi in chosen intervals
    pot_nz_chi = mix.get_slice(dd['potsc']['data'], ids_t, id_chi, ids_s)

    # averaging of the Phi on s
    pot_nz = np.mean(pot_nz_chi, axis=1)

    # plot data:
    cpr.plot_x1(t, pot_nz)


def plot_schi(dd, t1, oo={}):
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference
    r = dd['potsc']['r']  # create a new reference
    z = dd['potsc']['z']  # create a new reference

    id_t1 = np.argmax(t >= t1)
    pot_nz = np.array(dd['potsc']['data'][id_t1, :, :])  # actually copy data

    oop = {'xlabel': 's', 'ylabel': 'chi'}
    cpr.plot_x1x2(s, chi, pot_nz, oop)

    # oop = {'xlabel': 'r', 'ylabel': 'z'}
    # cpr.plot_x1x2(r, z, pot_nz, oop)


def calc_gamma_chi0(dd, s1, s2, chi1, oo={}):
    # coordinate axes
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    _, ids_s = mix.get_array(s, s1, s2)
    t, ids_t = mix.get_array_oo(oo, t, 't')
    id_chi,_ = mix.find(chi, chi1)

    # non-zonal (at phi = 0) Phi in chosen intervals
    pot_nz_chi = mix.get_slice(dd['potsc']['data'], ids_t, id_chi, ids_s)

    # averaging of the Phi on s
    pot_nz = np.mean(pot_nz_chi, axis=1)

    # estimation:
    wg_est = ymath.estimate_wg(t, pot_nz)

    # plotting: estimation: time evolution and peaks
    curves_est = crv.Curves().xlab('t').ylab('\Phi')
    curves_est.flag_semilogy = True
    curves_est.new('init')\
        .XS(t).YS(pot_nz).leg('init')
    curves_est.new('peaks')\
        .XS(wg_est['x_peaks']).YS(wg_est['y_peaks'])\
        .leg('peaks').sty('o').col('green')
    curves_est.new('fitting')\
        .XS(t).YS(wg_est['y_fit'])\
        .leg('fitting').col('red').sty('--')
    cpr.plot_curves(curves_est)

    # advanced w,g calculation
    ainf = {'est': wg_est,
            'x_start': wg_est['x_peaks'][0],
            'x_end': wg_est['x_peaks'][-1]
            }
    wg_adv = ymath.advanced_wg(t, pot_nz, ainf)
    curves_adv = crv.Curves().xlab('t').ylab('\Phi')
    curves_adv.flag_semilogy = True
    curves_adv.new('init') \
        .XS(t).YS(pot_nz).leg('init')
    curves_adv.new('fitting') \
        .XS(wg_adv['x_fit']).YS(wg_adv['y_fit'])\
        .leg('adv. fitting').col('red').sty('--')
    cpr.plot_curves(curves_adv)

    # print results:
    print('--- Estimation ---')
    print('w = {:0.3e}'.format(wg_est['w']))
    print('g = {:0.3e}'.format(wg_est['g']))

    print('--- Advanced ---')
    print('w = {:0.3e}'.format(wg_adv['w']))
    print('g = {:0.3e}'.format(wg_adv['g']))






