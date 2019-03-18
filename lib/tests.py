import Mix as mix
import ymath
import curve as crv
import ControlPlot as cpr
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(cpr)


def test_curves():
    x = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    fx = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

    curves = crv.Curves().xlab('\kappa').ylab('\mathbf{\Phi}')
    curves.flag_semilogy = True

    curves.new('c1').XS(x).YS(fx).sty('--').col('blue')\
        .leg('\Phi')
    curves.new('c2').XS(x).YS(100*fx).sty('o').col((1, 0, 0))\
        .leg('\overline{\Phi}')

    cpr.plot_curves(curves)

    curves.set_print()
    curves.map('c2').col('green')
    cpr.plot_curves(curves)


def test_advanced_wg():
    A0 = 2
    w = 0.5
    g = 0.05

    x = np.linspace(0.0, 100.0, 1001)
    y = A0 * np.cos(w * x) * np.exp(g * x)

    # estimation
    ainf = {'x_start': 60, 'x_end': 80}
    wg_est = ymath.estimate_wg(x, y, ainf)

    curves = crv.Curves().xlab('x').ylab('y')
    curves.flag_semilogy = True
    curves.new('init').XS(x).YS(y).leg('init')
    curves.new('peaks') \
        .XS(wg_est['x_peaks']).YS(wg_est['y_peaks']) \
        .leg('peaks').sty('o').col('green')
    curves.new('fitting') \
        .XS(wg_est['x_fit']).YS(wg_est['y_fit'])\
        .leg('est: fitting').col('red').sty('--')
    cpr.plot_curves(curves)

    # advanced fitting
    ainf = {'est': wg_est,
            'x_start': wg_est['x_peaks'][0],
            'x_end': wg_est['x_peaks'][-1]
            }
    wg_adv = ymath.advanced_wg(x, y, ainf)

    curves_adv = crv.Curves().xlab('t').ylab('\Phi')
    curves_adv.flag_semilogy = True
    curves_adv.new('init') \
        .XS(x).YS(y).leg('init')
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



