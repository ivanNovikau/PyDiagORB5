import Mix as mix
import ymath
import curve as crv
import ControlPlot as cpr
import numpy as np
import pyqtgraph.examples


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(cpr)


def test_curves():
    x = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    fx = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

    curves = crv.Curves().xlab('\kappa').ylab('\mathbf{\Phi}').tit('Some data')
    curves.flag_semilogy = True

    curves.new('c1').XS(x).YS(fx)\
        .leg('\mathbf{\Phi}').sty('--').col('blue')
    curves.new('c2').XS(x).YS(100*fx)\
        .leg('\mathbf{\overline{\Phi}}').sty('o').col('red')
    cpr.plot_curves(curves)

    curves.set_print().tit('The same data')
    curves.map('c2').col('grey').ms(16)
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


def test_animation_curves_2d():

    curves = crv.Curves().xlab('x').ylab('y')

    cpr.animation_curves_2d(curves)


def test_fft_gam_s1():
    T = 10
    w0 = 2 * np.pi / T
    print(w0)
    print(w0 / (2 * np.pi))
    x = np.linspace(0, 100, 201)
    y = np.cos(w0 * x)

    w, w2 = ymath.w_for_fft(x)
    f, f2 = ymath.fft_y(y)

    # initial singal
    curves = crv.Curves().xlab('x').ylab('y')
    curves.new('y').XS(x).YS(y).leg('init')
    cpr.plot_curves(curves)

    # fourier spectrum
    curves = crv.Curves().xlab('w').ylab('f')
    curves.new('f').XS(w).YS(f).leg('fft')
    cpr.plot_curves(curves)

    return


def test_fft_gam_2d():
    T = 10
    w0 = 2 * np.pi / T
    print(w0)
    print(w0 / (2 * np.pi))

    L = 0.5
    k0 = 2 * np.pi / L
    print(k0)
    print(k0 / (2 * np.pi))

    t = np.linspace(0, 100, 201)
    x = np.linspace(0, 2, 401)
    tt, xx = np.meshgrid(t, x)
    yy = np.cos(w0 * tt) * np.cos(k0 * xx)

    w, w2 = ymath.w_for_fft(t)
    ff, f2 = ymath.fft_y(yy, 1)

    # initial singal
    curves = crv.Curves().xlab('t').ylab('x').tit('yy')
    curves.new('yy').XS(t).YS(x).ZS(yy).leg('init')
    cpr.plot_curves_3d(curves)

    # fourier spectrum
    curves = crv.Curves().xlab('w').ylab('x').tit('ff')
    curves.new('f').XS(w).YS(x).ZS(ff).leg('ff')
    cpr.plot_curves_3d(curves)

    return


def test_pyqtgraph():
    pyqtgraph.examples.run()







