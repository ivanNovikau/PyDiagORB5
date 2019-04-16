import Mix as mix
import ymath
import curve as crv
import ControlPlot as cpr
import zf_gam as zf
import numpy as np
import pyqtgraph.examples
import pylab


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(cpr)
    mix.reload_module(zf)


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


def test_vorbar_1d():
    L = 0.25
    k0 = 2 * np.pi / L
    print(k0)
    print(k0 / (2 * np.pi))

    x = np.linspace(0, 0.4, 801)
    yy = np.zeros([4, np.size(x)])
    for ix in range(np.size(x)):
        x1 = x[ix]
        # yy[0, ix] = np.cos(k0 * x1)
        # yy[1, ix] = np.cos(k0 * x1)
        # yy[2, ix] = np.cos(k0 * x1)
        # yy[3, ix] = np.cos(k0 * x1)
        yy[0, ix] = x1**2 * np.cos(k0 * x1) + x1 * np.sin(3*k0*x1)
        yy[1, ix] = x1**2 * np.cos(k0 * x1) + x1 * np.sin(3*k0*x1)
        yy[2, ix] = x1**2 * np.cos(k0 * x1) + x1 * np.sin(3*k0*x1)
        yy[3, ix] = x1**2 * np.cos(k0 * x1) + x1 * np.sin(3*k0*x1)

    dd_test = {'phibar': {
        't': [0, 1, 2, 3],
        's': x,
        'data': yy}}

    zf.erbar(dd_test)
    zf.vorbar(dd_test)

    # # check
    # x1 = 0.49
    # id_x1, x1_check = mix.find(x, x1)
    # erbar_num = dd_test['erbar']['data'][0, id_x1]
    # erbar_check_x1 = k0 * np.sin(k0 * x1_check)
    #
    # print('x1 = {:0.3f}'.format(x1))
    # print('x1_check = {:0.3f}'.format(x1_check))
    #
    # print('---')
    # print(erbar_check_x1)
    # print(erbar_num)
    # print('---')

    # initial signal
    curves = crv.Curves().xlab('x').ylab('phibar')
    curves.new('yy').XS(x).YS(yy[0, :])
    cpr.plot_curves(curves)

    # erbar
    # erbar_check = k0 * np.sin(k0 * x)
    erbar_check = - 2*x * np.cos(k0 * x) + k0*x**2*np.sin(k0*x)\
                  - np.sin(3*k0*x) - 3*k0*x*np.cos(3*k0*x)

    curves = crv.Curves().xlab('x').ylab('erbar')
    curves.new('yy').XS(x).YS(dd_test['erbar']['data'][0, :]).leg('num')
    curves.new('yy_check').XS(x).YS(erbar_check) \
        .sty(':').col('red').leg('check')
    cpr.plot_curves(curves)

    # vorbar
    # vorbar_check = k0 * np.sin(k0*x) / x + k0**2*np.cos(k0*x)
    vorbar_check = - 4 * np.cos(k0*x) + 2*k0*x*np.sin(k0*x)\
        + 3*k0*x*np.sin(k0*x) + k0**2 * x**2*np.cos(k0*x)\
        - np.sin(3*k0*x)/x - 3*k0*np.cos(3*k0*x)\
        - 6*k0*np.cos(3*k0*x) + 9*k0**2*x*np.sin(3*k0*x)

    curves = crv.Curves().xlab('x').ylab('vorbar')
    curves.new('yy').XS(x).YS(dd_test['vorbar']['data'][0, :]).leg('num')
    curves.new('yy_check').XS(x).YS(vorbar_check)\
        .sty(':').col('red').leg('check')
    cpr.plot_curves(curves)


def test_vorbar():
    T = 1
    w0 = 2 * np.pi / T
    print(w0)
    print(w0 / (2 * np.pi))

    L = 0.25
    k0 = 2 * np.pi / L
    print(k0)
    print(k0 / (2 * np.pi))

    t = np.linspace(0, 1, 201)
    x = np.linspace(0, 1, 801)
    yy = np.zeros([np.size(t), np.size(x)])
    for ix in range(np.size(x)):
        for it in range(np.size(t)):
            x1 = x[ix]
            t1 = t[it]
            yy[it, ix] = np.cos(w0 * t1) * np.cos(k0 * x1)

    dd_test = {'phibar': {
        't': t,
        's': x,
        'data': yy}}

    zf.erbar(dd_test)
    zf.vorbar(dd_test)

    # check
    x1 = 0.49
    t1 = 0.83
    id_t1, t1_ref = mix.find(t, t1)
    id_x1, x1_ref = mix.find(x, x1)
    erbar_num = dd_test['erbar']['data'][id_t1, id_x1]
    erbar_check_t1_x1 = k0 * np.cos(w0 * t1_ref) * np.sin(k0 * x1_ref)

    print('---')
    print(erbar_check_t1_x1)
    print(erbar_num)
    print('---')

    # initial signal
    curves = crv.Curves().xlab('t').ylab('x').tit('phibar')
    curves.new('yy').XS(t).YS(x).ZS(yy)
    cpr.plot_curves_3d(curves)

    # erbar
    yy_check = np.zeros([np.size(t), np.size(x)])
    for ix in range(np.size(x)):
        for it in range(np.size(t)):
            x1 = x[ix]
            t1 = t[it]
            yy_check[it, ix] = k0 * np.cos(w0 * t1) * np.sin(k0 * x1)

    curves = crv.Curves().xlab('t').ylab('x').tit('erbar')
    curves.new('yy').XS(t).YS(x).ZS(dd_test['erbar']['data'])\
        .lev(30)
    cpr.plot_curves_3d(curves)

    curves = crv.Curves().xlab('t').ylab('x').tit('erbar-check')
    curves.new('yy').XS(t).YS(x).ZS(yy_check)\
        .lev(30)
    cpr.plot_curves_3d(curves)

    # vorbar
    yy_check = np.zeros([np.size(t), np.size(x)])
    for ix in range(np.size(x)):
        for it in range(np.size(t)):
            x1 = x[ix]
            t1 = t[it]
            yy_check[it, ix] = (k0 * np.sin(k0*x1) / x1
                + k0**2*np.cos(k0*x1)) * np.cos(w0 * t1)

    curves = crv.Curves().xlab('t').ylab('x').tit('vorbar')
    curves.new('yy').XS(t).YS(x).ZS(dd_test['vorbar']['data']).lev(30)
    cpr.plot_curves_3d(curves)

    curves = crv.Curves().xlab('t').ylab('x').tit('vorbar-check')
    curves.new('yy').XS(t).YS(x).ZS(yy_check).lev(30)
    cpr.plot_curves_3d(curves)


def smooth_demo2():
    x = np.linspace(-4,6,100)
    y = np.sin(x)
    yn = y + pylab.randn(len(y))*0.1

    yn_filt = ymath.smooth(yn, 15)

    curves = crv.Curves().xlab('x').ylab('y')
    curves.new('y').XS(x).YS(y).leg('init')
    curves.new('yn').XS(x).YS(yn).leg('noisy').sty('o')
    curves.new('yn_filt').XS(x).YS(yn_filt).leg('filt').sty(':')
    cpr.plot_curves(curves)


def test_pyqtgraph():
    pyqtgraph.examples.run()

