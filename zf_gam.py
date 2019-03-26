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


def fft_gam_s1(dd, s1, oo={}):
    rd.phibar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']

    cs = dd['cs']

    id_s1, _ = mix.find(s, s1)
    phibar_s1 = dd['phibar']['data'][:, id_s1]

    w, w2 = ymath.w_for_fft(t)
    f, f2 = ymath.fft_y(phibar_s1)

    # initial signal
    curves = crv.Curves().xlab('t').ylab('y')
    curves.new('y').XS(t).YS(y).leg('init')
    cpr.plot_curves(curves)

    # fourier spectrum
    curves = crv.Curves().xlab('w').ylab('f')
    curves.new('f').XS(w).YS(f).leg('fft')
    cpr.plot_curves(curves)

    return


def fft_gam_2d(oo={}):
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
