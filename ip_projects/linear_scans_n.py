import Mix as mix
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np
from scipy import interpolate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)


def aug_scan_n(oo={}):
    dd = oo.get('dd_aug', None)
    sel_norm = oo.get('sel_norm', 'wc')

    ns, ws, gs = [], [], []

    # ws[i] and gs[i] should have the same size
    def add_n(n, w, g):
        ns.append(n)
        ws.append(w)
        gs.append(g)

    # HERE, w, g are normalized to wci
    add_n(10, [1.614e-04], [1.111e-04])
    add_n(20, [3.977e-04], [3.050e-04])
    add_n(30, [6.201e-04], [4.684e-04])
    add_n(40, [9.123e-04], [5.417e-04])
    add_n(50, [1.000e-03], [5.348e-04])
    add_n(60, [1.480e-03, 1.461e-04], [4.498e-04, 7.173e-04])
    add_n(70, [2.732e-04], [8.455e-04])
    add_n(80, [4.500e-04], [9.582e-04])
    add_n(90, [6.900e-04], [1.027e-03])
    add_n(100, [1.083e-04], [1.088e-03])
    add_n(110, [2.380e-04], [1.240e-03])
    add_n(120, [4.284e-04], [1.326e-03])
    add_n(130, [6.363e-04], [1.364e-03])
    add_n(140, [8.390e-04], [1.413e-03])
    add_n(150, [1.047e-03], [1.439e-03])
    add_n(160, [1.266e-03], [1.445e-03])
    add_n(170, [1.274e-03], [1.469e-03])
    add_n(180, [1.714e-03], [1.407e-03])
    add_n(190, [1.858e-03], [1.373e-03])
    add_n(200, [2.049e-03], [1.371e-03])

    plot_scan_n(ns, ws, gs, sel_norm, dd)


def tcv_scan_n(oo={}):
    dd = oo.get('dd_tcv', None)
    sel_norm = oo.get('sel_norm', 'wc')

    ns, ws, gs = [], [], []

    # ws[i] and gs[i] should have the same size
    def add_n(n, w, g):
        ns.append(n)
        ws.append(w)
        gs.append(g)

    # n10 - no frequency, just growth
    # n128 - no growth

    # HERE, w, g are normalized to wci
    add_n(40, [2.997e-03], [1.091e-03])
    add_n(45, [1.256e-03], [1.099e-03])
    add_n(50, [1.525e-03, 1.794e-04, 1.705e-03],
              [1.063e-03, 1.063e-03, 1.063e-03])
    add_n(55, [6.344e-05], [1.275e-03])
    add_n(58, [3.4e-05],   [1.377e-03])
    add_n(60, [7.850e-05], [1.459e-03])
    add_n(62, [7.850e-05], [1.581e-03])
    add_n(65, [2.355e-04], [1.666e-03])
    add_n(70, [4.039e-04], [1.756e-03])
    add_n(75, [6.077e-04], [1.826e-03])
    add_n(78, [7.397e-04], [1.824e-03])
    add_n(80, [8.306e-04], [1.823e-03])
    add_n(82, [8.947e-04], [1.906e-03])
    add_n(85, [1.013e-03], [1.906e-03])
    add_n(90, [1.258e-03], [1.769e-03])
    add_n(95, [1.621e-03], [1.544e-03])
    add_n(98, [1.828e-03], [1.421e-03])
    add_n(100, [1.957e-03], [1.267e-03])
    add_n(102, [2.084e-03], [1.178e-03])
    add_n(105, [2.244e-03], [9.492e-04])
    add_n(110, [2.398e-03], [5.436e-04])

    plot_scan_n(ns, ws, gs, sel_norm, dd)


def plot_scan_n(ns, ws, gs, sel_norm, dd):
    # reorganise and save data as np.array
    ns_plot, ws_plot, gs_plot = [], [], []
    for id_n in range(len(ns)):
        for id_wg in range(len(ws[id_n])):
            ns_plot.append(ns[id_n])
            ws_plot.append(ws[id_n][id_wg])
            gs_plot.append(gs[id_n][id_wg])
    ns_plot = np.array(ns_plot)
    ws_plot = np.array(ws_plot)
    gs_plot = np.array(gs_plot)

    # normalization:
    coef_norm_w, coef_norm_g = np.NaN, np.NaN
    line_w, line_g = '', ''
    if sel_norm == 'khz':
        coef_norm_w = dd['wc'] / (2 * np.pi * 1.e3)
        coef_norm_g = dd['wc'] / 1.e3
        line_w = '\omega,\ kHz',
        line_g = '\gamma,\ 1e3/s',
    if sel_norm == 'wc':
        coef_norm_w = coef_norm_g = 1
        line_w = '\omega[\omega_c]'
        line_g = '\gamma[\omega_c]'
    if sel_norm == 'csa':
        coef_norm_w = coef_norm_g = dd['wc'] / (dd['cs'] / dd['a0'])
        line_w = '\omega[c_s/a_0]'
        line_g = '\gamma[c_s/a_0]'
    if sel_norm == 'csr':
        coef_norm_w = coef_norm_g = dd['wc'] / (dd['cs'] / dd['R0'])
        line_w = '\omega[c_s/R_0]'
        line_g = '\gamma[c_s/R_0]'
    ws_plot = ws_plot * coef_norm_w
    gs_plot = gs_plot * coef_norm_g

    # plotting
    curves_w = crv.Curves()\
        .xlab('n')\
        .ylab(line_w)\
        .tit('frequency').xsty('plain')
    curves_w.flag_legend = False
    curves_w.new('w').XS(ns_plot).YS(ws_plot).sty('o')
    cpr.plot_curves(curves_w)

    curves_g = crv.Curves()\
        .xlab('n')\
        .ylab(line_g)\
        .tit('growth\ rate').xsty('plain')
    curves_g.flag_legend = True
    curves_g.new('g').XS(ns_plot).YS(gs_plot).sty('o')
    cpr.plot_curves(curves_g)



