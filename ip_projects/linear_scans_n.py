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


def aug_scan_n():
    # w, g normalized to wci

    ns, ws, gs = [], [], []

    # ws[i] and gs[i] should have the same size
    def add_n(n, w, g):
        ns.append(n)
        ws.append(w)
        gs.append(g)

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

    plot_scan_n(ns, ws, gs)


def tcv_scan_n():
    # w, g normalized to wci
    ns, ws, gs = [], [], []

    # ws[i] and gs[i] should have the same size
    def add_n(n, w, g):
        ns.append(n)
        ws.append(w)
        gs.append(g)

    # n10 - no frequency, just growth
    # n128 - no growth

    add_n(50,  [1.173e-03], [1.074e-03])
    add_n(60,  [9.848e-05], [1.451e-03])
    add_n(70,  [4.039e-04], [1.756e-03])
    add_n(80,  [8.544e-04], [1.818e-03])
    add_n(85,  [1.013e-03], [1.906e-03])
    add_n(90,  [1.258e-03], [1.769e-03])
    add_n(95,  [1.628e-03], [1.542e-03])
    add_n(100, [1.958e-04], [1.268e-03])
    add_n(105, [2.244e-03], [9.492e-04])
    add_n(110, [2.398e-03], [5.436e-04])

    plot_scan_n(ns, ws, gs)


def plot_scan_n(ns, ws, gs):
    ns_plot, ws_plot, gs_plot = [], [], []
    for id_n in range(len(ns)):
        for id_wg in range(len(ws[id_n])):
            ns_plot.append(ns[id_n])
            ws_plot.append(ws[id_n][id_wg])
            gs_plot.append(gs[id_n][id_wg])

    curves_w = crv.Curves().xlab('n').ylab('\omega[\omega_{ci}]')
    curves_w.new('w').XS(ns_plot).YS(ws_plot).sty('o')
    cpr.plot_curves(curves_w)

    curves_g = crv.Curves().xlab('n').ylab('\gamma[\omega_{ci}]')
    curves_g.new('g').XS(ns_plot).YS(gs_plot).sty('o')
    cpr.plot_curves(curves_g)



