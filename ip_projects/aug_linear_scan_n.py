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
    add_n(30, [6.201e-04], [4.684e-04])
    add_n(60, [1.480e-03, 1.461e-04], [4.498e-04, 7.173e-04])

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



