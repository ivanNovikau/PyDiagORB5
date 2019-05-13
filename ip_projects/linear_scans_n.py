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
    ## --- ATTENTION ---
    # Be carefull when you chose normalization, since
    # you use only one project (one dictionary dd) to define normalization,
    # different scans can have different normalization, e.g., for the case
    # of different scans with different temperature profiles

    dd = oo.get('dd_aug', None)
    sel_norm = oo.get('sel_norm', 'wc')

    ns, ws, gs = [], [], []
    data_plot = {'n': [], 'w': [], 'g': [], 'leg': []}

    # ws[i] and gs[i] should have the same size
    def add_n(n, w, g):
        ns.append(n)
        ws.append(w)
        gs.append(g)

    # HERE, w, g are normalized to wci
    add_n(10, [2.092e-04], [1.097e-04])
    add_n(20, [3.530e-04], [3.057e-04])
    add_n(30, [6.193e-04], [4.681e-04])
    add_n(40, [9.129e-04], [5.417e-04])
    add_n(50, [1.046e-04, 1.255e-03],
              [5.386e-04, 5.386e-04])
    add_n(60, [1.793e-04], [7.117e-04])
    add_n(70, [2.718e-04], [8.464e-04])
    add_n(80, [4.488e-04], [9.597e-04])
    add_n(90, [6.910e-04], [1.026e-03])
    add_n(100, [8.958e-05], [1.101e-03])
    add_n(110, [2.090e-04], [1.250e-03])
    add_n(120, [4.295e-04], [1.327e-03])
    add_n(130, [6.352e-04], [1.369e-03])
    add_n(140, [8.415e-04], [1.415e-03])
    add_n(150, [1.047e-03], [1.442e-03])
    add_n(160, [1.263e-03], [1.445e-03])
    add_n(170, [1.274e-03], [1.470e-03])
    add_n(180, [1.714e-03], [1.407e-03])
    add_n(190, [1.858e-03], [1.374e-03])
    add_n(200, [2.049e-03], [1.369e-03])

    data_plot['leg'].append('initial')
    plot_scan_n(ns, ws, gs, sel_norm, dd, data_plot)


def tcv_scan_n(oo={}):
    dd = oo.get('dd_tcv', None)
    sel_norm = oo.get('sel_norm', 'wc')

    ns, ws, gs = [], [], []
    data_plot = {'n': [], 'w': [], 'g': [], 'leg': []}

    ## --- ATTENTION ---
    # Be carefull when you chose normalization, since
    # you use only one project (one dictionary dd) to define normalization,
    # different scans can have different normalization, e.g., for the case
    # of different scans with different temperature profiles

    # ws[i] and gs[i] should have the same size
    def add_n(n, w, g):
        ns.append(n)
        ws.append(w)
        gs.append(g)

    # --- mu_T = 12, mu_n = 5: w, g are normalized to wci ---
    ns, ws, gs = [], [], []
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

    data_plot['leg'].append('\mu_T = 12,\ \mu_n = 5')
    reorganise_data(ns, ws, gs, sel_norm, dd, data_plot)
    plot_scan_n(data_plot['n'][-1],
                data_plot['w'][-1],
                data_plot['g'][-1], sel_norm)

    # --- mu_T = 6, mu_n = 5: w, g are normalized to wci ---
    ns, ws, gs = [], [], []
    add_n(10, [3.640e-03], [9.856e-04])
    add_n(20, [7.065e-04], [1.135e-03])
    add_n(30, [1.963e-03, 1.806e-03, 1.492e-03],
              [1.060e-03, 1.060e-03, 1.060e-03])
    add_n(40, [1.021e-03], [1.164e-03])
    add_n(50, [1.479e-03], [1.212e-03])
    add_n(60, [1.973e-03], [1.100e-03])
    add_n(75, [2.601e-03], [7.232e-04])
    add_n(80, [7.850e-05, 2.748e-03],
              [5.764e-04, 5.764e-04])
    add_n(95, [6.531e-04], [6.054e-04])

    data_plot['leg'].append('\mu_T = 6,\ \mu_n = 5')
    reorganise_data(ns, ws, gs, sel_norm, dd, data_plot)

    # --- mu_T = 12, mu_n = 3: w, g are normalized to wci ---
    ns, ws, gs = [], [], []
    add_n(40, [4.780e-04], [1.538e-03])
    add_n(50, [8.113e-04], [1.871e-03])
    add_n(60, [1.235e-03], [2.032e-03])
    add_n(70, [1.730e-03], [2.227e-03])
    add_n(75, [1.962e-03], [2.211e-03])
    add_n(80, [2.180e-03], [2.131e-03])
    add_n(90, [2.643e-03], [1.969e-03])
    add_n(95, [2.933e-03], [1.664e-03])
    add_n(100, [3.191e-03], [1.319e-03])
    add_n(110, [3.453e-03], [5.247e-04])

    data_plot['leg'].append('\mu_T = 12,\ \mu_n = 3')
    reorganise_data(ns, ws, gs, sel_norm, dd, data_plot)

    # --- mu_T = 6, mu_n = 3: w, g are normalized to wci ---
    ns, ws, gs = [], [], []
    add_n(50, [1.418e-03], [1.240e-03])
    add_n(60, [7.200e-04], [1.251e-03])
    add_n(75, [1.237e-03], [1.367e-03])
    add_n(80, [1.405e-03], [1.285e-03])
    add_n(95, [1.918e-03], [8.790e-04])

    data_plot['leg'].append('\mu_T = 6,\ \mu_n = 3')
    reorganise_data(ns, ws, gs, sel_norm, dd, data_plot)

    # --- mu_T = 9, mu_n = 5: w, g are normalized to wci ---
    ns, ws, gs = [], [], []
    add_n(50, [1.435e-03, 1.705e-03],
              [1.143e-03, 1.143e-03])
    add_n(60, [1.974e-03, 2.691e-03, 2.422e-03],
              [1.030e-03, 1.030e-03, 1.030e-03])
    add_n(75, [3.452e-04], [1.383e-03])
    add_n(80, [5.359e-04], [1.405e-03])
    add_n(95, [1.241e-03], [1.186e-03])

    data_plot['leg'].append('\mu_T = 9,\ \mu_n = 5')
    reorganise_data(ns, ws, gs, sel_norm, dd, data_plot)

    # *** combined plot ***
    plot_several_scans(data_plot, sel_norm)


def reorganise_data(ns, ws, gs, sel_norm, dd, data_plot):
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
    if sel_norm == 'khz':
        coef_norm_w = dd['wc'] / (2 * np.pi * 1.e3)
        coef_norm_g = dd['wc'] / 1.e3
    if sel_norm == 'wc':
        coef_norm_w = coef_norm_g = 1
    if sel_norm == 'csa':
        coef_norm_w = coef_norm_g = dd['wc'] / (dd['cs'] / dd['a0'])
    if sel_norm == 'csr':
        coef_norm_w = coef_norm_g = dd['wc'] / (dd['cs'] / dd['R0'])
    ws_plot = ws_plot * coef_norm_w
    gs_plot = gs_plot * coef_norm_g

    data_plot['n'].append(ns_plot)
    data_plot['w'].append(ws_plot)
    data_plot['g'].append(gs_plot)


def plot_scan_n(ns_plot, ws_plot, gs_plot, sel_norm):
    # normalization:
    line_w, line_g = '', ''
    if sel_norm == 'khz':
        line_w = '\omega,\ kHz',
        line_g = '\gamma,\ 1e3/s',
    if sel_norm == 'wc':
        line_w = '\omega[\omega_c]'
        line_g = '\gamma[\omega_c]'
    if sel_norm == 'csa':
        line_w = '\omega[c_s/a_0]'
        line_g = '\gamma[c_s/a_0]'
    if sel_norm == 'csr':
        line_w = '\omega[c_s/R_0]'
        line_g = '\gamma[c_s/R_0]'

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
    curves_g.flag_legend = False
    curves_g.new('g').XS(ns_plot).YS(gs_plot).sty('o')
    cpr.plot_curves(curves_g)

    return


def plot_several_scans(data_plot, sel_norm):
    # normalization:
    line_w, line_g = '', ''
    if sel_norm == 'khz':
        line_w = '\omega,\ kHz',
        line_g = '\gamma,\ 1e3/s',
    if sel_norm == 'wc':
        line_w = '\omega[\omega_c]'
        line_g = '\gamma[\omega_c]'
    if sel_norm == 'csa':
        line_w = '\omega[c_s/a_0]'
        line_g = '\gamma[c_s/a_0]'
    if sel_norm == 'csr':
        line_w = '\omega[c_s/R_0]'
        line_g = '\gamma[c_s/R_0]'

    # plotting
    number_scans = np.shape(data_plot['n'])[0]
    marker_styles = ['o', 's', "x", '*', "v", "<"]

    curves_w = crv.Curves() \
        .xlab('n') \
        .ylab(line_w) \
        .tit('frequency').xsty('plain')\
        .xlim([0, 150]) # !!!
    curves_g = crv.Curves() \
        .xlab('n') \
        .ylab(line_g) \
        .tit('growth\ rate').xsty('plain')\
        .xlim([0, 150]) # !!!
    for i_scan in range(number_scans):
        n = data_plot['n'][i_scan]
        w = data_plot['w'][i_scan]
        g = data_plot['g'][i_scan]
        leg = data_plot['leg'][i_scan]
        msty = marker_styles[i_scan]
        curves_w.new('w').XS(n).YS(w).sty(msty).leg(leg)
        curves_g.new('g').XS(n).YS(g).sty(msty).leg(leg)
    cpr.plot_curves(curves_w)
    cpr.plot_curves(curves_g)

    return


def tcv_egam(oo):
    dd = oo.get('dd_tcv', None)
    sel_norm = oo.get('sel_norm', 'wc')

    sp, ws, gs, ls = [], [], [], []
    data_plot = {'n': [], 'w': [], 'g': [], 'leg': []}

    # --- ATTENTION : be carefull when you chose normalization

    # len(n) = len(ls), len(w[is]) = len(g[is])
    def add_n(n, w, g, leg):
        ns.append(n)
        ws.append(w)
        gs.append(g)
        ls.append(leg)

    # --- mu_T = 12, mu_n = 5: w, g are normalized to wci ---
    ns, ws, gs, ls = [], [], [], []
    add_n(40, [2.997e-03], [1.091e-03])

    # reorganise_data(ns, ws, gs, sel_norm, dd, data_plot)
    # plot_scan_n(data_plot['n'][-1],
    #             data_plot['w'][-1],
    #             data_plot['g'][-1], sel_norm)

    # *** combined plot ***
    # plot_several_scans(data_plot, sel_norm)




