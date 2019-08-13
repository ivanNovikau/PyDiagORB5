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


def es_egam_fpart_scan(oo):
    # ES EGAMb:  v|| = 8.0: scan on fpart

    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')

    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': []
    }

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # --- B: s = 0.50: w, g are normalized to wci ---
    fs, ws, gs = [], [], []
    add_f(0.01,
          [2.877e-03], [7.720e-05],
          [1.548e-04], [9.219e-06]
        )
    add_f(0.02,
          [2.757e-03], [7.000e-05],
          [2.455e-04], [7.076e-06]
          )
    add_f(0.05,
          [2.538e-03], [6.496e-05],
          [3.245e-04], [5.588e-06]
          )
    add_f(0.07,
          [2.444e-03], [5.002e-05],
          [3.216e-04], [1.827e-05]
          )
    add_f(0.09,
          [2.400e-03], [1.171e-04],
          [3.250e-04], [3.813e-05]
          )

    data_plot['leg'].append('ES\ EGAMb:\ s = 0.50')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)
    plot_scan_n(data_plot['f'][-1],
                data_plot['w'][-1], data_plot['w_err'][-1],
                data_plot['g'][-1], data_plot['g_err'][-1],
                sel_norm)

    # --- B: s = 0.60: w, g are normalized to wci ---
    fs, ws, gs = [], [], []
    add_f(0.01,
          [2.876e-03], [7.636e-05],
          [1.521e-04], [4.387e-05]
          )
    add_f(0.02,
          [2.756e-03], [6.211e-05],
          [2.461e-04], [3.318e-05]
          )
    add_f(0.05,
          [2.545e-03], [6.293e-05],
          [3.293e-04], [2.115e-05]
          )
    add_f(0.07,
          [2.461e-03], [5.497e-05],
          [3.388e-04], [1.889e-05]
          )
    add_f(0.09,
          [2.368e-03], [9.121e-05],
          [3.440e-04], [7.004e-05]
          )
    data_plot['leg'].append('ES\ EGAMb:\ s = 0.60')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- M: s = 0.50: w, g are normalized to wci ---
    fs, ws, gs = [], [], []
    add_f(0.01,
          [2.814e-03], [9.739e-05],
          [9.944e-05], [5.117e-06]
          )
    add_f(0.02,
          [2.690e-03], [8.337e-05],
          [1.541e-04], [5.834e-06]
          )
    add_f(0.05,
          [2.473e-03], [6.150e-05],
          [1.884e-04], [1.077e-05]
          )
    add_f(0.07,
          [2.400e-03], [1.159e-04],
          [1.934e-04], [3.701e-05]
          )
    add_f(0.09,
          [2.341e-03], [9.016e-05],
          [1.945e-04], [2.617e-05]
          )

    data_plot['leg'].append('ES\ EGAMm:\ s = 0.50')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- M: s = 0.60: w, g are normalized to wci ---
    fs, ws, gs = [], [], []
    add_f(0.01,
          [2.816e-03], [6.068e-05],
          [9.734e-05], [3.544e-05]
          )
    add_f(0.02,
          [2.686e-03], [6.128e-05],
          [1.578e-04], [1.798e-05]
          )
    add_f(0.05,
          [2.495e-03], [8.152e-05],
          [1.936e-04], [1.206e-05]
          )
    add_f(0.07,
          [2.444e-03], [4.644e-05],
          [1.828e-04], [2.972e-05]
          )
    add_f(0.09,
          [2.402e-03], [6.177e-05],
          [1.867e-04], [9.143e-06]
          )

    data_plot['leg'].append('ES\ EGAMm:\ s = 0.60')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # # *** combined plot ***
    # plot_several_scans(data_plot, sel_norm, xlim=[0.0, 0.1])


def es_egam_vp_scan(oo):
    # ES EGAMb: f = 0.01, rho_f = 0.25, scan on v_parallel

    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')

    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': []
    }

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # --- Er(s1): w, g are normalized to wci ---
    fs, ws, gs = [], [], []
    add_f(7.0,
          [2.652e-03], [1.884e-05],
          [1.316e-04], [2.818e-06]
    )
    add_f(7.5,
          [2.763e-03], [3.188e-06],
          [1.483e-04], [2.606e-06]
          )
    add_f(7.8,
          [2.830e-03], [3.203e-06],
          [1.537e-04], [4.429e-06]
          )
    add_f(8.0,
          [2.876e-03], [2.247e-06],
          [1.508e-04], [3.596e-06]
    )
    add_f(8.2,
          [2.918e-03], [1.951e-06],
          [1.487e-04], [3.689e-06]
          )
    add_f(8.5,
          [2.988e-03], [3.114e-05],
          [1.382e-04], [4.353e-06]
          )
    add_f(9.0,
          [3.102e-03], [5.984e-06],
          [1.052e-04], [4.224e-06]
          )

    data_plot['leg'].append('\overline{E}_r(s = 0.5)')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)
    # plot_scan_n(data_plot['f'][-1],
    #             data_plot['w'][-1], data_plot['w_err'][-1],
    #             data_plot['g'][-1], data_plot['g_err'][-1],
    #             sel_norm, xlab='v_{\parallel, f}')

    # --- MPR: TOT: g are normalized to wci ---
    fs, ws, gs = [], [], []
    add_f(7.0, [np.nan], [np.nan],
          [1.263e-04], [4.267e-06]
    )
    add_f(7.5, [np.nan], [np.nan],
          [1.455e-04], [2.866e-06]
          )
    add_f(7.8, [np.nan], [np.nan],
          [1.509e-04], [5.397e-06]
          )
    add_f(8.0, [np.nan], [np.nan],
          [1.487e-04], [2.796e-06]
    )
    add_f(8.2, [np.nan], [np.nan],
          [1.462e-04], [4.858e-06]
          )
    add_f(8.5, [np.nan], [np.nan],
          [1.366e-04], [6.874e-06]
          )
    add_f(9.0, [np.nan], [np.nan],
          [1.025e-04], [5.031e-06]
          )

    data_plot['leg'].append('MPR:\ Total')
    reorganise_data(fs, [], [], gs, gs_err, sel_norm, dd, data_plot)
    # plot_scan_n(data_plot['f'][-1],
    #             data_plot['w'][-1], data_plot['w_err'][-1],
    #             data_plot['g'][-1], data_plot['g_err'][-1],
    #             sel_norm, xlab='v_{\parallel, f}')

    # --- MPR: DEUTERIUM: g are normalized to wci ---
    fs, ws, gs = [], [], []
    add_f(7.0, [np.nan], [np.nan],
          [-5.304e-04], [7.938e-06]
          )
    add_f(7.5, [np.nan], [np.nan],
          [-5.001e-04], [5.246e-06]
          )
    add_f(7.8, [np.nan], [np.nan],
          [-4.719e-04], [4.041e-06]
          )
    add_f(8.0, [np.nan], [np.nan],
          [-4.550e-04], [4.624e-06]
          )
    add_f(8.2, [np.nan], [np.nan],
          [-4.288e-04], [3.504e-06]
          )
    add_f(8.5, [np.nan], [np.nan],
          [-3.923e-04], [5.716e-06]
          )
    add_f(9.0, [np.nan], [np.nan],
          [-3.223e-04], [5.491e-06]
          )

    gs = np.abs(gs)

    data_plot['leg'].append('MPR:\ Deuterium: |\gamma|')
    reorganise_data(fs, [], [], gs, gs_err, sel_norm, dd, data_plot)
    # plot_scan_n(data_plot['f'][-1],
    #             data_plot['w'][-1], data_plot['w_err'][-1],
    #             data_plot['g'][-1], data_plot['g_err'][-1],
    #             sel_norm, xlab='v_{\parallel, f}')

    # --- MPR: FAST: g are normalized to wci ---
    fs, ws, gs = [], [], []
    add_f(7.0, [np.nan], [np.nan],
          [6.569e-04], [4.802e-06]
          )
    add_f(7.5, [np.nan], [np.nan],
          [6.460e-04], [5.626e-06]
          )
    add_f(7.8, [np.nan], [np.nan],
          [6.219e-04], [4.478e-06]
          )
    add_f(8.0, [np.nan], [np.nan],
          [6.038e-04], [1.533e-06]
          )
    add_f(8.2, [np.nan], [np.nan],
          [5.755e-04], [2.573e-06]
          )
    add_f(8.5, [np.nan], [np.nan],
          [5.282e-04], [2.108e-06]
          )
    add_f(9.0, [np.nan], [np.nan],
          [4.250e-04], [1.707e-06]
          )

    data_plot['leg'].append('MPR:\ Fast')
    reorganise_data(fs, [], [], gs, gs_err, sel_norm, dd, data_plot)
    # plot_scan_n(data_plot['f'][-1],
    #             data_plot['w'][-1], data_plot['w_err'][-1],
    #             data_plot['g'][-1], data_plot['g_err'][-1],
    #             sel_norm, xlab='v_{\parallel, f}')

    # # *** combined plot ***
    plot_several_scans(data_plot, sel_norm, xlim=[6.8, 9.2],
                       xlab='v_{\parallel,f}', ylab_g='dynamic\ rates')


def es_egam_fpart_scan_rho025(oo):
    # ES EGAMb: rho_f = 0.25, v|| = 8.0: scan on fpart

    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')

    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': []
    }

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # --- B: s = 0.50: w, g are normalized to wci ---
    fs, ws, gs = [], [], []
    add_f(0.004,
          [2.978e-03], [6.040e-06],
          [2.210e-05], [5.217e-06]
          )
    add_f(0.006,
          [2.939e-03], [2.530e-06],
          [8.123e-05], [5.416e-06]
          )
    add_f(0.008,
          [2.904e-03], [8.635e-06],
          [1.214e-04], [3.571e-06]
          )
    add_f(0.01,
          [2.875e-03], [2.414e-06],
          [1.520e-04], [4.687e-06]
        )
    add_f(0.02,
          [2.754e-03], [3.013e-06],
          [2.452e-04], [6.985e-06]
          )
    add_f(0.05,
          [2.539e-03], [5.065e-06],
          [3.251e-04], [7.383e-06]
          )
    add_f(0.07,
          [2.453e-03], [9.048e-06],
          [3.364e-04], [9.816e-06]
          )
    add_f(0.09,
          [2.375e-03], [4.630e-05],
          [3.477e-04], [4.980e-05]
          )

    data_plot['leg'].append('ES\ EGAMb:\ s = 0.50')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)
    plot_scan_n(data_plot['f'][-1],
                data_plot['w'][-1], data_plot['w_err'][-1],
                data_plot['g'][-1], data_plot['g_err'][-1],
                sel_norm, xlab='f_{EP}')

    # # *** combined plot ***
    # plot_several_scans(data_plot, sel_norm, xlim=[0.0, 0.1])


def reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot):
    flag_no_freq = False
    if len(ws) == 0:
        ws = np.array(gs)
        ws_err = np.array(gs_err)
        flag_no_freq = True

    # reorganise and save data as np.array
    fs_plot, ws_plot, ws_err_plot, gs_plot, gs_err_plot = [], [], [], [], []
    for id_n in range(len(fs)):
        for id_wg in range(len(ws[id_n])):
            fs_plot.append(fs[id_n])
            ws_plot.append(ws[id_n][id_wg])
            ws_err_plot.append(ws_err[id_n][id_wg])
            gs_plot.append(gs[id_n][id_wg])
            gs_err_plot.append(gs_err[id_n][id_wg])
    fs_plot = np.array(fs_plot)
    ws_plot = np.array(ws_plot)
    ws_err_plot = np.array(ws_err_plot)
    gs_plot = np.array(gs_plot)
    gs_err_plot = np.array(gs_err_plot)

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
    ws_err_plot = ws_err_plot * coef_norm_w
    gs_plot = gs_plot * coef_norm_g
    gs_err_plot = gs_err_plot * coef_norm_g

    if flag_no_freq:
        ws_plot = []
        ws_err_plot = []

    data_plot['f'].append(fs_plot)
    data_plot['w'].append(ws_plot)
    data_plot['w_err'].append(ws_err_plot)
    data_plot['g'].append(gs_plot)
    data_plot['g_err'].append(gs_err_plot)


def plot_scan_n(fs_plot, ws_plot, ws_err_plot, gs_plot, gs_err_plot,
                sel_norm, xlab='n'):
    # normalization:
    line_w, line_g = '', ''
    if sel_norm == 'khz':
        line_w = '\omega,\ kHz'
        line_g = '\gamma,\ 10^3/s'
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
    if len(ws_plot) > 0:
        curves_w = crv.Curves()\
            .xlab(xlab)\
            .ylab(line_w)\
            .tit('frequency').xsty('plain')
        curves_w.new('w')\
            .XS(fs_plot)\
            .YS(ws_plot)\
            .set_errorbar(True, ys=ws_err_plot)\
            .sty('o')
        cpr.plot_curves(curves_w)

    if len(gs_plot) > 0:
        curves_g = crv.Curves()\
            .xlab(xlab)\
            .ylab(line_g)\
            .tit('growth\ rate').xsty('plain')
        curves_g.new('g')\
            .XS(fs_plot)\
            .YS(gs_plot) \
            .set_errorbar(True, ys=gs_err_plot) \
            .sty('o')
        cpr.plot_curves(curves_g)

    return


def plot_several_scans(data_plot, sel_norm, xlim=[0, 1], xlab='f_{part}', ylab_g='growth\ rate'):
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
    number_scans = np.shape(data_plot['f'])[0]
    marker_styles = ['o', 's', "x", '*', "v", "<"]

    curves_w = crv.Curves() \
        .xlab(xlab) \
        .ylab(line_w) \
        .tit('frequency').xsty('plain')\
        .xlim(xlim)
    curves_g = crv.Curves() \
        .xlab(xlab) \
        .ylab(line_g) \
        .tit(ylab_g).xsty('plain')\
        .xlim(xlim)
    for i_scan in range(number_scans):
        f = data_plot['f'][i_scan]
        w = data_plot['w'][i_scan]
        w_err = data_plot['w_err'][i_scan]
        g = data_plot['g'][i_scan]
        g_err = data_plot['g_err'][i_scan]
        leg = data_plot['leg'][i_scan]
        msty = marker_styles[i_scan]
        if len(w) > 0:
            curves_w.new('w')\
                .XS(f).YS(w).set_errorbar(True, ys=w_err)\
                .sty(msty).leg(leg)
        curves_g.new('g')\
            .XS(f).YS(g).set_errorbar(True, ys=g_err)\
            .sty(msty).leg(leg)
    cpr.plot_curves(curves_w)
    cpr.plot_curves(curves_g)

    return
