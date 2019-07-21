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
    dd = oo.get('dd_tcv', None)
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
    # plot_scan_n(data_plot['f'][-1],
    #             data_plot['w'][-1], data_plot['w_err'][-1],
    #             data_plot['g'][-1], data_plot['g_err'][-1],
    #             sel_norm)

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

    # *** combined plot ***
    plot_several_scans(data_plot, sel_norm, xlim=[0.0, 0.1])


def reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot):
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

    data_plot['f'].append(fs_plot)
    data_plot['w'].append(ws_plot)
    data_plot['w_err'].append(ws_err_plot)
    data_plot['g'].append(gs_plot)
    data_plot['g_err'].append(gs_err_plot)


def plot_scan_n(fs_plot, ws_plot, ws_err_plot, gs_plot, gs_err_plot, sel_norm):
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
    curves_w.new('w')\
        .XS(fs_plot)\
        .YS(ws_plot)\
        .set_errorbar(True, ys=ws_err_plot)\
        .sty('o')
    cpr.plot_curves(curves_w)

    curves_g = crv.Curves()\
        .xlab('n')\
        .ylab(line_g)\
        .tit('growth\ rate').xsty('plain')
    curves_g.new('g')\
        .XS(fs_plot)\
        .YS(gs_plot) \
        .set_errorbar(True, ys=gs_err_plot) \
        .sty('o')
    cpr.plot_curves(curves_g)

    return


def plot_several_scans(data_plot, sel_norm, xlim=[0, 1]):
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
        .xlab('f_{part}') \
        .ylab(line_w) \
        .tit('frequency').xsty('plain')\
        .xlim(xlim)
    curves_g = crv.Curves() \
        .xlab('f_{part}') \
        .ylab(line_g) \
        .tit('growth\ rate').xsty('plain')\
        .xlim(xlim)
    for i_scan in range(number_scans):
        f = data_plot['f'][i_scan]
        w = data_plot['w'][i_scan]
        w_err = data_plot['w_err'][i_scan]
        g = data_plot['g'][i_scan]
        g_err = data_plot['g_err'][i_scan]
        leg = data_plot['leg'][i_scan]
        msty = marker_styles[i_scan]
        curves_w.new('w')\
            .XS(f).YS(w).set_errorbar(True, ys=w_err)\
            .sty(msty).leg(leg)
        curves_g.new('g')\
            .XS(f).YS(g).set_errorbar(True, ys=g_err)\
            .sty(msty).leg(leg)
    cpr.plot_curves(curves_w)
    cpr.plot_curves(curves_g)

    return