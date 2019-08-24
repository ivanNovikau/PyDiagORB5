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
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
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
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
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
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
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
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
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


def es_egam_Tf_scan(oo):
    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # *** ES EGAMb: f = 0.01, rho_f = 0.25, v_f = 6.0, scan on Tf ***
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': []
    }

    # --- Er(s1): w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.4,
          [2.431e-03], [2.676e-06],
          [1.843e-04], [2.755e-06]
    )
    add_f(0.6,
          [2.431e-03], [2.676e-06],
          [1.351e-04], [2.423e-06]
        )
    add_f(0.8,
          [2.434e-03], [3.979e-06],
          [9.558e-05], [3.320e-06]
          )
    add_f(1.0,
          [2.441e-03], [4.283e-06],
          [6.304e-05], [2.952e-06]
          )
    data_plot['leg'].append('\overline{E}_r(s = 0.55)')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- MPR: TOT: g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.4, [np.nan], [np.nan],
          [1.819e-04], [4.647e-06]
    )
    add_f(0.6, [np.nan], [np.nan],
          [1.313e-04], [4.938e-06]
          )
    add_f(0.8, [np.nan], [np.nan],
          [9.128e-05], [5.794e-06]
          )
    add_f(1.0, [np.nan], [np.nan],
          [5.624e-05], [7.724e-06]
          )

    data_plot['leg'].append('MPR:\ Total')
    reorganise_data(fs, [], [], gs, gs_err, sel_norm, dd, data_plot)

    # --- MPR: DEUTERIUM: g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.4, [np.nan], [np.nan],
          [-6.449e-04], [7.634e-06]
          )
    add_f(0.6, [np.nan], [np.nan],
          [-5.860e-04], [7.488e-06]
          )
    add_f(0.8, [np.nan], [np.nan],
          [-5.291e-04], [8.612e-06]
          )
    add_f(1.0, [np.nan], [np.nan],
          [-4.779e-04], [8.381e-06]
          )

    gs = np.abs(gs)

    data_plot['leg'].append('MPR:\ Deuterium: |\gamma|')
    reorganise_data(fs, [], [], gs, gs_err, sel_norm, dd, data_plot)

    # --- MPR: FAST: g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.4, [np.nan], [np.nan],
          [8.259e-04], [1.014e-05]
          )
    add_f(0.6, [np.nan], [np.nan],
          [7.186e-04], [1.322e-05]
          )
    add_f(0.8, [np.nan], [np.nan],
          [6.196e-04], [1.449e-05]
          )
    add_f(1.0, [np.nan], [np.nan],
          [5.337e-04], [1.311e-05]
          )

    data_plot['leg'].append('MPR:\ Fast')
    reorganise_data(fs, [], [], gs, gs_err, sel_norm, dd, data_plot)

    # # # *** combined plot ***
    # plot_several_scans(data_plot, sel_norm, xlim=[0.38, 1.02],
    #                    xlab='T_{f}', ylab_g='dynamic\ rates')

    # *** Tf scan for different velocities ***
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': []
    }

    # --- v = 8.0 : w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.4,
          [2.956e-03], [1.829e-06],
          [2.328e-04], [2.355e-06]
          )
    add_f(0.6,
          [2.921e-03], [6.949e-07],
          [2.042e-04], [3.196e-06]
          )
    add_f(0.8,
          [2.894e-03], [1.324e-06],
          [1.782e-04], [2.877e-06]
          )
    add_f(1.0,
          [2.874e-03], [2.032e-06],
          [1.525e-04], [4.064e-06]
          )
    data_plot['leg'].append('\overline{E}_r(s = 0.50): v = 8.0, f = 0.010')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- v = 6.0 : w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.4,
          [2.431e-03], [2.676e-06],
          [1.843e-04], [2.755e-06]
          )
    add_f(0.6,
          [2.431e-03], [2.676e-06],
          [1.351e-04], [2.423e-06]
          )
    add_f(0.8,
          [2.434e-03], [3.979e-06],
          [9.558e-05], [3.320e-06]
          )
    add_f(1.0,
          [2.441e-03], [4.283e-06],
          [6.304e-05], [2.952e-06]
          )
    data_plot['leg'].append('\overline{E}_r(s = 0.55): v = 6.0, f = 0.010')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- v = 3.5 : w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.15,
          [1.276e-03], [1.076e-06],
          [1.451e-04], [2.436e-06]
          )
    add_f(0.20,
          [1.271e-03], [1.074e-06],
          [1.173e-04], [2.003e-06]
          )
    add_f(0.22,
          [1.269e-03], [1.119e-06],
          [1.065e-04], [2.010e-06]
          )
    add_f(0.25,
          [1.268e-03], [1.107e-06],
          [9.063e-05], [2.111e-06]
          )

    data_plot['leg'].append('\overline{E}_r(s = 0.53): v = 3.5, f = 0.095')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- v = 3.0 : w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.10,
          [1.100e-03], [1.200e-06],
          [1.092e-04], [4.093e-06]
          )
    add_f(0.15,
          [1.091e-03], [1.420e-06],
          [7.304e-05], [2.555e-06]
          )
    add_f(0.20,
          [1.086e-03], [4.661e-06],
          [3.819e-05], [6.519e-06]
          )

    data_plot['leg'].append('\overline{E}_r(s = 0.53): v = 3.0, f = 0.095')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # *** combined plot ***
    plot_several_scans(data_plot, sel_norm, xlim=[0.08, 1.02],
                       xlab='T_{f}', ylab_g='dynamic\ rates')


def egam_kin_mie(oo):
    # EM EGAMb: f = 0.01, rho_f = 0.25, v_f = 8.0, scan on ion/mass ratio

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
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(500,
          [2.884e-03], [7.434e-06],
          [3.854e-05], [5.484e-06]
          )
    add_f(1000,
          [2.878e-03], [1.974e-05],
          [4.350e-05], [4.736e-06]
          )
    add_f(2000,
          [2.870e-03], [3.221e-05],
          [5.176e-05], [5.327e-06]
          )
    add_f(3676,
          [2.870e-03], [3.612e-06],
          [5.853e-05], [5.223e-06]
    )

    data_plot['leg'].append('\overline{E}_r(s = 0.50)')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- MPR: TOT: g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(500, [np.nan], [np.nan],
          [], []
    )
    add_f(1000, [np.nan], [np.nan],
          [], []
          )
    add_f(2000, [np.nan], [np.nan],
          [], []
          )
    add_f(3676, [np.nan], [np.nan],
          [], []
          )

    data_plot['leg'].append('MPR:\ Total')
    reorganise_data(fs, [], [], gs, gs_err, sel_norm, dd, data_plot)

    # # --- MPR: DEUTERIUM: g are normalized to wci ---
    # fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    # add_f(500, [np.nan], [np.nan],
    #       [], []
    #       )
    # add_f(1000, [np.nan], [np.nan],
    #       [], []
    #       )
    # add_f(2000, [np.nan], [np.nan],
    #       [], []
    #       )
    # add_f(3676, [np.nan], [np.nan],
    #       [], []
    #       )
    #
    # gs = np.abs(gs)
    #
    # data_plot['leg'].append('MPR:\ Deuterium: |\gamma|')
    # reorganise_data(fs, [], [], gs, gs_err, sel_norm, dd, data_plot)

    # # --- MPR: FAST: g are normalized to wci ---
    # fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    # add_f(500, [np.nan], [np.nan],
    #       [], []
    #       )
    # add_f(1000, [np.nan], [np.nan],
    #       [], []
    #       )
    # add_f(2000, [np.nan], [np.nan],
    #       [], []
    #       )
    # add_f(3676, [np.nan], [np.nan],
    #       [], []
    #       )
    #
    # data_plot['leg'].append('MPR:\ Fast')
    # reorganise_data(fs, [], [], gs, gs_err, sel_norm, dd, data_plot)

    # # *** combined plot ***
    plot_several_scans(data_plot, sel_norm, xlim=[490, 3680],
                       xlab='T_{f}', ylab_g='dynamic\ rates')


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
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
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
    add_f(0.03,
          [2.666e-03], [3.814e-06],
          [2.873e-04], [7.442e-06]
          )
    add_f(0.04,
          [2.597e-03], [3.295e-06],
          [3.121e-04], [8.077e-06]
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

    data_plot['leg'].append('ES\ EGAMb:\ s_{EP} = 0.50')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)
    plot_scan_n(data_plot['f'][-1],
                data_plot['w'][-1], data_plot['w_err'][-1],
                data_plot['g'][-1], data_plot['g_err'][-1],
                sel_norm, xlab='n_{EP}/n_e')

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
    gs_plot     = gs_plot * coef_norm_g
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
            .tit('EGAM\ frequency').xsty('plain')
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
            .tit('EGAM\ growth\ rate').xsty('plain')
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
        line_w = '\omega,\ kHz'
        line_g = '\gamma,\ 1e3/s'
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
    marker_styles = ['o', 's', "x", 'v', "*", "<"]

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


def es_egam_satur_scan(oo):
    # NL ES EGAMb:  scan on EGAM growth rates
    dd = oo.get('dd', None)
    res_data = []
    sel_norm_x = oo.get('sel_norm_x', 'none')

    # x, y normalization
    line_x, coef_x  = '', None
    if sel_norm_x == 'inv-s':
        line_x = '(10^3/s)'
        coef_x = dd['wc'] / 1.e3

    def create_init_dict():
        data_save = {'xs': [], 'ys': [], 'xs_err': [], 'ys_err': [], 'leg': '', 'sty': 'o', 'col': 'blue'}
        return data_save

    def add_f(x, y):
        data_plot['xs'].append(x[0])
        data_plot['xs_err'].append(x[1])
        data_plot['ys'].append(y[0])
        data_plot['ys_err'].append(y[1])

    def plot_several_data(res_data_plot, xlab='x', ylab='y', tit=''):
        curves = crv.Curves() \
            .xlab(xlab + line_x) \
            .ylab(ylab) \
            .tit(tit).xsty('plain')
        for id_data in range(len(res_data_plot)):
            data_plot_one = res_data_plot[id_data]
            curves.new() \
                .XS(np.array(data_plot_one['xs']) * coef_x) \
                .YS(np.array(data_plot_one['ys'])) \
                .set_errorbar(True,
                              xs=data_plot_one['xs_err'],
                              ys=data_plot_one['ys_err']
                              ) \
                .leg(data_plot_one['leg']) \
                .sty(data_plot_one['sty']) \
                .col(data_plot_one['col'])
        cpr.plot_curves(curves)

        return

    # scan on fpart
    data_plot = create_init_dict()

    add_f(
        [2.210e-05, 5.217e-06],
        [2.472e+00, 3.226e-02])  # f = 0.0040, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [1.520e-04, 4.687e-06],
        [9.628e+00, 9.795e-01])  # f = 0.0100, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [2.455e-04, 6.985e-06],
        [3.360e+01, 4.055e+00])  # f = 0.0200, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [2.873e-04, 7.442e-06],
        [3.639e+01, 1.798e+01])  # f = 0.0300, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [3.121e-04, 8.077e-06],
        [3.521e+01, 1.517e+01])  # f = 0.0400, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [3.251e-04, 7.383e-06],
        [3.465e+01, 1.775e+01])  # f = 0.0500, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [3.477e-04, 4.980e-05],
        [2.989e+01, 1.233e+01])  # f = 0.0949, vp = 8.0, T = 1.0, sf = 0.25

    data_plot['leg'] = ['n_{EP}/n_e = [0.004, 0.01, 0.02, 0.03, 0.04, 0.05, 0.095],\ ',
                       'v_{\parallel, EP} = 8.0,\ T = 1.0']
    data_plot['sty'] = 'o'
    data_plot['col'] = 'blue'
    res_data.append(data_plot)

    # scan on vp
    data_plot = create_init_dict()

    add_f(
        [1.052e-04, 4.224e-06],
        [7.079e+00, 2.106e-01]
    )  # f = 0.0100, vp = 9.0, T = 1.0, sf = 0.25
    add_f(
        [1.316e-04, 2.818e-06],
        [5.861e+00, 1.991e+00]
    )  # f = 0.0100, vp = 7.0, T = 1.0, sf = 0.25
    add_f(
        [1.483e-04, 2.606e-06],
        [7.744e+00, 1.266e+00]
    )  # f = 0.0100, vp = 7.5, T = 1.0, sf = 0.25

    data_plot['leg'] = 'v_{\parallel, EP} = [9.0, 7.0, 7.5],\ n_{EP}/n_e = 0.01,\ T_{EP} = 1.0'
    data_plot['sty'] = 's'
    data_plot['col'] = 'red'
    res_data.append(data_plot)

    # different vp and fpart
    data_plot = create_init_dict()

    add_f(
        [1.173e-04, 2.003e-06],
        [8.121e+00, 1.322e+00]
    )  # f = 0.0945, vp = 3.5, T = 0.2, sf = 0.25

    data_plot['leg'] = 'v_{\parallel, EP} = 3.5,\ n_{EP}/n_e = 0.095,\ T_{EP} = 0.20'
    data_plot['sty'] = '^'
    data_plot['col'] = 'green'
    res_data.append(data_plot)

    # scan on T
    data_plot = create_init_dict()

    add_f(
        [1.843e-04, 2.755e-06],
        [1.310e+01, 2.307e+00]
    )  # T = 0.4, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.351e-04, 2.423e-06],
        [9.375e+00, 1.820e+00]
    )  # T = 0.6, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [9.558e-05, 3.320e-06],
        [4.269e+00, 5.464e-01]
    )  # T = 0.8, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [6.304e-05, 2.952e-06],
        [1.441e+00, 4.716e-02]
    )  # T = 1.0, f = 0.01, vp = 6.0, sf = 0.25

    data_plot['leg'] = ['T_{EP} = [0.4, 0.6, 0.8, 1.0],\ ',
                       'v_{\parallel, EP} = 6.0,\ n_{EP}/n_e = 0.01']
    data_plot['sty'] = '*'
    data_plot['col'] = 'black'
    res_data.append(data_plot)

    # # plot results
    # plot_several_data(res_data, xlab='\gamma', ylab='\sqrt{<\overline{E}^2>_{s,t}}',
    #           tit='t = [{:0.3e}, {:0.3e}]'.format(5e4, 1.37e5))
    plot_several_data(res_data, xlab='\gamma', ylab='\sqrt{<\overline{E}^2>_{s,t}}',
                      tit='s_{EP} = 0.50:\ Saturation\ levels')


