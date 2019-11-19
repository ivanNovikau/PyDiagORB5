import Mix as mix
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np
from scipy import interpolate
from scipy import constants


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
    # data_plot['leg'].append('\overline{E}_r(s = 0.50): v_{\parallel, EP} = 8.0, n_{EP}/n_e = 0.010')
    data_plot['leg'].append('v_{\parallel, EP} = 8.0, n_{EP}/n_e = 0.010')
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
    # data_plot['leg'].append('\overline{E}_r(s = 0.55): v_{\parallel, EP} = 6.0, n_{EP}/n_e = 0.010')
    data_plot['leg'].append('v_{\parallel, EP} = 6.0, n_{EP}/n_e = 0.010')
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

    # data_plot['leg'].append('\overline{E}_r(s = 0.53): v_{\parallel, EP} = 3.5, n_{EP}/n_e = 0.095')
    data_plot['leg'].append('v_{\parallel, EP} = 3.5, n_{EP}/n_e = 0.095')
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

    # data_plot['leg'].append('\overline{E}_r(s = 0.53): v_{\parallel, EP} = 3.0, n_{EP}/n_e = 0.095')
    data_plot['leg'].append('v_{\parallel, EP} = 3.0, n_{EP}/n_e = 0.095')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # *** combined plot ***
    plot_several_scans(data_plot, sel_norm, xlim=[0.08, 1.02],
                       xlab='T_{EP}', ylab_g='EGAM\ growth\ rate', tit_w='EGAM\ frequency')


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


def egam_fpart_AE_KE_wg(oo):
    # EGAMb: adiabatic vs drift-kinetic electrons: rho_f = 0.25, v|| = 8.0, T = 1.0:
    # SCAN ON n_{EP}/n_{e}

    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')

    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': [], 'sty': [], 'col': []
    }

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # --- ADIAB. ELE: DIRECT: w, g are normalized to wci ---
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

    data_plot['leg'].append('AE:\ from\ \overline{E}_r(s\ of\ EGAM\ localisation)')
    data_plot['sty'].append('o')
    data_plot['col'].append('blue')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # ***
    # ******
    # --- DK. ELE: DIRECT: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [2.868e-03], [5.402e-06],
          [5.525e-05], [4.799e-06])
    add_f(0.0200,
          [2.767e-03], [7.330e-06],
          [1.124e-04], [5.080e-06])
    add_f(0.0300,
          [2.689e-03], [9.091e-06],
          [1.323e-04], [9.729e-06])
    add_f(0.0500,
          [2.602e-03], [1.673e-05],
          [1.266e-04], [3.506e-05])
    add_f(0.0700,
          [2.582e-03], [1.785e-05],
          [1.545e-04], [5.912e-05])
    add_f(0.0900,
          [2.593e-03], [2.102e-05],
          [1.341e-04], [1.286e-05])

    data_plot['leg'].append('KE:\ from\ \overline{E}_r(s\ of\ EGAM\ localisation)')
    data_plot['sty'].append('s')
    data_plot['col'].append('red')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # *** combined plot ***
    plot_several_scans_adv(data_plot, sel_norm, xlim=[0.0, 0.1],
                           ylab_g='EGAM\ growth\ rate', tit_w='EGAM\ frequency')


def egam_fpart_AE_KE_wg_ele(oo):
    # EGAMb: adiabatic vs drift-kinetic electrons: rho_f = 0.25, v|| = 8.0, T = 1.0:
    # SCAN ON n_{EP}/n_{e}

    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')

    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': [], 'sty': [], 'col': []
    }

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # --- ADIAB. ELE: DIRECT: w, g are normalized to wci ---
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

    data_plot['leg'].append('AE:\ from\ \overline{E}_r(s\ of\ EGAM\ localisation)')
    data_plot['sty'].append('o')
    data_plot['col'].append('blue')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # ***
    # ******
    # --- DK. ELE: DIRECT: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [2.868e-03], [5.402e-06],
          [5.525e-05], [4.799e-06])
    add_f(0.0200,
          [2.767e-03], [7.330e-06],
          [1.124e-04], [5.080e-06])
    add_f(0.0300,
          [2.689e-03], [9.091e-06],
          [1.323e-04], [9.729e-06])
    add_f(0.0500,
          [2.602e-03], [1.673e-05],
          [1.266e-04], [3.506e-05])
    add_f(0.0700,
          [2.582e-03], [1.785e-05],
          [1.545e-04], [5.912e-05])
    add_f(0.0900,
          [2.593e-03], [2.102e-05],
          [1.341e-04], [1.286e-05])

    data_plot['leg'].append('KE:\ from\ \overline{E}_r(s\ of\ EGAM\ localisation)')
    data_plot['sty'].append('s')
    data_plot['col'].append('red')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- DK. ELE: MPR THERMAL ELECTRONS: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [np.nan], [np.nan],
          [-6.439e-05], [7.856e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [-6.974e-05], [4.502e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [-6.475e-05], [2.567e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [-4.938e-05], [4.053e-06])
    add_f(0.0700,
          [np.nan], [np.nan],
          [-5.120e-05], [2.323e-06])
    add_f(0.0900,
          [np.nan], [np.nan],
          [-5.094e-05], [1.876e-06])

    data_plot['leg'].append('KE:\ MPR:\ thermal\ electrons: |\gamma|')
    data_plot['sty'].append('s')
    data_plot['col'].append('orange')
    reorganise_data(fs, ws, ws_err, np.abs(gs), gs_err, sel_norm, dd, data_plot)

    # *** combined plot ***
    plot_several_scans_adv(data_plot, sel_norm, xlim=[0.0, 0.1],
                           ylab_g='EGAM\ growth\ rate', tit_w='EGAM\ frequency')


def egam_fpart_AE_KE_g_mpr(oo):
    # EGAMb: adiabatic vs drift-kinetic electrons: rho_f = 0.25, v|| = 8.0, T = 1.0:
    # SCAN ON n_{EP}/n_{e}

    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')

    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': [], 'sty': [], 'col': []
    }

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # --- ADIAB. ELE.: MPR TOTAL: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0060,
          [np.nan], [np.nan],
          [7.735e-05], [6.706e-06])
    add_f(0.0100,
          [np.nan], [np.nan],
          [1.503e-04], [4.945e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [2.412e-04], [3.550e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [2.845e-04], [2.703e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [3.165e-04], [2.586e-06])
    add_f(0.0700,
          [np.nan], [np.nan],
          [3.135e-04], [7.597e-06])
    add_f(0.0900,
          [np.nan], [np.nan],
          [3.116e-04], [2.222e-05])

    # data_plot['leg'].append('AE:\ total')
    data_plot['leg'].append('AE: \gamma_{TOT}')
    data_plot['sty'].append('o')
    data_plot['col'].append('red')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- ADIAB. ELE.: MPR THERMAL DEUTERIUM: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0060,
          [np.nan], [np.nan],
          [-3.906e-04], [6.946e-06])
    add_f(0.0100,
          [np.nan], [np.nan],
          [-4.512e-04], [2.937e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [-5.355e-04], [4.854e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [-5.788e-04], [4.660e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [-6.139e-04], [1.244e-05])
    add_f(0.0700,
          [np.nan], [np.nan],
          [-6.053e-04], [2.832e-05])
    add_f(0.0900,
          [np.nan], [np.nan],
          [-6.285e-04], [5.390e-05])

    # data_plot['leg'].append('AE:\ therm.\ deut.: |\gamma|')
    data_plot['leg'].append('AE: |\gamma_D|')
    data_plot['sty'].append('o')
    data_plot['col'].append('green')
    reorganise_data(fs, ws, ws_err, np.abs(gs), gs_err, sel_norm, dd, data_plot)

    # --- ADIAB. ELE.: MPR EP: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0060,
          [np.nan], [np.nan],
          [4.684e-04], [1.275e-06])
    add_f(0.0100,
          [np.nan], [np.nan],
          [6.012e-04], [4.036e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [7.770e-04], [4.770e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [8.637e-04], [5.102e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [9.295e-04], [1.132e-05])
    add_f(0.0700,
          [np.nan], [np.nan],
          [9.195e-04], [3.621e-05])
    add_f(0.0900,
          [np.nan], [np.nan],
          [9.489e-04], [6.361e-05])

    # data_plot['leg'].append('AE:\ EP')
    data_plot['leg'].append('AE: \gamma_{EP}')
    data_plot['sty'].append('o')
    data_plot['col'].append('black')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # ***
    # ******
    # --- DK. ELE: MPR TOTAL: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [np.nan],    [np.nan],
          [5.621e-05], [4.360e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [1.287e-04], [2.505e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [1.549e-04], [4.018e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [1.568e-04], [1.437e-05])
    add_f(0.0700,
          [np.nan], [np.nan],
          [1.470e-04], [1.528e-05])
    add_f(0.0900,
          [np.nan], [np.nan],
          [1.582e-04], [1.171e-05])

    # data_plot['leg'].append('KE:\ total')
    data_plot['leg'].append('KE: \gamma_{TOT}')
    data_plot['sty'].append('*')
    data_plot['col'].append('red')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)
    #
    # --- DK. ELE: MPR THERMAL DEUTERIUM: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [np.nan], [np.nan],
          [-5.415e-04], [4.849e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [-6.531e-04], [1.794e-05])
    add_f(0.0300,
          [np.nan], [np.nan],
          [-7.107e-04], [3.389e-05])
    add_f(0.0500,
          [np.nan], [np.nan],
          [-7.622e-04], [7.616e-05])
    add_f(0.0700,
          [np.nan], [np.nan],
          [-8.681e-04], [6.217e-05])
    add_f(0.0900,
          [np.nan], [np.nan],
          [-9.322e-04], [7.633e-05])

    # data_plot['leg'].append('KE:\ therm.\ deut.: |\gamma|')
    data_plot['leg'].append('KE: |\gamma_D|')
    data_plot['sty'].append('*')
    data_plot['col'].append('green')
    reorganise_data(fs, ws, ws_err, np.abs(gs), gs_err, sel_norm, dd, data_plot)
    #
    # --- DK. ELE: MPR THERMAL ELECTRONS: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [np.nan], [np.nan],
          [-6.439e-05], [7.856e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [-6.974e-05], [4.502e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [-6.475e-05], [2.567e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [-4.938e-05], [4.053e-06])
    add_f(0.0700,
          [np.nan], [np.nan],
          [-5.120e-05], [2.323e-06])
    add_f(0.0900,
          [np.nan], [np.nan],
          [-5.094e-05], [1.876e-06])

    # data_plot['leg'].append('KE:\ therm.\ elec.: |\gamma|')
    data_plot['leg'].append('KE: |\gamma_E|')
    data_plot['sty'].append('*')
    data_plot['col'].append('orange')
    reorganise_data(fs, ws, ws_err, np.abs(gs), gs_err, sel_norm, dd, data_plot)
    #
    # --- DK. ELE: MPR EP: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [np.nan], [np.nan],
          [6.609e-04], [5.312e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [8.514e-04], [1.699e-05])
    add_f(0.0300,
          [np.nan], [np.nan],
          [9.311e-04], [3.851e-05])
    add_f(0.0500,
          [np.nan], [np.nan],
          [9.716e-04], [1.057e-04])
    add_f(0.0700,
          [np.nan], [np.nan],
          [1.069e-03], [9.174e-05])
    add_f(0.0900,
          [np.nan], [np.nan],
          [1.148e-03], [9.418e-05])

    # data_plot['leg'].append('KE:\ EP')
    data_plot['leg'].append('KE: \gamma_{EP}')
    data_plot['sty'].append('*')
    data_plot['col'].append('black')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # *** combined plot ***
    plot_several_scans_adv(data_plot, sel_norm, xlim=[0.0, 0.1],
                           ylab_g='EGAM\ growth\ rate\ from\ MPR', tit_w='EGAM\ frequency')


def egam_fpart_AE_KE_g_mpr_tot(oo):
    # EGAMb: adiabatic vs drift-kinetic electrons: rho_f = 0.25, v|| = 8.0, T = 1.0:
    # SCAN ON n_{EP}/n_{e}

    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')

    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': [], 'sty': [], 'col': []
    }

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # --- ADIAB. ELE.: MPR TOTAL: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0060,
          [np.nan], [np.nan],
          [7.735e-05], [6.706e-06])
    add_f(0.0100,
          [np.nan], [np.nan],
          [1.503e-04], [4.945e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [2.412e-04], [3.550e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [2.845e-04], [2.703e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [3.165e-04], [2.586e-06])
    add_f(0.0700,
          [np.nan], [np.nan],
          [3.135e-04], [7.597e-06])
    add_f(0.0900,
          [np.nan], [np.nan],
          [3.116e-04], [2.222e-05])

    # data_plot['leg'].append('AE:\ total')
    data_plot['leg'].append('AE: \gamma')
    data_plot['sty'].append('o')
    data_plot['col'].append('blue')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # ***
    # ******
    # --- DK. ELE: MPR TOTAL: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [np.nan],    [np.nan],
          [5.621e-05], [4.360e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [1.287e-04], [2.505e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [1.549e-04], [4.018e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [1.568e-04], [1.437e-05])
    add_f(0.0700,
          [np.nan], [np.nan],
          [1.470e-04], [1.528e-05])
    add_f(0.0900,
          [np.nan], [np.nan],
          [1.582e-04], [1.171e-05])

    # data_plot['leg'].append('KE:\ total')
    data_plot['leg'].append('KE: \gamma')
    data_plot['sty'].append('*')
    data_plot['col'].append('red')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # *** combined plot ***
    plot_several_scans_adv(data_plot, sel_norm, xlim=[0.0, 0.1],
                           ylab_g='EGAM\ growth\ rate\ from\ MPR', tit_w='EGAM\ frequency')


def egam_fpart_AE_KE_g_mpr_ele(oo):
    # EGAMb: adiabatic vs drift-kinetic electrons: rho_f = 0.25, v|| = 8.0, T = 1.0:
    # SCAN ON n_{EP}/n_{e}

    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')

    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': [], 'sty': [], 'col': []
    }

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # --- ADIAB. ELE.: MPR TOTAL: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0060,
          [np.nan], [np.nan],
          [7.735e-05], [6.706e-06])
    add_f(0.0100,
          [np.nan], [np.nan],
          [1.503e-04], [4.945e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [2.412e-04], [3.550e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [2.845e-04], [2.703e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [3.165e-04], [2.586e-06])
    add_f(0.0700,
          [np.nan], [np.nan],
          [3.135e-04], [7.597e-06])
    add_f(0.0900,
          [np.nan], [np.nan],
          [3.116e-04], [2.222e-05])

    data_plot['leg'].append('AE:\ total')
    data_plot['sty'].append('o')
    data_plot['col'].append('blue')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # ***
    # ******
    # --- DK. ELE: MPR TOTAL: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [np.nan],    [np.nan],
          [5.621e-05], [4.360e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [1.287e-04], [2.505e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [1.549e-04], [4.018e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [1.568e-04], [1.437e-05])
    add_f(0.0700,
          [np.nan], [np.nan],
          [1.470e-04], [1.528e-05])
    add_f(0.0900,
          [np.nan], [np.nan],
          [1.582e-04], [1.171e-05])

    data_plot['leg'].append('KE:\ total')
    data_plot['sty'].append('s')
    data_plot['col'].append('red')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)
    #
    # --- DK. ELE: MPR THERMAL ELECTRONS: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [np.nan], [np.nan],
          [-6.439e-05], [7.856e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [-6.974e-05], [4.502e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [-6.475e-05], [2.567e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [-4.938e-05], [4.053e-06])
    add_f(0.0700,
          [np.nan], [np.nan],
          [-5.120e-05], [2.323e-06])
    add_f(0.0900,
          [np.nan], [np.nan],
          [-5.094e-05], [1.876e-06])

    data_plot['leg'].append('KE:\ therm.\ elec.: |\gamma|')
    data_plot['sty'].append('s')
    data_plot['col'].append('orange')
    reorganise_data(fs, ws, ws_err, np.abs(gs), gs_err, sel_norm, dd, data_plot)

    # *** combined plot ***
    plot_several_scans_adv(data_plot, sel_norm, xlim=[0.0, 0.1],
                           ylab_g='EGAM\ growth\ rate\ from\ MPR', tit_w='EGAM\ frequency')


def scan_sigma_rho025(oo):
    # EGAMb: adiabatic vs drift-kinetic electrons: rho_f = 0.25, v|| = 8.0, T = 1.0:
    # SCAN ON n_{EP}/n_{e}

    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')

    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': [], 'sty': [], 'col': []
    }

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # --- ADIAB. ELE: DIRECT: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.05,
          [2.772e-03], [1.256e-06],
          [2.406e-04], [2.732e-06]
          )
    add_f(0.10,
          [2.874e-03], [1.905e-06],
          [1.530e-04], [3.896e-06]
          )
    add_f(0.15,
          [2.951e-03], [6.691e-06],
          [9.768e-05], [1.012e-05]
          )
    add_f(0.20,
          [3.016e-03], [4.911e-06],
          [8.796e-05], [1.138e-05]
          )

    data_plot['leg'].append('n_{EP}/n_e = 0.01')
    data_plot['sty'].append('o')
    data_plot['col'].append('blue')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- ADIAB. ELE: DIRECT: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.05,
          [2.291e-03], [1.275e-05],
          [3.112e-04], [6.794e-06])
    add_f(0.10,
          [2.362e-03], [1.049e-05],
          [3.669e-04], [7.301e-06])
    add_f(0.15,
          [2.513e-03], [1.992e-06],
          [2.893e-04], [1.119e-05])
    add_f(0.20,
          [2.590e-03], [1.083e-06],
          [3.370e-04], [8.098e-06])

    data_plot['leg'].append('n_{EP}/n_e = 0.095')
    data_plot['sty'].append('o')
    data_plot['col'].append('red')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # *** combined plot ***
    plot_several_scans_adv(data_plot, sel_norm, xlim=[0.04, 0.22], xlab='\sigma_{EP}',
                           ylab_g='EGAM\ growth\ rate', tit_w='EGAM\ frequency')


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


def plot_several_scans(data_plot, sel_norm, xlim=[0, 1], xlab='n_{EP}/n_e', ylab_g='growth\ rate', tit_w='frequency'):
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
        .tit(tit_w).xsty('plain')\
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


def plot_several_scans_adv(data_plot, sel_norm, oo_texts=[], xlim=[0, 1], xlab='n_{EP}/n_e',
                           ylab_g='EGAM\ growth\ rate', tit_w='frequency'):
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
    color_freq = ['blue', 'red']
    number_scans = np.shape(data_plot['f'])[0]

    curves_w = crv.Curves() \
        .xlab(xlab) \
        .ylab(line_w) \
        .tit(tit_w).xsty('plain')\
        .xlim(xlim)
    curves_g = crv.Curves() \
        .xlab(xlab) \
        .ylab(line_g) \
        .tit(ylab_g).xsty('plain')\
        .xlim(xlim)
    for oo_text in oo_texts:
        oText = crv.PlText(oo_text)
        curves_g.newt(oText)

    count_w = -1
    for i_scan in range(number_scans):
        f = data_plot['f'][i_scan]
        w = data_plot['w'][i_scan]
        w_err = data_plot['w_err'][i_scan]
        g = data_plot['g'][i_scan]
        g_err = data_plot['g_err'][i_scan]
        leg = data_plot['leg'][i_scan]
        msty = data_plot['sty'][i_scan]
        col = data_plot['col'][i_scan]
        if len(w) > 0:
            curves_w.new('w')\
                .XS(f).YS(w).set_errorbar(True, ys=w_err)\
                .sty(msty)
            if not np.isnan(w[0]):
                count_w += 1
                curves_w.list_curves[-1].leg(leg).col(color_freq[count_w])
        curves_g.new('g')\
            .XS(f).YS(g).set_errorbar(True, ys=g_err)\
            .sty(msty).col(col).leg(leg)
    cpr.plot_curves(curves_w)
    cpr.plot_curves(curves_g)

    return


def es_egam_satur_scan(oo):
    # NL ES EGAMb:  scan on EGAM growth rates
    dd = oo.get('dd', None)
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

    def plot_several_data(res_data_plot, xlab='x', ylab='y', tit='', ylim=None):
        curves = crv.Curves() \
            .xlab(xlab + line_x) \
            .ylab(ylab) \
            .tit(tit).xsty('plain').ylim(ylim)
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

    # ------------------------------------------------
    # --- SATURATION LEVELS OF OSCILLATING SIGNALS ---
    # ------------------------------------------------
    res_data = []

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
    # add_f(
    #     [2.873e-04, 7.442e-06],
    #     [3.639e+01, 1.798e+01])  # f = 0.0300, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [3.121e-04, 8.077e-06],
    #     [3.521e+01, 1.517e+01])  # f = 0.0400, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [3.251e-04, 7.383e-06],
    #     [3.465e+01, 1.775e+01])  # f = 0.0500, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [3.477e-04, 4.980e-05],
    #     [2.989e+01, 1.233e+01])  # f = 0.0949, vp = 8.0, T = 1.0, sf = 0.25

    # data_plot['leg'] = ['n_{EP}/n_e = [0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 9.5]\cdot 10^{-2},\ ',
    #                     'v_{\parallel, EP} = 8.0,\ T_{EP} = 1.0', '---------------']
    data_plot['leg'] = ['n_{EP}/n_e = [0.4, 1.0, 2.0]\cdot 10^{-2},\ v_{\parallel, EP} = 8.0,\ T_{EP} = 1.0']
    data_plot['sty'] = 'o'
    data_plot['col'] = 'blue'
    res_data.append(data_plot)

    # scan on vp
    data_plot = create_init_dict()

    add_f(
        [1.052e-04, 4.224e-06],
        [7.079e+00, 2.106e-01])  # f = 0.0100, vp = 9.0, T = 1.0, sf = 0.25
    add_f(
        [1.316e-04, 2.818e-06],
        [5.861e+00, 1.991e+00])  # f = 0.0100, vp = 7.0, T = 1.0, sf = 0.25
    add_f(
        [1.483e-04, 2.606e-06],
        [7.744e+00, 1.266e+00])  # f = 0.0100, vp = 7.5, T = 1.0, sf = 0.25

    # data_plot['leg'] = ['v_{\parallel, EP} = [9.0, 7.0, 7.5],',
    #                     'n_{EP}/n_e = 0.01,\ T_{EP} = 1.0', '---------------']
    data_plot['leg'] = ['v_{\parallel, EP} = [9.0, 7.0, 7.5],\ n_{EP}/n_e = 0.01,\ T_{EP} = 1.0']
    data_plot['sty'] = 's'
    data_plot['col'] = 'red'
    res_data.append(data_plot)

    # different vp and fpart
    data_plot = create_init_dict()

    add_f(
        [1.451e-04, 2.436e-06],
        [1.234e+01, 4.366e+00])  # f = 0.0945, vp = 3.5, T = 0.15, sf = 0.25
    add_f(
        [1.173e-04, 2.003e-06],
        [9.332e+00, 5.342e-01])  # f = 0.0945, vp = 3.5, T = 0.20, sf = 0.25
    add_f(
        [1.065e-04, 2.010e-06],
        [9.374e+00, 5.542e-01])  # f = 0.0945, vp = 3.5, T = 0.22, sf = 0.25
    add_f(
        [9.063e-05, 2.111e-06],
        [7.429e+00, 2.944e-01])  # f = 0.0945, vp = 3.5, T = 0.25, sf = 0.25

    # data_plot['leg'] = ['T_{EP} = [0.25, 0.22, 0.20, 0.15],',
    #                     'v_{\parallel, EP} = 3.5,\ n_{EP}/n_e = 0.095', '---------------']
    data_plot['leg'] = ['T_{EP} = [0.25, 0.22, 0.20, 0.15],\ v_{\parallel, EP} = 3.5,\ n_{EP}/n_e = 0.095']
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

    # data_plot['leg'] = ['T_{EP} = [1.0, 0.8, 0.6, 0.4],\ ',
    #                    'v_{\parallel, EP} = 6.0,\ n_{EP}/n_e = 0.010']
    data_plot['leg'] = ['T_{EP} = [1.0, 0.8, 0.6, 0.4],\ v_{\parallel, EP} = 6.0,\ n_{EP}/n_e = 0.010']
    data_plot['sty'] = '*'
    data_plot['col'] = 'black'
    res_data.append(data_plot)

    # plot results: saturation levels of oscillating signals
    plot_several_data(res_data,
                      xlab='\gamma_{egam, lin}',
                      ylab='\sqrt{<\overline{E}^2>_{s,t}}\ \ (a.u.)',
                      tit='s_{EP} = 0.50:\ Saturation\ levels')

    # ---------------------------------------------
    # --- SATURATION LEVELS OF SMOOTHED SIGNALS ---
    # ---------------------------------------------
    res_data = []

    # scan on fpart: vp = 8.0, T = 1.0, sf = 0.25
    data_plot = create_init_dict()

    add_f(
        [2.210e-05, 5.217e-06],
        [1.901e+00, 1.532e-02])  # f = 0.0040
    add_f(
        [1.520e-04, 4.687e-06],
        [6.612e+00, 3.132e-01])  # f = 0.0100
    add_f(
        [2.455e-04, 6.985e-06],
        [2.648e+01, 1.069e+01])  # f = 0.0200
    add_f(
        [2.873e-04, 7.442e-06],
        [2.453e+01, 5.137e+00])  # f = 0.0300
    add_f(
        [3.121e-04, 8.077e-06],
        [2.555e+01, 6.931e+00])  # f = 0.0400
    add_f(
        [3.251e-04, 7.383e-06],
        [2.362e+01, 7.496e-01])  # f = 0.0500
    add_f(
        [3.477e-04, 4.980e-05],
        [2.590e+01, 3.141e+00])  # f = 0.0949

    data_plot['leg'] = ['n_{EP}/n_e = ', '[0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 9.5]\cdot 10^{-2},\ ',
                        'v_{\parallel, EP} = 8.0,\ T_{EP} = 1.0',
                        '------------']
    data_plot['sty'] = 'o'
    data_plot['col'] = 'blue'
    res_data.append(data_plot)

    # scan on vp: f = 0.0100, T = 1.0, sf = 0.25
    data_plot = create_init_dict()

    add_f(
        [1.052e-04, 4.224e-06],
        [6.338e+00, 4.949e-02])  # vp = 9.0
    add_f(
        [1.316e-04, 2.818e-06],
        [3.253e+00, 4.578e-01])  # vp = 7.0
    add_f(
        [1.483e-04, 2.606e-06],
        [4.780e+00, 2.052e-01])  # vp = 7.5

    data_plot['leg'] = ['v_{\parallel, EP} = [9.0, 7.0, 7.5],',
                        'n_{EP}/n_e = 0.01,\ T_{EP} = 1.0',
                        '------------']
    data_plot['sty'] = 's'
    data_plot['col'] = 'red'
    res_data.append(data_plot)

    # different T: sf = 0.25, f = 0.0945, vp = 3.5
    data_plot = create_init_dict()

    add_f(
        [1.451e-04, 2.436e-06],
        [9.044e+00, 2.656e+00])  # T = 0.15
    add_f(
        [1.173e-04, 2.003e-06],
        [7.152e+00, 7.933e-02])  # T = 0.20
    add_f(
        [1.065e-04, 2.010e-06],
        [6.768e+00, 1.566e-01])  # T = 0.22
    add_f(
        [9.063e-05, 2.111e-06],
        [5.583e+00, 8.869e-02])  # T = 0.25

    data_plot['leg'] = ['T_{EP} = [0.25, 0.22, 0.20, 0.15],',
                        'v_{\parallel, EP} = 3.5,\ n_{EP}/n_e = 0.095',
                        '------------']
    data_plot['sty'] = '^'
    data_plot['col'] = 'green'
    res_data.append(data_plot)

    # scan on T: f = 0.01, vp = 6.0, sf = 0.25
    data_plot = create_init_dict()

    add_f(
        [6.304e-05, 2.952e-06],
        [1.129e+00, 5.064e-02])  # T = 1.0
    add_f(
        [9.558e-05, 3.320e-06],
        [3.360e+00, 2.376e-01])  # T = 0.8
    add_f(
        [1.351e-04, 2.423e-06],
        [6.739e+00, 5.434e-01])  # T = 0.6
    add_f(
        [1.843e-04, 2.755e-06],
        [9.912e+00, 1.118e+00])  # T = 0.4

    data_plot['leg'] = ['T_{EP} = [1.0, 0.8, 0.6, 0.4],\ ',
                        'v_{\parallel, EP} = 6.0,\ n_{EP}/n_e = 0.010']
    data_plot['sty'] = '*'
    data_plot['col'] = 'black'
    res_data.append(data_plot)

    # plot results: saturation levels of smoothed signals
    plot_several_data(res_data, xlab='\gamma_{egam, lin}', ylab='\sqrt{<\overline{E}^2>_{s,t}}',
                      tit='s_{EP} = 0.50:\ Saturation\ levels')


def es_egam_satur_scan_overshoot(oo):
    # NL ES EGAMb:  scan on EGAM growth rates
    dd = oo.get('dd', None)
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

    def plot_several_data(res_data_plot, xlab='x', ylab='y', tit='', ylim=None):
        curves = crv.Curves() \
            .xlab(xlab + line_x) \
            .ylab(ylab) \
            .tit(tit).xsty('plain').ylim(ylim)
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

    # ------------------------------------------------
    # --- SATURATION LEVELS OF OSCILLATING SIGNALS ---
    # ------------------------------------------------
    res_data = []

    # scan on fpart
    data_plot = create_init_dict()

    add_f(
        [2.210e-05, 5.217e-06],
        [2.270e+00, 0])  # f = 0.0040, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [1.520e-04, 4.687e-06],
        [1.689e+01, 0])  # f = 0.0100, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [2.455e-04, 6.985e-06],
        [3.427e+01, 0])  # f = 0.0200, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [2.873e-04, 7.442e-06],
    #     [4.085e+01, 0])  # f = 0.0300, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [3.121e-04, 8.077e-06],
    #     [4.664e+01, 0])  # f = 0.0400, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [3.251e-04, 7.383e-06],
    #     [5.176e+01, 0])  # f = 0.0500, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [3.477e-04, 4.980e-05],
    #     [6.665e+01, 0])  # f = 0.0949, vp = 8.0, T = 1.0, sf = 0.25

    # data_plot['leg'] = \
    #     ['n_{EP}/n_e = [0.4, 1, 2, 3, 4, 5, 9]\cdot 10^{-2},\ v_{\parallel, EP} = 8.0,\ T_{EP} = 1.0']
    data_plot['leg'] = \
        ['n_{EP}/n_e = [0.4, 1.0, 2.0]\cdot 10^{-2},\ v_{\parallel, EP} = 8.0,\ T_{EP} = 1.0']
    data_plot['sty'] = 'o'
    data_plot['col'] = 'blue'
    res_data.append(data_plot)

    # scan on vp
    data_plot = create_init_dict()

    add_f(
        [1.052e-04, 4.224e-06],
        [1.103e+01, 0])  # f = 0.0100, vp = 9.0, T = 1.0, sf = 0.25
    add_f(
        [1.316e-04, 2.818e-06],
        [1.246e+01, 0])  # f = 0.0100, vp = 7.0, T = 1.0, sf = 0.25
    add_f(
        [1.483e-04, 2.606e-06],
        [1.593e+01, 0])  # f = 0.0100, vp = 7.5, T = 1.0, sf = 0.25

    data_plot['leg'] = ['v_{\parallel, EP} = [9.0, 7.0, 7.5],\ n_{EP}/n_e = 0.01,\ T_{EP} = 1.0']
    data_plot['sty'] = 's'
    data_plot['col'] = 'red'
    res_data.append(data_plot)

    # different vp and fpart
    data_plot = create_init_dict()

    add_f(
        [9.063e-05, 2.111e-06],
        [7.809e+00, 0])  # f = 0.0945, vp = 3.5, T = 0.25, sf = 0.25
    add_f(
        [1.065e-04, 2.010e-06],
        [9.929e+00, 0])  # f = 0.0945, vp = 3.5, T = 0.22, sf = 0.25
    add_f(
        [1.173e-04, 2.003e-06],
        [1.121e+01, 0])  # f = 0.0945, vp = 3.5, T = 0.20, sf = 0.25
    add_f(
        [1.451e-04, 2.436e-06],
        [1.636e+01, 0])  # f = 0.0945, vp = 3.5, T = 0.15, sf = 0.25

    data_plot['leg'] = ['T_{EP} = [0.25, 0.22, 0.20, 0.15],\ v_{\parallel, EP} = 3.5,\ n_{EP}/n_e = 0.095']
    data_plot['sty'] = '^'
    data_plot['col'] = 'green'
    res_data.append(data_plot)

    # scan on T
    data_plot = create_init_dict()

    add_f(
        [6.304e-05, 2.952e-06],
        [3.509e+00, 0])  # T = 1.0, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [9.558e-05, 3.320e-06],
        [6.822e+00, 0])  # T = 0.8, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.351e-04, 2.423e-06],
        [1.173e+01, 0])  # T = 0.6, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.843e-04, 2.755e-06],
        [1.754e+01, 0])  # T = 0.4, f = 0.01, vp = 6.0, sf = 0.25

    data_plot['leg'] = ['T_{EP} = [1.0, 0.8, 0.6, 0.4],\ v_{\parallel, EP} = 6.0,\ n_{EP}/n_e = 0.010']
    data_plot['sty'] = '*'
    data_plot['col'] = 'black'
    res_data.append(data_plot)

    # plot results: saturation levels of oscillating signals
    plot_several_data(res_data,
                      xlab='\gamma_{egam, lin}',
                      ylab='max.\ of\ \sqrt{<\overline{E}^2>_s}',
                      tit='Saturation\ levels')


def es_egam_energy_transfer_scan(oo):
    # NL ES EGAMb:  scan on EGAM growth rates
    dd = oo.get('dd', None)
    sel_norm_x = oo.get('sel_norm_x', 'none')

    # x, y normalization
    line_x, coef_x  = '', 1
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

    def plot_several_data(res_data_plot, xlab='x', ylab='y', tit='', ylim=None):
        curves = crv.Curves() \
            .xlab(xlab + line_x) \
            .ylab(ylab) \
            .tit(tit).xsty('plain').ylim(ylim)
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

    # ------------------------------------------------
    # --- TRANSFERRED ENERGY VS SATURATION LEVELS  ---
    # ------------------------------------------------
    res_data = []

    # scan on fpart
    data_plot = create_init_dict()

    add_f(
        [2.210e-05, 5.217e-06],
        [2.126e-03, 0])  # f = 0.0040, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [1.520e-04, 4.687e-06],
        [4.212e-02, 0])  # f = 0.0100, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [2.455e-04, 6.985e-06],
        [4.251e-01, 0])  # f = 0.0200, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [2.873e-04, 7.442e-06],
    #     [9.161e-01, 0])  # f = 0.0300, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [3.121e-04, 8.077e-06],
    #     [1.061e+00, 0])  # f = 0.0400, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [3.251e-04, 7.383e-06],
    #     [6.843e-01, 0])  # f = 0.0500, vp = 8.0, T = 1.0, sf = 0.25
    # add_f(
    #     [3.477e-04, 4.980e-05],
    #     [1.823e+00, 0])  # f = 0.0949, vp = 8.0, T = 1.0, sf = 0.25

    # data_plot['leg'] = ['scan\ on\ n_{EP}/n_e', 'v_{\parallel, EP} = 8.0,\ T = 1.0', '---------------']
    data_plot['leg'] = ['scan\ on\ n_{EP}/n_e,\ v_{\parallel, EP} = 8.0,\ T_{EP} = 1.0']
    data_plot['sty'] = 'o'
    data_plot['col'] = 'blue'
    res_data.append(data_plot)

    # scan on vp: f = 0.0100, T = 1.0, sf = 0.25
    data_plot = create_init_dict()

    add_f(
        [1.052e-04, 4.224e-06],
        [1.456e-02, 0])  # vp = 9.0
    add_f(
        [1.316e-04, 2.818e-06],
        [3.149e-02, 0])  # vp = 7.0
    add_f(
        [1.483e-04, 2.606e-06],
        [4.156e-02, 0])  # vp = 7.5

    # data_plot['leg'] = ['scan\ on\ v_{\parallel, EP}', 'n_{EP}/n_e = 0.01,\ T_{EP} = 1.0', '------------']
    data_plot['leg'] = ['scan\ on\ v_{\parallel, EP},\ n_{EP}/n_e = 0.01,\ T_{EP} = 1.0']
    data_plot['sty'] = 's'
    data_plot['col'] = 'red'
    res_data.append(data_plot)

    # scan on T, v = 3.5
    data_plot = create_init_dict()

    add_f(
        [9.063e-05, 2.111e-06],
        [1.057e-01, 0])  # f = 0.0945, vp = 3.5, T = 0.25, sf = 0.25
    add_f(
        [1.065e-04, 2.010e-06],
        [1.561e-01, 0])  # f = 0.0945, vp = 3.5, T = 0.22, sf = 0.25
    add_f(
        [1.173e-04, 2.003e-06],
        [1.927e-01, 0])  # f = 0.0945, vp = 3.5, T = 0.20, sf = 0.25
    add_f(
        [1.451e-04, 2.436e-06],
        [2.955e-01, 0])  # f = 0.0945, vp = 3.5, T = 0.15, sf = 0.25

    # data_plot['leg'] = ['scan\ on\ T_{EP}', 'v_{\parallel, EP} = 3.5,\ n_{EP}/n_e = 0.095', '---------------']
    data_plot['leg'] = ['scan\ on\ T_{EP},\ v_{\parallel, EP} = 3.5,\ n_{EP}/n_e = 0.095']
    data_plot['sty'] = '^'
    data_plot['col'] = 'green'
    res_data.append(data_plot)

    # scan on T, v = 6.0
    data_plot = create_init_dict()

    add_f(
        [6.304e-05, 2.952e-06],
        [4.203e-03, 0])  # T = 1.0, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [9.558e-05, 3.320e-06],
        [1.342e-02, 0])  # T = 0.8, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.351e-04, 2.423e-06],
        [3.507e-02, 0])  # T = 0.6, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.843e-04, 2.755e-06],
        [7.596e-02, 0])  # T = 0.4, f = 0.01, vp = 6.0, sf = 0.25

    data_plot['leg'] = ['scan\ on\ T_{EP},\ v_{\parallel, EP} = 6.0,\ n_{EP}/n_e = 0.010']
    data_plot['sty'] = '*'
    data_plot['col'] = 'black'
    res_data.append(data_plot)

    # plot results: saturation levels of oscillating signals
    plot_several_data(res_data,
                      xlab='\gamma_{egam, lin}',
                      ylab='\\int \\mathcal{P}_D dt\ \ \ [J/m^3]',
                      tit='Energy\ transferred\ to\ deuterium\ plasma')


def es_egam_energy_transfer_scan_diff_t(oo):
    # NL ES EGAMb:  scan on EGAM growth rates
    dd = oo.get('dd', None)
    sel_norm_x = oo.get('sel_norm_x', 'none')

    # x, y normalization
    line_x, coef_x  = '', 1
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

    def plot_several_data(res_data_plot, xlab='x', ylab='y', tit='', ylim=None):
        curves = crv.Curves() \
            .xlab(xlab + line_x) \
            .ylab(ylab) \
            .tit(tit).xsty('plain').ylim(ylim)
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

    # ------------------------------------------------
    # --- TRANSFERRED ENERGY VS SATURATION LEVELS  ---
    # ------------------------------------------------
    res_data = []

    T_speak_eV = dd['T_speak'] / constants.elementary_charge
    coef_norm_inv = 1. / (dd['cs'] * T_speak_eV / (dd['a0'] * dd['wc']))
    coef_norm_inv *= 1e2

    # scan on fpart
    data_plot = create_init_dict()

    add_f(
        [2.210e-05, 5.217e-06],
        [1.997e-03  * coef_norm_inv, 0])  # f = 0.0040, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [1.520e-04, 4.687e-06],
        [3.957e-02 * coef_norm_inv, 0])  # f = 0.0100, vp = 8.0, T = 1.0, sf = 0.25
    add_f(
        [2.455e-04, 6.985e-06],
        [2.068e-01 * coef_norm_inv, 0])  # f = 0.0200, vp = 8.0, T = 1.0, sf = 0.25
    data_plot['leg'] = ['scan\ on\ n_{EP}/n_e,\ v_{\parallel, EP} = 8.0,\ T_{EP} = 1.0']
    data_plot['sty'] = 'o'
    data_plot['col'] = 'blue'
    res_data.append(data_plot)

    # scan on vp: f = 0.0100, T = 1.0
    data_plot = create_init_dict()
    add_f(
        [1.052e-04, 4.224e-06],
        [1.046e-02 * coef_norm_inv, 0])  # vp = 9.0
    add_f(
        [1.316e-04, 2.818e-06],
        [2.695e-02 * coef_norm_inv, 0])  # vp = 7.0
    add_f(
        [1.483e-04, 2.606e-06],
        [3.535e-02 * coef_norm_inv, 0])  # vp = 7.5
    data_plot['leg'] = ['scan\ on\ v_{\parallel, EP},\ n_{EP}/n_e = 0.01,\ T_{EP} = 1.0']
    data_plot['sty'] = 's'
    data_plot['col'] = 'red'
    res_data.append(data_plot)

    # scan on T, v = 3.5
    data_plot = create_init_dict()
    add_f(
        [9.063e-05, 2.111e-06],
        [8.111e-02 * coef_norm_inv, 0])  # f = 0.0945, vp = 3.5, T = 0.25, sf = 0.25
    add_f(
        [1.065e-04, 2.010e-06],
        [1.154e-01 * coef_norm_inv, 0])  # f = 0.0945, vp = 3.5, T = 0.22, sf = 0.25
    add_f(
        [1.173e-04, 2.003e-06],
        [1.433e-01 * coef_norm_inv, 0])  # f = 0.0945, vp = 3.5, T = 0.20, sf = 0.25
    add_f(
        [1.451e-04, 2.436e-06],
        [2.252e-01 * coef_norm_inv, 0])  # f = 0.0945, vp = 3.5, T = 0.15, sf = 0.25
    data_plot['leg'] = ['scan\ on\ T_{EP},\ v_{\parallel, EP} = 3.5,\ n_{EP}/n_e = 0.095']
    data_plot['sty'] = '^'
    data_plot['col'] = 'green'
    res_data.append(data_plot)

    # scan on T, v = 6.0
    data_plot = create_init_dict()

    add_f(
        [6.304e-05, 2.952e-06],
        [3.978e-03 * coef_norm_inv, 0])  # T = 1.0, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [9.558e-05, 3.320e-06],
        [1.185e-02 * coef_norm_inv, 0])  # T = 0.8, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.351e-04, 2.423e-06],
        [3.030e-02 * coef_norm_inv, 0])  # T = 0.6, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.843e-04, 2.755e-06],
        [6.345e-02 * coef_norm_inv, 0])  # T = 0.4, f = 0.01, vp = 6.0, sf = 0.25

    data_plot['leg'] = ['scan\ on\ T_{EP},\ v_{\parallel, EP} = 6.0,\ n_{EP}/n_e = 0.010']
    data_plot['sty'] = '*'
    data_plot['col'] = 'black'
    res_data.append(data_plot)

    # plot results: saturation levels of oscillating signals
    plot_several_data(res_data,
                      xlab='\gamma_{egam, lin}',
                      ylab='\\int \\mathcal{P}_D dt (10^{-2})',
                      tit='t(ms)= [0.0, 1.0]', ylim=[-3e-1, 1.4e1])


def es_egam_satur_scan_overshoot_TALK(oo):
    # NL ES EGAMb:  scan on EGAM growth rates
    dd = oo.get('dd', None)
    sel_norm_x = oo.get('sel_norm_x', 'none')
    oo_texts = oo.get('text', [])

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

    def plot_several_data(res_data_plot, oo_texts, xlab='x', ylab='y', tit='', ylim=None):
        curves = crv.Curves() \
            .xlab(xlab + line_x) \
            .ylab(ylab) \
            .tit(tit).xsty('plain').ylim(ylim)
        # additional text:
        for oo_text in oo_texts:
            oText = crv.PlText(oo_text)
            curves.newt(oText)
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

    # ------------------------------------------------
    # --- SATURATION LEVELS OF OSCILLATING SIGNALS ---
    # ------------------------------------------------
    res_data = []

    # scan on T, v = 3.5
    data_plot = create_init_dict()

    add_f(
        [9.063e-05, 2.111e-06],
        [7.809e+00, 0])  # f = 0.0945, vp = 3.5, T = 0.25, sf = 0.25
    add_f(
        [1.065e-04, 2.010e-06],
        [9.929e+00, 0])  # f = 0.0945, vp = 3.5, T = 0.22, sf = 0.25
    add_f(
        [1.173e-04, 2.003e-06],
        [1.121e+01, 0])  # f = 0.0945, vp = 3.5, T = 0.20, sf = 0.25
    add_f(
        [1.451e-04, 2.436e-06],
        [1.636e+01, 0])  # f = 0.0945, vp = 3.5, T = 0.15, sf = 0.25

    data_plot['leg'] = None
    data_plot['sty'] = 'o'
    data_plot['col'] = 'blue'
    res_data.append(data_plot)

    # scan on T, v = 6.0
    data_plot = create_init_dict()

    add_f(
        [6.304e-05, 2.952e-06],
        [3.509e+00, 0])  # T = 1.0, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [9.558e-05, 3.320e-06],
        [6.822e+00, 0])  # T = 0.8, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.351e-04, 2.423e-06],
        [1.173e+01, 0])  # T = 0.6, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.843e-04, 2.755e-06],
        [1.754e+01, 0])  # T = 0.4, f = 0.01, vp = 6.0, sf = 0.25

    data_plot['leg'] = None
    data_plot['sty'] = '*'
    data_plot['col'] = 'red'
    res_data.append(data_plot)

    # plot results: saturation levels of oscillating signals
    plot_several_data(res_data, oo_texts,
                      xlab='\gamma_{egam, lin}',
                      ylab='max.\ of\ \sqrt{<\overline{E}^2>_s}',
                      tit=None)


def es_egam_energy_transfer_scan_diff_t_TALK(oo):
    # NL ES EGAMb:  scan on EGAM growth rates
    dd = oo.get('dd', None)
    sel_norm_x = oo.get('sel_norm_x', 'none')
    oo_texts = oo.get('text', [])

    # x, y normalization
    line_x, coef_x  = '', 1
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

    def plot_several_data(res_data_plot, oo_texts, xlab='x', ylab='y', tit='', ylim=None):
        curves = crv.Curves() \
            .xlab(xlab + line_x) \
            .ylab(ylab) \
            .tit(tit).xsty('plain').ylim(ylim)
        # additional text:
        for oo_text in oo_texts:
            oText = crv.PlText(oo_text)
            curves.newt(oText)
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

    # ------------------------------------------------
    # --- TRANSFERRED ENERGY VS SATURATION LEVELS  ---
    # ------------------------------------------------
    res_data = []

    # scan on T, v = 3.5
    data_plot = create_init_dict()
    add_f(
        [9.063e-05, 2.111e-06],
        [8.111e-02, 0])  # f = 0.0945, vp = 3.5, T = 0.25, sf = 0.25
    add_f(
        [1.065e-04, 2.010e-06],
        [1.154e-01, 0])  # f = 0.0945, vp = 3.5, T = 0.22, sf = 0.25
    add_f(
        [1.173e-04, 2.003e-06],
        [1.433e-01, 0])  # f = 0.0945, vp = 3.5, T = 0.20, sf = 0.25
    add_f(
        [1.451e-04, 2.436e-06],
        [2.252e-01, 0])  # f = 0.0945, vp = 3.5, T = 0.15, sf = 0.25
    data_plot['leg'] = None
    data_plot['sty'] = 'o'
    data_plot['col'] = 'blue'
    res_data.append(data_plot)

    # scan on T, v = 6.0
    data_plot = create_init_dict()

    add_f(
        [6.304e-05, 2.952e-06],
        [3.978e-03, 0])  # T = 1.0, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [9.558e-05, 3.320e-06],
        [1.185e-02, 0])  # T = 0.8, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.351e-04, 2.423e-06],
        [3.030e-02, 0])  # T = 0.6, f = 0.01, vp = 6.0, sf = 0.25
    add_f(
        [1.843e-04, 2.755e-06],
        [6.345e-02, 0])  # T = 0.4, f = 0.01, vp = 6.0, sf = 0.25

    data_plot['leg'] = None
    data_plot['sty'] = '*'
    data_plot['col'] = 'red'
    res_data.append(data_plot)

    # plot results: saturation levels of oscillating signals
    plot_several_data(res_data, oo_texts,
                      xlab='\gamma_{egam, lin}',
                      ylab='\\int \\mathcal{P}_D dt\ \ \ [J/m^3]',
                      tit=None)


def egam_fpart_AE_KE_g_mpr_tot_TALK(oo):
    # EGAMb: adiabatic vs drift-kinetic electrons: rho_f = 0.25, v|| = 8.0, T = 1.0:
    # SCAN ON n_{EP}/n_{e}

    dd = oo.get('dd', None)
    sel_norm = oo.get('sel_norm', 'wc')
    oo_texts = oo.get('text', [])

    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    data_plot = {
        'f': [],
        'w': [], 'w_err': [],
        'g': [], 'g_err': [],
        'leg': [], 'sty': [], 'col': []
    }

    # ws[i] and gs[i] should have the same size
    def add_f(f, w, w_err, g, g_err):
        fs.append(f)
        ws.append(w)
        ws_err.append(w_err)
        gs.append(g)
        gs_err.append(g_err)

    # --- ADIAB. ELE.: MPR TOTAL: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0060,
          [np.nan], [np.nan],
          [7.735e-05], [6.706e-06])
    add_f(0.0100,
          [np.nan], [np.nan],
          [1.503e-04], [4.945e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [2.412e-04], [3.550e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [2.845e-04], [2.703e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [3.165e-04], [2.586e-06])
    add_f(0.0700,
          [np.nan], [np.nan],
          [3.135e-04], [7.597e-06])
    add_f(0.0900,
          [np.nan], [np.nan],
          [3.116e-04], [2.222e-05])

    # data_plot['leg'].append('AE:\ total')
    data_plot['leg'].append(None)
    data_plot['sty'].append('o')
    data_plot['col'].append('blue')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # --- DK. ELE: MPR TOTAL: w, g are normalized to wci ---
    fs, ws, ws_err, gs, gs_err = [], [], [], [], []
    add_f(0.0100,
          [np.nan],    [np.nan],
          [5.621e-05], [4.360e-06])
    add_f(0.0200,
          [np.nan], [np.nan],
          [1.287e-04], [2.505e-06])
    add_f(0.0300,
          [np.nan], [np.nan],
          [1.549e-04], [4.018e-06])
    add_f(0.0500,
          [np.nan], [np.nan],
          [1.568e-04], [1.437e-05])
    add_f(0.0700,
          [np.nan], [np.nan],
          [1.470e-04], [1.528e-05])
    add_f(0.0900,
          [np.nan], [np.nan],
          [1.582e-04], [1.171e-05])

    # data_plot['leg'].append('KE:\ total')
    data_plot['leg'].append(None)
    data_plot['sty'].append('*')
    data_plot['col'].append('red')
    reorganise_data(fs, ws, ws_err, gs, gs_err, sel_norm, dd, data_plot)

    # *** combined plot ***
    plot_several_scans_adv(data_plot, sel_norm, oo_texts, xlim=[0.0, 0.1],
                           ylab_g=None, tit_w='EGAM\ frequency')


