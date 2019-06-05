import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import gam_theory
import general as gn
import numpy as np
import copy
from scipy import interpolate
from scipy import signal


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(gam_theory)
    mix.reload_module(gn)


def erbar(dd):
    if 'erbar' in dd:
        return
    rd.phibar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']
    data = - np.gradient(dd['phibar']['data'], s, axis=1)
    dd['erbar'] = {
        't': t,
        's': s,
        'data': data}


def vorbar(dd):
    # shearing rate, vorticity
    if 'vorbar' in dd:
        return
    erbar(dd)
    t = dd['erbar']['t']
    s = dd['erbar']['s']

    data = dd['erbar']['data'] * s
    data = np.gradient(data, s, axis=1)
    if s[0] > 0:
        data = data / s
    else:
        data[:, 1:-1] = data[:, 1:-1] / s[1:-1]
        f = interpolate.interp2d(s[1:-1], t, data[:, 1:-1], kind='cubic')
        data = f(s, t)

    dd['vorbar'] = {
        't': t,
        's': s,
        'data': data}


def choose_var(dd, oo):
    vorbar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']

    opt_var = oo.get('opt_var', 'erbar')
    vvar = dd[opt_var]['data']

    tit_var = ''
    if opt_var == 'phibar':
        tit_var = '\overline{\Phi}'
    if opt_var == 'erbar':
        tit_var = '\overline{E}_r'
    if opt_var == 'vorbar':
        tit_var = '\overline{\Omega}_r'

    res = {
        'var': vvar,
        's': s,
        't': t,
        'tit': tit_var
    }

    return res


# NEW: take signal (t,s)
def choose_one_var_ts(ovar, dd):
    vorbar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']

    opt_var = ovar[0]

    vvar = dd[opt_var]['data']

    tit_var = ''
    if opt_var == 'phibar':
        tit_var = '\overline{\Phi}'
    if opt_var == 'erbar':
        tit_var = '\overline{E}_r'
    if opt_var == 'vorbar':
        tit_var = '\overline{\Omega}_r'

    res = {
        'data': vvar,
        's': s,
        't': t,
        'tit': tit_var
    }

    return res


def plot_fft(dd, oo):
    out = choose_var(dd, oo)
    vvar, r, t, tit_var = out['var'], out['s'], out['t'], out['tit']

    # radial coordinate normalization
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'
    line_r = ''
    if sel_r == 's':
        r = r
        line_r = 's = \sqrt{\psi/\psi_{edge}}'
    if sel_r == 'psi':
        r = r ** 2
        line_r = '\psi/\psi_{edge}'

    # plotting
    oo_fft = dict(oo)
    oo_fft.update({
        'var': vvar, 't_wci': t, 'r': r,
        'tit': tit_var,
        'labr': line_r
    })
    gn.plot_fft(dd, oo_fft)


def plot_fft_1d(dd, oo):
    out = choose_var(dd, oo)
    vvar, s, t, tit_var = out['var'], out['s'], out['t'], out['tit']

    # plotting
    oo_fft = dict(oo)
    oo_fft.update({
        'vars': [vvar], 'ts_wci': [t], 'ss': [s],
        'tit': tit_var,
        'vars_names': ['']
    })
    gn.plot_fft_1d(dd, oo_fft)


def find_gamma(dd, oo):
    out = choose_var(dd, oo)
    vvar, s, t, tit_var = out['var'], out['s'], out['t'], out['tit']

    oo_gamma = oo  # oo must have some 'opt_av'
    oo_gamma.update({
        'var': vvar, 't': t, 's': s,
        'tit': tit_var, 'var_name': ''
    })
    gn.find_gamma(dd, oo_gamma)


# time evolution of zonal Phi, Er, Vorticity at different radial points:
def t_evol(dd, oo):
    sel_norm = oo.get('sel_norm', 'wci')

    vorbar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']

    # time normalization
    if sel_norm == 'ms':
        coef_norm = 1.e3 / dd['wc']
        line_t = 't,\ ms'
    if sel_norm == 'wci':
        coef_norm = 1
        line_t = 't[\omega_c^{-1}]'
    if sel_norm == 'csa':
        coef_norm = (dd['cs'] / dd['a0']) / dd['wc']
        line_t = 't[a_0/c_s]'
    if sel_norm == 'csr':
        coef_norm = (dd['cs'] / dd['R0']) / dd['wc']
        line_t = 't[R_0/c_s]'
    t = t * coef_norm

    # radial interval
    t, ids_t = mix.get_array_oo(oo, t, 't')

    # time points:
    ss = np.array(oo.get('s_points', [0.5]))
    ids_s = np.zeros(np.size(ss))
    for id_s in range(np.size(ss)):
        ids_s[id_s], ss[id_s] = mix.find(s, ss[id_s])

    # plotting
    curves_phi = crv.Curves().xlab(line_t).ylab('\overline{\Phi}') \
        .tit('\overline{\Phi}')
    curves_er = crv.Curves().xlab(line_t).ylab('\overline{E}_r') \
        .tit('\overline{E}_r')
    curves_vor = crv.Curves().xlab(line_t).ylab('\overline{\Omega}_r') \
        .tit('\overline{\Omega}_r')

    for count_s in range(np.size(ids_s)):
        id_s = int(ids_s[count_s])
        s1 = ss[count_s]
        curves_phi.new('{:d}'.format(count_s)) \
            .XS(t) \
            .YS(dd['phibar']['data'][ids_t[0]:ids_t[-1] + 1, id_s]) \
            .leg('s = {:0.3e}'.format(s1)).new_sty(count_s)
        curves_er.new('{:d}'.format(count_s)) \
            .XS(t) \
            .YS(dd['erbar']['data'][ids_t[0]:ids_t[-1] + 1, id_s]) \
            .leg('s = {:0.3e}'.format(s1)).new_sty(count_s)
        curves_vor.new('{:d}'.format(count_s)) \
            .XS(t) \
            .YS(dd['vorbar']['data'][ids_t[0]:ids_t[-1] + 1, id_s]) \
            .leg('s = {:0.3e}'.format(s1)).new_sty(count_s)

    cpr.plot_curves(curves_phi)
    cpr.plot_curves(curves_er)
    cpr.plot_curves(curves_vor)


# frequency(s) and dynamic rate(s) of zonal Er
def wg_s(dd, oo):
    flag_save = oo.get('flag_save', False)
    save_name = oo.get('save_name', None)
    sel_norm = oo.get('sel_norm', 'wci')

    erbar(dd)
    s = dd['phibar']['s']
    s, _ = mix.get_array_oo(oo, s, 's')
    ns = np.size(s)
    ws_est, gs_est, ws_adv, gs_adv = np.zeros(ns), np.zeros(ns), np.zeros(ns), np.zeros(ns)
    line_w, line_g, out = None, None, None
    for id_s in range(ns):
        oo_s1 = dict(oo)
        oo_s1.update({'flag_output': True, 's_point': s[id_s], 'flag_print': False})
        out = wg_s1(dd, oo_s1)
        if id_s == 0:
            line_w = out['line_w']
            line_g = out['line_g']
        if out is None:
            ws_est[id_s], ws_adv[id_s] = np.NaN, np.NaN
            gs_est[id_s], gs_adv[id_s] = np.NaN, np.NaN
        else:
            # with the corresponding normalization
            ws_est[id_s], ws_adv[id_s] = out['w-est'], out['w-adv']
            gs_est[id_s], gs_adv[id_s] = out['g-est'], out['g-adv']
    del out

    # parameters for the analytical and experimental data:
    rd.krpert(dd, 'pf')
    kr = dd['pf'].krpert * np.pi / dd['Lwork']
    line_k = 'k = {:0.3f}'.format(kr * dd['rhoL_speak'])
    oo_th = {'sel_norm': sel_norm, 'r': s, 'kr': kr}

    # --- frequencies ---
    curves_w = crv.Curves().xlab('s').ylab(line_w)\
        .tit('\omega\ of\ \overline{E}_r') \
        .xlim([s[0], s[-1]])

    curves_w.new('f')\
        .XS(s).YS(ws_est)\
        .leg('Estimation:\ Lin.\ GK\ ' + line_k)\
        .sty('o').col('red').mfc('red')
    curves_w.new('f')\
        .XS(s).YS(ws_adv)\
        .leg('Fitting:\ Lin.\ GK\ ' + line_k) \
        .sty('o').col('green').mfc('green')

    oo_th.update({'curves': curves_w, 'col': 'green'})
    curves_w = gam_theory.get_gao(dd, oo_th)
    # oo_th.update({'col': 'green'})
    # curves_w = gam_theory.get_gk_fit(dd, oo_th)

    cpr.plot_curves(curves_w)

    # # --- dynamic rate ---
    # curves_g = crv.Curves().xlab('s').ylab(line_g) \
    #     .tit('\gamma\ of\ \overline{E}_r:\ ') \
    #     .xlim([s[0], s[-1]])
    #
    # curves_g.new('f').XS(s).YS(gs).leg('Lin.\ GK\ ' + line_k)\
    #     .sty('o').col('red').mfc('red')
    #
    # # oo_th.update({'curves': curves_g, 'sel_res': 'g'})
    # # curves_g = gam_theory.get_gao(dd, oo_th)
    # # curves_g = gam_theory.get_gk_fit(dd, oo_th)
    #
    # cpr.plot_curves(curves_g)

    # save data
    if flag_save:
        if 'saved_data' not in dd:
            dd['saved_data'] = {}
        dd['saved_data'][save_name + '-w'] = curves_w
        # dd['saved_data'][save_name + '-g'] = curves_g


# w,g in several s-points:
def wg_several_s1(dd, oo):
    s_points = oo.get('s_points', [0.5])

    oo_points = dict(oo)
    oo_points['flag_output'] = True
    out_s = []
    for s_point in s_points:
        oo_points['s_point'] = s_point
        out_s.append(wg_s1(dd, oo_points))

    print('--- RESULTS ---')
    for out_s1 in out_s:
        print('-- ' + out_s1['line_s'] + ', ' + out_s1['line_t'] + ' --')
        print('E -> ' + out_s1['line_w_pr'] +
              '{:0.3e}'.format(out_s1['w-est']))
        print('E -> ' + out_s1['line_g_pr'] +
              '{:0.3e}'.format(out_s1['g-est']))
        print('A -> ' + out_s1['line_w_pr'] +
              '{:0.3e}'.format(out_s1['w-adv']))
        print('A -> ' + out_s1['line_g_pr'] +
              '{:0.3e}'.format(out_s1['g-adv']))


# frequency and dynamic rate of zonal Er at a radial point s1
def wg_s1(dd, oo):
    sel_norm = oo.get('sel_norm', 'wci')  # -> 'wci', 'kHz', 'csa', 'csr'
    flag_output = oo.get('flag_output', False)
    flag_print = oo.get('flag_print', True)
    s_point = oo.get('s_point', 0.5)
    oo_filter = oo.get('filter', [None])

    erbar(dd)
    t = dd['phibar']['t']  # normalized to wc
    s = dd['phibar']['s']

    # normalization:
    coef_norm_w, coef_norm_g = np.NaN, np.NaN
    line_w, line_w_pr = '', ''
    line_g, line_g_pr = '', ''
    if sel_norm == 'khz':
        coef_norm_w = dd['wc'] / (2 * np.pi * 1.e3)
        coef_norm_g = dd['wc'] / 1.e3
        line_w, line_w_pr = '\omega,\ kHz',   'w(kHz) = '
        line_g, line_g_pr = '\gamma,\ 1e3/s', 'g(1e3/s) = '
    if sel_norm == 'wci':
        coef_norm_w = coef_norm_g = 1
        line_w, line_w_pr = '\omega[\omega_c]', 'w[wc] = '
        line_g, line_g_pr = '\gamma[\omega_c]', 'g[wc] = '
    if sel_norm == 'csa':
        coef_norm_w = coef_norm_g = dd['wc'] / (dd['cs'] / dd['a0'])
        line_w, line_w_pr = '\omega[c_s/a_0]', 'w[cs/a0] = '
        line_g, line_g_pr = '\gamma[c_s/a_0]', 'g[cs/a0] = '
    if sel_norm == 'csr':
        coef_norm_w = coef_norm_g = dd['wc'] / (dd['cs'] / dd['R0'])
        line_w, line_w_pr = '\omega[c_s/R_0]', 'w[cs/R0] = '
        line_g, line_g_pr = '\gamma[c_s/R_0]', 'g[cs/R0] = '

    # radial interval
    ids_s1, s_point = mix.find(s, s_point)
    line_s1 = 's = {:0.3f}'.format(s_point)

    # zonal radial electric field at the chosen radial point:
    erbar_s1 = dd['erbar']['data'][:, ids_s1]

    # radial wavelength:
    rd.krpert(dd, 'pf')
    kr = dd['pf'].krpert * np.pi / dd['Lwork']
    line_k = 'k = {:0.3f}'.format(kr * dd['rhoL_speak'])

    # filtering of the signal:
    filt = ymath.filtering(t, erbar_s1, oo_filter)

    # plot filtered signals:
    if flag_print:
        curves = crv.Curves()\
            .xlab('t[\omega_c^{-1}]')\
            .ylab('\overline{E}_r')\
            .tit('\overline{E}_r:\ ' + line_s1)
        curves.new()\
            .XS(t)\
            .YS(erbar_s1)\
            .leg('initial')
        curves.new()\
            .XS(filt['x'])\
            .YS(filt['filt']) \
            .leg('filtered').sty(':')
        cpr.plot_curves(curves)

        curves = crv.Curves()\
            .xlab(line_w)\
            .ylab('FFT:\ \overline{E}_r')\
            .tit('FFT:\ \overline{E}_r:\ ' + line_s1)
        curves.new('y')\
            .XS(2 * np.pi * filt['w2'] * coef_norm_w)\
            .YS(filt['fft_init_2'])\
            .leg('initial')
        curves.new('filt')\
            .XS(2 * np.pi * filt['w2'] * coef_norm_w)\
            .YS(filt['fft_filt_2']) \
            .leg('filtered').sty(':')
        cpr.plot_curves(curves)

    # chose time domain to calculate the frequency and damping rate:
    t_work, ids_t = mix.get_array_oo(oo, t, 't')
    line_t_wci = 't = [{:0.2e}, {:0.2e}]'.format(t_work[0], t_work[-1])
    erbar_s1_work = filt['filt'][ids_t[0]:ids_t[-1]+1]

    # --------------------------------------------------------------
    # --- ESTIMATION ---
    wg_est = ymath.estimate_wg(t_work, erbar_s1_work)
    if wg_est is None:
        return None
    w_est = wg_est['w'] * coef_norm_w
    g_est = wg_est['g'] * coef_norm_g

    # plotting of the fitting
    if flag_print:
        curves = crv.Curves().xlab('t[wci^{-1}]').ylab('\overline{E}_r')\
            .tit('\overline{E}_r:\ ' + line_s1 + ':\ ' + line_k)\
            .titn(line_w + ' = {:0.3e}'.format(w_est))\
            .titn(line_g + ' = {:0.3e}'.format(g_est))
        curves.flag_semilogy = True
        curves.new('filt')\
            .XS(t_work)\
            .YS(erbar_s1_work)\
            .leg('init')
        curves.new('peaks')\
            .XS(wg_est['x_peaks'])\
            .YS(wg_est['y_peaks'])\
            .leg('peaks').sty('o').col('green')
        curves.new('fitting') \
            .XS(wg_est['x_fit']) \
            .YS(wg_est['y_fit']) \
            .leg('fitting').col('red').sty('--')
        cpr.plot_curves(curves)

        # print the results
        print('--- ESTIMATION: ' + line_t_wci + ', ' + line_s1 + '---')
        print(line_w_pr + '{:0.3e}'.format(w_est))
        print(line_g_pr + '{:0.3e}'.format(g_est))

    # --------------------------------------------------------------
    # --- ADVANCED FITTING ---
    ainf = {'est': wg_est,
            'x_start': wg_est['x_peaks'][0],
            'x_end':   wg_est['x_peaks'][-1]
            }
    wg_adv = ymath.advanced_wg(t_work, erbar_s1_work, ainf)
    if wg_adv is None:
        return None
    w_adv = wg_adv['w'] * coef_norm_w
    g_adv = wg_adv['g'] * coef_norm_g

    # plotting of the fitting
    if flag_print:
        curves = crv.Curves().xlab('t[wci^{-1}]').ylab('\overline{E}_r') \
            .tit('\overline{E}_r:\ ' + line_s1 + ':\ ' + line_k) \
            .titn(line_w + ' = {:0.3e}'.format(w_adv)) \
            .titn(line_g + ' = {:0.3e}'.format(g_adv))
        curves.flag_semilogy = True
        curves.new('filt') \
            .XS(t_work) \
            .YS(erbar_s1_work) \
            .leg('init')
        curves.new('fitting') \
            .XS(wg_adv['x_fit']) \
            .YS(wg_adv['y_fit']) \
            .leg('adv.\ fitting').col('red').sty('--')
        cpr.plot_curves(curves)

        # print the results
        print('--- ADVANCED FITTING: ' + line_t_wci + ', ' + line_s1 + '---')
        print(line_w_pr + '{:0.3e}'.format(w_adv))
        print(line_g_pr + '{:0.3e}'.format(g_adv))

    # results
    if flag_output:
        out = {'w-est': w_est, 'g-est': g_est,
               'w-adv': w_adv, 'g-adv': g_adv,
               'line_w': line_w, 'line_g': line_g,
               'line_w_pr': line_w_pr, 'line_g_pr': line_g_pr,
               'line_s': line_s1, 'line_t': line_t_wci}
        return out

















