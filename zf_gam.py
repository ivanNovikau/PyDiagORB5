import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import gam_theory
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


def fft_gam_s1(dd, ss, oo={}):
    # Fourier spectra at different radial points ss:
    rd.phibar(dd)
    erbar(dd)
    vorbar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']

    # intervals
    t, ids_t = mix.get_array_oo(oo, t, 't')

    phibar_s = []
    erbar_s = []
    vorbar_s = []
    for s1 in ss:
        id_s1, _ = mix.find(s, s1)
        phibar_s.append(mix.get_slice(dd['phibar']['data'], ids_t, id_s1))
        erbar_s.append(mix.get_slice(dd['erbar']['data'], ids_t, id_s1))
        vorbar_s.append(mix.get_slice(dd['vorbar']['data'], ids_t, id_s1))

    # change normalization of the time grid
    cs = dd['cs']
    wc = dd['wc']
    a0 = dd['a0']
    t = t * (cs / a0) / wc

    # find FFT of the signal:
    w, _, nw = ymath.w_for_fft(t)
    w = w * 2 * np.pi

    fs_phibar = []
    fs_erbar = []
    fs_vorbar = []
    for id_s1 in range(len(ss)):
        f_phibar, _ = ymath.fft_y(phibar_s[id_s1], nw)
        f_erbar,  _ = ymath.fft_y(erbar_s[id_s1], nw)
        f_vorbar, _ = ymath.fft_y(vorbar_s[id_s1], nw)
        fs_phibar.append(f_phibar)
        fs_erbar.append(f_erbar)
        fs_vorbar.append(f_vorbar)
    del f_phibar, f_erbar, f_vorbar

    # intervals for the frequency:
    w, ids_w = mix.get_array_oo(oo, w, 'w')
    for id_s1 in range(len(ss)):
        fs_phibar[id_s1] = mix.get_slice(fs_phibar[id_s1], ids_w)
        fs_erbar[id_s1]  = mix.get_slice(fs_erbar[id_s1],  ids_w)
        fs_vorbar[id_s1] = mix.get_slice(fs_vorbar[id_s1], ids_w)

    # initial signal
    curves_phi = crv.Curves().xlab('t[a/cs]').ylab('\overline{\Phi}')
    curves_er = crv.Curves().xlab('t[a/cs]').ylab('\overline{E}_r')
    curves_vor = crv.Curves().xlab('t[a/cs]').ylab('\overline{\Omega}_r')
    for id_s1 in range(len(ss)):
        name_y = 'y{:d}'.format(id_s1)
        leg_s1 = 's = {:0.3f}'.format(ss[id_s1])
        curves_phi.new(name_y).XS(t).YS(phibar_s[id_s1])\
            .leg('\overline{\Phi}' + '({})'.format(leg_s1))
        curves_er.new(name_y).XS(t).YS(erbar_s[id_s1])\
            .leg('\overline{E}_r' + '({})'.format(leg_s1))
        curves_vor.new(name_y).XS(t).YS(vorbar_s[id_s1])\
            .leg('\overline{\Omega}_r' + '({})'.format(leg_s1))
    cpr.plot_curves(curves_phi)
    cpr.plot_curves(curves_er)
    cpr.plot_curves(curves_vor)

    # fourier spectrum
    curves_phi = crv.Curves().xlab('w[cs/a]').ylab('FFT:\ \overline{\Phi}')
    curves_er = crv.Curves().xlab('w[cs/a]').ylab('FFT:\ \overline{E}_r')
    curves_vor = crv.Curves().xlab('w[cs/a]').ylab('FFT:\ \overline{\Omega}_r')
    for id_s1 in range(len(ss)):
        name_y = 'fft:\ y{:d}'.format(id_s1)
        leg_s1 = 's = {:0.3f}'.format(ss[id_s1])
        curves_phi.new(name_y).XS(w).YS(fs_phibar[id_s1])\
            .leg('FFT:\ \overline{\Phi}' + '({})'.format(leg_s1))
        curves_er.new(name_y).XS(w).YS(fs_erbar[id_s1])\
            .leg('FFT:\ \overline{E}_r' + '({})'.format(leg_s1))
        curves_vor.new(name_y).XS(w).YS(fs_vorbar[id_s1])\
            .leg('FFT:\ \overline{\Omega}_r' + '({})'.format(leg_s1))
    cpr.plot_curves(curves_phi)
    cpr.plot_curves(curves_er)
    cpr.plot_curves(curves_vor)


# FFT of zonal radial eletric field and vorticity as (s,t)
def fft_gam_2d(dd, oo={}):
    sel_norm = oo.get('sel_norm', 'wci')  # -> 'wci', 'kHz', 'csa', 'csr'
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'
    sel_cmp = oo.get('sel_cmp', 'pink_r')  # -> 'jet', 'hot' etc.
    curves_to_load = oo.get('curves_to_load', None)

    rd.phibar(dd)
    erbar(dd)
    vorbar(dd)
    t = dd['phibar']['t']
    r = dd['phibar']['s']

    # time interval
    t, ids_t = mix.get_array_oo(oo, t, 't')

    # radial coordinate and its normalization
    if sel_r == 's':
        r = r
        line_r = '\sqrt{\psi/\psi_{edge}}'
    if sel_r == 'psi':
        r = r**2
        line_r = '\psi/\psi_{edge}'
    r, ids_r = mix.get_array_oo(oo, r, sel_r)

    # zonal radial electric field
    erbar_st  = mix.get_slice(dd['erbar']['data'],  ids_t, ids_r)

    # zonal vorticity
    vorbar_st = mix.get_slice(dd['vorbar']['data'], ids_t, ids_r)

    # information about the time-interval where FFT is found
    line_t_wci = mix.test_array(t, 't[wci^{-1}]', ':0.2e')

    # change normalization of the time grid -> seconds
    t = t / dd['wc']

    # frequency grid with a chosen normalization
    w, _, nw = ymath.w_for_fft(t)
    if sel_norm == 'khz':
        coef_norm = 1. / 1.e3
        line_w = '\omega,\ kHz'
    if sel_norm == 'wci':
        coef_norm = 2 * np.pi / dd['wc']
        line_w = '\omega[\omega_c]'
    if sel_norm == 'csa':
        coef_norm = 2 * np.pi / (dd['cs'] / dd['a0'])
        line_w = '\omega[c_s/a_0]'
    if sel_norm == 'csr':
        coef_norm = 2 * np.pi / (dd['cs'] / dd['R0'])
        line_w = '\omega[c_s/R_0]'
    w = w * coef_norm

    # FFT of the signal
    f_erbar,  _ = ymath.fft_y(erbar_st,  nw, 0)
    f_vorbar, _ = ymath.fft_y(vorbar_st, nw, 0)

    # intervals for the frequency:
    w, ids_w = mix.get_array_oo(oo, w, 'w')
    f_erbar  = mix.get_slice(f_erbar,  ids_w)
    f_vorbar = mix.get_slice(f_vorbar, ids_w)

    # --- PLOTTING ---
    curves_er = crv.Curves().xlab(line_r).ylab(line_w)\
        .tit('FFT:\ \overline{E}_r:\ ' + line_t_wci)\
        .xlim([r[0], r[-1]]).ylim([w[0], w[-1]]) \
        .leg_pos('lower left')
    curves_vor = crv.Curves().xlab(line_r).ylab(line_w)\
        .tit('FFT:\ \overline{\Omega}_r:\ ' + line_t_wci)\
        .xlim([r[0], r[-1]]).ylim([w[0], w[-1]])\
        .leg_pos('lower left')

    curves_er.new('f').XS(r).YS(w).ZS(f_erbar.T).lev(20).cmp(sel_cmp)
    curves_vor.new('f').XS(r).YS(w).ZS(f_vorbar.T).lev(20).cmp(sel_cmp)

    # parameters for the analytical and experimental data:
    oo_th = {'curves': curves_er, 'sel_norm': sel_norm,
             'sel_r': sel_r, 's': r, 'col': 'red'}

    curves_er = exp_AUG20787(dd, oo_th)
    curves_er = gam_theory.get_gao(dd, oo_th)
    curves_er = gam_theory.get_gk_fit(dd, oo_th)

    oo_th.update({'curves': curves_vor}) # !!!
    curves_vor = exp_AUG20787(dd, oo_th)
    curves_vor = gam_theory.get_gao(dd, oo_th)
    curves_vor = gam_theory.get_gk_fit(dd, oo_th)

    # load saved data:
    curves_er.load(curves_to_load)
    curves_vor.load(curves_to_load)

    # plot curves
    cpr.plot_curves_3d(curves_er)
    cpr.plot_curves_3d(curves_vor)

    return


# Zonal Er, Vorticity (s,t)
def time_evol_st(dd, oo={}):
    vorbar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']

    # intervals
    t, ids_t = mix.get_array_oo(oo, t, 't')
    s, ids_s = mix.get_array_oo(oo, s, 's')
    # phibar_st = mix.get_slice(dd['phibar']['data'], ids_t, ids_s)
    erbar_st = mix.get_slice(dd['erbar']['data'], ids_t, ids_s)
    vorbar_st = mix.get_slice(dd['vorbar']['data'], ids_t, ids_s)

    # plotting
    # curves_phi = crv.Curves().xlab('t[wci^{-1}]').ylab('s') \
    #     .tit('\overline{\Phi}')
    curves_er = crv.Curves().xlab('t[wci^{-1}]').ylab('s') \
        .tit('\overline{E}_r')
    curves_vor = crv.Curves().xlab('t[wci^{-1}]').ylab('s') \
        .tit('\overline{\Omega}_r')
    # curves_phi.new('y').XS(t).YS(s).ZS(phibar_st)
    curves_er.new('y').XS(t).YS(s).ZS(erbar_st).lev(60)
    curves_vor.new('y').XS(t).YS(s).ZS(vorbar_st).lev(60)
    # cpr.plot_curves_3d(curves_phi)
    cpr.plot_curves_3d(curves_er)
    cpr.plot_curves_3d(curves_vor)


# radial strcuture of zonal Phi, Er, Vorticity at different time points:
def rad_structure(dd, oo={}):
    vorbar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']

    # radial interval
    s, ids_s = mix.get_array_oo(oo, s, 's')

    # time points:
    ts = np.array(oo.get('t_points', [0.0]))
    ids_t = np.zeros(np.size(ts))
    for id_t in range(np.size(ts)):
        t1 = ts[id_t]
        ids_t[id_t], ts[id_t] = mix.find(t, t1)

    # plotting
    curves_phi = crv.Curves().xlab('s').ylab('\overline{\Phi}') \
        .tit('\overline{\Phi}')
    curves_er = crv.Curves().xlab('s').ylab('\overline{E}_r') \
        .tit('\overline{E}_r')
    curves_vor = crv.Curves().xlab('s').ylab('\overline{\Omega}_r') \
        .tit('\overline{\Omega}_r')

    for it in range(np.size(ids_t)):
        id_t = int(ids_t[it])
        curves_phi.new('{:d}'.format(it))\
            .XS(s)\
            .YS(dd['phibar']['data'][id_t, ids_s[0]:ids_s[-1]+1])\
            .leg('t = {:0.3e}'.format(ts[it])).new_sty(it)
        curves_er.new('{:d}'.format(it)) \
            .XS(s) \
            .YS(dd['erbar']['data'][id_t, ids_s[0]:ids_s[-1]+1]) \
            .leg('t = {:0.3e}'.format(ts[it])).new_sty(it)
        curves_vor.new('{:d}'.format(it)) \
            .XS(s) \
            .YS(dd['vorbar']['data'][id_t, ids_s[0]:ids_s[-1]+1]) \
            .leg('t = {:0.3e}'.format(ts[it])).new_sty(it)

    cpr.plot_curves(curves_phi)
    cpr.plot_curves(curves_er)
    cpr.plot_curves(curves_vor)


# time evolution of zonal Phi, Er, Vorticity at different radial points:
def t_evol(dd, oo={}):
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
def wg_s(dd, oo={}):
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


# frequency and dynamic rate of zonal Er at a radial point s1
def wg_s1(dd, oo={}):
    sel_norm = oo.get('sel_norm', 'wci')  # -> 'wci', 'kHz', 'csa', 'csr'
    s_point = oo.get('s_point', 0.5)
    w_interval = oo.get('w_interval', None)
    flag_output = oo.get('flag_output', False)
    flag_print = oo.get('flag_print', True)

    erbar(dd)
    t = dd['phibar']['t']  # normalized to wc
    s = dd['phibar']['s']

    # radial interval
    ids_s1, s_point = mix.find(s, s_point)
    line_s1 = 's = {:0.3f}'.format(s_point)

    # radial wavelength:
    rd.krpert(dd, 'pf')
    kr = dd['pf'].krpert * np.pi / dd['Lwork']
    line_k = 'k = {:0.3f}'.format(kr * dd['rhoL_speak'])

    # filtering of the signal:
    t_filt, ids_t_filt = mix.get_array_oo(oo, t, 't_filt')
    erbar_s1 = mix.get_slice(dd['erbar']['data'], ids_t_filt, ids_s1)
    filt = None
    if w_interval is not None:
        filt = ymath.rough_filter(t_filt, erbar_s1, {'w_interval': w_interval})
        erbar_s1_filt = filt['filt']
    else:
        erbar_s1_filt = erbar_s1
    del ids_t_filt, ids_s1, t

    # plot filtered signals:
    if flag_print and w_interval is not None:
        curves = crv.Curves().xlab('t[1/wci]').ylab('\overline{E}_r')
        curves.new('y').XS(t_filt).YS(erbar_s1).leg('init')
        curves.new('filt').XS(t_filt).YS(erbar_s1_filt) \
            .leg('filt').sty(':')
        cpr.plot_curves(curves)

        curves = crv.Curves().xlab('w[wci]').ylab('FFT:\ \overline{E}_r')
        curves.new('y').XS(filt['w2']).YS(filt['fft_init_2']).leg('init')
        curves.new('filt').XS(filt['w2']).YS(filt['fft_filt_2']) \
            .leg('filt').sty(':')
        cpr.plot_curves(curves)

    # chose time domain to calculate the frequency and damping rate:
    t_work, ids_t = mix.get_array_oo(oo, t_filt, 't')
    line_t_wci = mix.test_array(t_work, 't[1/wci]', ':0.2e')
    erbar_s1_work = erbar_s1_filt[ids_t[0]:ids_t[-1]+1]
    del erbar_s1_filt, t_filt, ids_t

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
        print('--- ESTIMATION ' + line_t_wci + '---')
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
        print('--- ADVANCED FITTING ' + line_t_wci + '---')
        print(line_w_pr + '{:0.3e}'.format(w_adv))
        print(line_g_pr + '{:0.3e}'.format(g_adv))

    # results
    if flag_output:
        out = {'w-est': w_est, 'g-est': g_est,
               'w-adv': w_adv, 'g-adv': g_adv,
               'line_w': line_w, 'line_g': line_g}
        return out


def heat_vorticity(species_name, dd, oo={}):
    # species heat flux
    rd.radial_heat_flux(dd)
    efluxw_rad = dd[species_name].efluxw_rad
    t_wci_heat, ids_t_heat = mix.get_array_oo(
        oo, efluxw_rad['t'], 't')

    # vorticity
    vorbar(dd)
    t_wci_vorbar, ids_t_vorbar = mix.get_array_oo(
        oo, dd['vorbar']['t'], 't')

    # time evolution in different domains:
    ns_av = oo.get('ns_av', 1)
    heat_s_av = {}
    vorbar_s_av = {}
    lines_s_av_heat = {}
    lines_s_av_vorbar = {}
    for is_av in range(ns_av):
        # heat flux
        s_av, ids_s_av = mix.get_array_oo(oo, efluxw_rad['s'],
                                    's_av{:d}_heat'.format(is_av + 1))
        temp = mix.get_slice(
            efluxw_rad['data'], ids_t_heat, ids_s_av)
        heat_s_av[is_av] = np.mean(temp, axis=1)
        lines_s_av_heat[is_av] = mix.test_array(s_av, 's', ':0.3f')

        # vorticity
        s_av, ids_s_av = mix.get_array_oo(oo, dd['vorbar']['s'],
                                    's_av{:d}_vorbar'.format(is_av + 1))
        temp = mix.get_slice(
            dd['vorbar']['data'], ids_t_vorbar, ids_s_av)
        vorbar_s_av[is_av] = np.mean(temp, axis=1)
        lines_s_av_vorbar[is_av] = mix.test_array(s_av, 's', ':0.3f')
    del temp, s_av, ids_s_av

    # plotting
    for is_av in range(ns_av):
        curves = crv.Curves().xlab('t[wci^{-1}]') \
            .tit('<heat\ flux>_s\ vs\ |<\overline{\Omega}_r>_s|')
        curves.flag_norm = True
        curves.new('y{:d}'.format(is_av)) \
            .XS(t_wci_heat).YS(heat_s_av[is_av])\
            .leg('<heat\ flux>_s:\ ' + lines_s_av_heat[is_av])
        curves.new('y{:d}'.format(is_av)) \
            .XS(t_wci_vorbar).YS(np.abs(vorbar_s_av[is_av])) \
            .leg('|<\overline{\Omega}_r>_s|:\ '
                            + lines_s_av_vorbar[is_av]) \
            .sty('-.')
        cpr.plot_curves(curves)


def exp_AUG20787(dd, oo={}):
    # curves, where to add the data
    curves = oo.get('curves', crv.Curves())
    sel_norm = oo.get('sel_norm', 'wci')  # -> 'wci', 'khz', 'csa', 'csr'
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'

    # experimental GAM frequency:
    rho_fr = np.array([0.863, 0.879, 0.891, 0.902, 0.912, 0.912,
              0.922, 0.922, 0.932, 0.932, 0.941, 0.950,
              0.959, 0.967, 0.975, 0.984])
    rho_err = 0.005 * np.ones(np.size(rho_fr))
    fr_GAM_kHz = np.array([20.2, 20.1, 18.9, 18.9, 17.7, 15.5, 16.5,
                  14.0, 14.5, 13.4, 14.0, 13.4, 13.2, 12.2,
                  12.2, 12.5])  # kHz
    fr_err = 0.4 * np.ones(np.size(rho_fr))

    # experimental profiles:
    rho_T_AUG = np.array([0.863, 0.879, 0.891, 0.902, 0.912,
                 0.922, 0.932, 0.941, 0.950, 0.959,
                 0.967, 0.975, 0.984])
    Te_AUG = np.array([240,   210,   185,   165,   140,
              122,   105,    85,   70,     52,
               40,    25,    15])  # eV
    Ti_AUG = np.array([240,   210,   185,   165,   140,
              122,   105,    85,    70,    57,
               50,    38,    30])  # eV

    # choose normalization:
    coef_norm = None
    if sel_norm == 'khz':
        coef_norm = 1
    if sel_norm == 'wci':
        coef_norm = 1.e3 * 2 * np.pi / dd['wc']
    if sel_norm == 'csa':
        coef_norm = 1.e3 * 2 * np.pi / (dd['cs']/dd['a0'])
    if sel_norm == 'csr':
        coef_norm = 1.e3 * 2 * np.pi / (dd['cs']/dd['R0'])
    fr_res     = fr_GAM_kHz * coef_norm
    fr_err_res = fr_err     * coef_norm

    if sel_r == 's':
        r = np.sqrt(rho_fr)
        r_err = np.zeros([2, np.size(rho_fr)])
        r_err[0] = np.sqrt(rho_fr) - np.sqrt(rho_fr - rho_err)
        r_err[1] = np.sqrt(rho_fr + rho_err) - np.sqrt(rho_fr)
    if sel_r == 'psi':
        r = rho_fr
        r_err = rho_err

    curves.new('aug20787').XS(r).YS(fr_res)\
        .XS_ERR(r_err).YS_ERR(fr_err_res).leg('EXP.\ AUG20787')\
        .sty('o').col('red').ms(1)
    # curves.new('aug20787').XS(s_fr).YS(fr_res) \
    #     .leg('EXP.\ AUG20787') \
    #     .sty('o').col('red')

    return curves


def test_low_pass_filtering():
    L = 1
    x = np.linspace(0, L, 101)
    k1 = 2*np.pi / L
    k2 = 10 * 2*np.pi / L
    # y = np.cos(k1*x) + np.cos(k2*x)
    y = 2 + np.cos(k2 * x)

    order = 10
    # fr_threshold = k2 / (2*np.pi) / 3.
    fr_threshold = k2 / (2 * np.pi) / 1.
    sos = signal.butter(order, fr_threshold,
        'hp', fs=np.size(x), output='sos')
    y_filt = signal.sosfilt(sos, y)

    curves = crv.Curves().xlab('x').ylab('y')
    curves.new('y').XS(x).YS(y).leg('init')
    curves.new('filt').XS(x).YS(y_filt).leg('filt')
    cpr.plot_curves(curves)


def test_low_pass_filtering_fft():
    L = 1
    t = np.linspace(0, L, 369)
    k1 = 2*np.pi / L
    k2 = 4 * 2*np.pi / L
    k3 = 10 * 2*np.pi / L
    k4 = 14 * 2 * np.pi / L
    y = np.cos(k1*t) + np.cos(k2*t) + np.cos(k3*t) + np.cos(k4*t)
    # y = 2 + np.cos(k2 * t)

    # --- test two-sided FFT ---
    # w, w2, nw = ymath.w_for_fft(t)
    # y_work, _ = ymath.prepare_y_for_fft(t, y, nw)
    # f, _, f2_arranged = ymath.fft_y(y_work, nw, oo={'flag_f2_arranged': True})
    # curves = crv.Curves().xlab('w').ylab('yw')
    # curves.new('yw').XS(w2).YS(f2_arranged)
    # cpr.plot_curves(curves)

    # --- test filtering ---
    oo_filt = {'w_interval': [0., 2 * k1 / (2*np.pi)]}
    out = ymath.high_pass_filter(t, y, oo_filt)

    curves = crv.Curves().xlab('t').ylab('y')
    curves.new('y').XS(t).YS(y).leg('init')
    curves.new('filt').XS(t).YS(out['filt'])\
        .leg('filt').sty(':')
    cpr.plot_curves(curves)

    curves = crv.Curves().xlab('w').ylab('yw')
    curves.new('y').XS(out['w']).YS(out['fft_init']).leg('init')
    curves.new('filt').XS(out['w']).YS(out['fft_filt']) \
        .leg('filt').sty(':')
    cpr.plot_curves(curves)

    curves = crv.Curves().xlab('w').ylab('yw')
    curves.new('y').XS(out['w2']).YS(out['fft_init_2']).leg('init')
    curves.new('filt').XS(out['w2']).YS(out['fft_filt_2']) \
        .leg('filt').sty(':')
    cpr.plot_curves(curves)






















