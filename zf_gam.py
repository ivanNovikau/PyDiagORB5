import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import gam_theory
import numpy as np
import copy
from scipy import interpolate


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


def fft_gam_2d(dd, oo={}):
    sel_norm = oo.get('sel_norm', 'wci')  # -> 'wci', 'kHz', 'csa', 'csr'
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'
    sel_cmp = oo.get('sel_cmp', 'pink-r')  # -> 'jet', 'hot' etc.

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

    # analytical and experimental data:
    oo_th = {'curves': curves_er, 'sel_norm': sel_norm,
             'sel_r': sel_r, 's': r}
    # curves_er = exp_AUG20787(dd, oo_th)
    curves_er = gam_theory.get_gao(dd, oo_th)
    # curves_er = gam_theory.get_gk_fit(dd, oo_th)

    oo_th.update({'curves': curves_vor})
    # curves_vor = exp_AUG20787(dd, oo_th)
    curves_vor = gam_theory.get_gao(dd, oo_th)
    # curves_vor = gam_theory.get_gk_fit(dd, oo_th)

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
        id_s = int(ids_t[count_s])
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
    sel_norm = oo.get('sel_norm', 'wci')  # -> 'wci', 'kHz', 'csa', 'csr'

    rd.phibar(dd)
    erbar(dd)
    t = dd['phibar']['t']  # normalized to wc
    s = dd['phibar']['s']

    # time interval
    t, ids_t = mix.get_array_oo(oo, t, 't')
    line_t_wci = mix.test_array(t, 't[wci^{-1}]', ':0.2e')

    # radial interval
    s, ids_r = mix.get_array_oo(oo, s, 's')

    # zonal radial electric field
    erbar_st = mix.get_slice(dd['erbar']['data'], ids_t, ids_r)

    ns = np.size(s)
    ws, gs = np.zeros(ns), np.zeros(ns)
    for id_s in range(ns):
        out = ymath.estimate_wg(t, erbar_st[:, id_s])
        if out is None:
            ws[id_s] = np.NaN
            gs[id_s] = np.NaN
        else:
            ws[id_s] = out['w']  # ~ 2pi / wc
            gs[id_s] = out['g']  # ~ 1/wc

    # normalization:
    if sel_norm == 'khz':
        coef_norm_w = 1. / (2 * np.pi * 1.e3)
        coef_norm_g = 1. / 1.e3
        line_w = '\omega,\ kHz'
        line_g = '\gamma,\ 1e3/s'
    if sel_norm == 'wci':
        coef_norm_w = coef_norm_g = 1 / dd['wc']
        line_w = '\omega[\omega_c]'
        line_g = '\gamma[\omega_c]'
    if sel_norm == 'csa':
        coef_norm_w = coef_norm_g = 1 / (dd['cs'] / dd['a0'])
        line_w = '\omega[c_s/a_0]'
        line_g = '\gamma[c_s/a_0]'
    if sel_norm == 'csr':
        coef_norm_w = coef_norm_g = 1 / (dd['cs'] / dd['R0'])
        line_w = '\omega[c_s/R_0]'
        line_g = '\gamma[c_s/R_0]'
    ws = ws * coef_norm_w
    gs = gs * coef_norm_g

    # plotting frequencies
    curves_w = crv.Curves().xlab('s').ylab(line_w)\
        .tit('\omega\ of\ \overline E_r:\ ' + line_t_wci)
    curves_w.new('f').XS(s).YS(ws)
    cpr.plot_curves(curves_w)

    # plotting dynamic rate
    curves_g = crv.Curves().xlab('s').ylab(line_g) \
        .tit('\gamma\ of\ \overline E_r:\ ' + line_t_wci)
    curves_g.new('f').XS(s).YS(gs)
    cpr.plot_curves(curves_g)


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



