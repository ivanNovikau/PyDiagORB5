import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np
from scipy import interpolate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)


# radial derivative of the full potential at a particular chi:
def ersc_old(dd, oo):
    chi_s = oo.get('chi_s', [0.0])

    rd.potsc_ref(dd, oo)
    t   = dd['potsc']['t']
    s   = dd['potsc']['s']
    chi = dd['potsc']['chi']

    nchi = np.size(chi_s)
    names_ersc_chi = []
    for count_chi in range(nchi):
        one_chi = chi_s[count_chi]
        name_ersc_chi = 'ersc-chi-' + '{:0.3f}'.format(one_chi)
        names_ersc_chi.append(name_ersc_chi)

        if name_ersc_chi in dd:
            continue

        id_chi, chi1 = mix.find(chi, one_chi)

        data = - np.gradient(dd['potsc']['data'][:, id_chi, :], s, axis=1)
        dd[name_ersc_chi] = {
            't': t,
            's': s,
            'chi_1': chi1, 'id_chi_1': id_chi,
            'data': data}
    return names_ersc_chi


# radial derivative of the full potential at a particular chi:
def ersc(dd, oo):
    chi_s = oo.get('chi_s', [0.0])

    names_potsc = rd.potsc_chi(dd, oo)

    nchi = np.size(chi_s)
    names_ersc_chi = []
    for count_chi in range(nchi):
        one_chi       = chi_s[count_chi]
        name_ersc_chi = 'ersc-chi-' + '{:0.3f}'.format(one_chi)
        names_ersc_chi.append(name_ersc_chi)
        if name_ersc_chi in dd:
            continue

        one_name_potsc = names_potsc[count_chi]
        data = - np.gradient(
            dd[one_name_potsc]['data'], dd['potsc_grids']['s'], axis=1)
        dd[name_ersc_chi] = {
            'chi_1':    dd[one_name_potsc]['chi_1'],
            'id_chi_1': dd[one_name_potsc]['id_chi_1'],
            'data': data}
    return names_ersc_chi


# phibar along potsc (t,s) grids:
def phibar_interp(dd):
    if 'phibar_interp' in dd:
        return

    rd.phibar(dd)
    rd.potsc_grids(dd)

    t = dd['phibar']['t']
    s = dd['phibar']['s']

    f_interp = interpolate.interp2d(t, s, dd['phibar']['data'].T)
    dd['phibar_interp'] = {
        'data': f_interp(dd['potsc_grids']['t'], dd['potsc_grids']['s']).T
    }


# non-zonal potential:
def phinz(dd, oo):
    chi_s = oo.get('chi_s', [0.0])

    phibar_interp(dd)
    names_potsc = rd.potsc_chi(dd, oo)

    nchi = np.size(chi_s)
    names_phinz_chi = []
    for count_chi in range(nchi):
        one_chi        = chi_s[count_chi]
        name_phinz_chi = 'phinz-chi-' + '{:0.3f}'.format(one_chi)
        names_phinz_chi.append(name_phinz_chi)
        if name_phinz_chi in dd:
            continue

        one_name_potsc = names_potsc[count_chi]
        data = dd[one_name_potsc]['data'] - dd['phibar_interp']['data']
        dd[name_phinz_chi] = {
            'chi_1':    dd[one_name_potsc]['chi_1'],
            'id_chi_1': dd[one_name_potsc]['id_chi_1'],
            'data': data}
    return names_phinz_chi


# non-zonal radial derivative of the full potential at a particular chi:
def ernz(dd, oo):
    chi_s = oo.get('chi_s', [0.0])

    names_phinz = phinz(dd, oo)

    nchi = np.size(chi_s)
    names_ernz_chi = []
    for count_chi in range(nchi):
        one_chi       = chi_s[count_chi]
        name_ernz_chi = 'ernz-chi-{:0.3f}'.format(one_chi)
        names_ernz_chi.append(name_ernz_chi)
        if name_ernz_chi in dd:
            continue

        one_name_phinz = names_phinz[count_chi]
        data = - np.gradient(
            dd[one_name_phinz]['data'], dd['potsc_grids']['s'], axis=1)
        dd[name_ernz_chi] = {
            'chi_1':    dd[one_name_phinz]['chi_1'],
            'id_chi_1': dd[one_name_phinz]['id_chi_1'],
            'data': data}
    return names_ernz_chi


def plot_t(dd, oo={}):
    sel_norm = oo.get('sel_norm', 'wci')
    s_intervals = oo.get('s_intervals', [[0.0, 1.0]])
    chi1 = oo.get('chi1', 0.0)

    rd.potsc(dd)
    t = dd['potsc']['t']
    s = dd['potsc']['s']
    chi = dd['potsc']['chi']

    # radial interval;
    ns_intervals = np.shape(s_intervals)[0]
    ids_s_intervals, s_final_intervals = [], []
    lines_s = []
    for count_s in range(ns_intervals):
        one_s, one_ids_s = mix.get_array(
                s, s_intervals[count_s][0], s_intervals[count_s][-1]
            )
        s_final_intervals.append(one_s)
        ids_s_intervals.append(one_ids_s)
        lines_s.append('s = [{:.3f}, {:.3f}]'.format(
            one_s[0], one_s[-1])
        )
    del one_s, one_ids_s, count_s

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

    # time interval
    t, ids_t = mix.get_array_oo(oo, t, 't')

    # poloidal angle
    id_chi, chi1 = mix.find(chi, chi1)
    line_chi1 = '\chi = {:0.1f}'.format(chi1)

    # plotting:
    line_Phi = '\Phi({:s})'.format(line_chi1)
    curves = crv.Curves().xlab(line_t).ylab(line_Phi)
    curves.flag_semilogy = True
    for count_s in range(ns_intervals):
        ids_s = ids_s_intervals[count_s]
        Phi = np.mean(
            dd['potsc']['data'][
                ids_t[0]:ids_t[-1]+1, id_chi, ids_s[0]:ids_s[-1]+1
            ], axis=1)
        curves.new(str(count_s)).XS(t).YS(Phi)\
            .leg(lines_s[count_s]).new_sty(count_s)
    cpr.plot_curves(curves)


def plot_schi(dd, t1, oo={}):
    rd.potsc(dd)
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    s, ids_s     = mix.get_array_oo(oo, s, 's')
    chi, ids_chi = mix.get_array_oo(oo, chi, 'chi')
    id_t1, _ = mix.find(t, t1)

    # signal in the chosen intervals
    pot_nz = mix.get_slice(dd['potsc']['data'], id_t1, ids_chi, ids_s)

    # form 3d curve:
    curves = crv.Curves().xlab('s').ylab('\chi').tit('\Phi')
    curves.new('Phi_schi').XS(s).YS(chi).ZS(pot_nz) \
        .leg('\Phi').cmp('hot')
    cpr.plot_curves_3d(curves, 'mat')


def plot_rz(dd, t1, oo={}):
    rd.potsc(dd)
    t = dd['potsc']['t']  # create a new reference
    r = dd['potsc']['r']  # create a new reference
    z = dd['potsc']['z']  # create a new reference

    # intervals
    id_t1, _ = mix.find(t, t1)

    # signal in the chosen intervals
    pot_nz = np.array(dd['potsc']['data'][id_t1, :, :])  # actually copy data

    # form 3d curve:
    curves = crv.Curves().xlab('r').ylab('z').tit('\Phi')
    curves.new('Phi_schi').XS(r).YS(z).ZS(pot_nz) \
        .leg('\Phi').cmp('jet')
    cpr.plot_curves_3d(curves, 'mat')


def calc_gamma_chi0(dd, oo):
    # initial signal and grids
    rd.potsc_grids(dd)
    t = dd['potsc_grids']['t']
    s = dd['potsc_grids']['s']
    chi = dd['potsc_grids']['chi']

    # parameters to treat the results
    sel_norm = oo.get('sel_norm', 'wci')
    s_intervals = oo.get('s_intervals', [[0.0, 1.0]])
    t_intervals = oo.get('t_intervals', [[t[0], t[-1]]])  # taking into account sel_norm
    filters = oo.get('filters', [None])
    chi1 = oo.get('chi1', 0.0)
    flag_est = oo.get('flag_est', True)
    # flag_adv = oo.get('flag_adv', True)

    # number of t- and s- intervals has to be the same:
    if np.shape(s_intervals)[0] != np.shape(t_intervals)[0]:
        return None
    if np.shape(s_intervals)[0] != np.shape(filters)[0]:
        return None

    # full potential at a chosen poloidal angle
    oo_phi = {'chi_s': [chi1]}
    name_potsc = rd.potsc_chi(dd, oo_phi)[0]

    # time normalization
    coef_norm = None
    line_t, line_w, line_g = None, None, None
    if sel_norm == 'ms':
        coef_norm = 1.e3 / dd['wc']
        line_t, line_w, line_g = 't,\ ms', 'w(kHz) = ', 'g(1e3/s) = '
    if sel_norm == 'wci':
        coef_norm = 1
        line_t, line_w, line_g = 't[\omega_c^{-1}]', 'w[wci] = ', 'g[wci] = '
    if sel_norm == 'csa':
        coef_norm = (dd['cs'] / dd['a0']) / dd['wc']
        line_t, line_w, line_g = 't[a_0/c_s]', 'w[cs/a0] = ', 'g[cs/a0 ] = '
    if sel_norm == 'csr':
        coef_norm = (dd['cs'] / dd['R0']) / dd['wc']
        line_t, line_w, line_g = 't[R_0/c_s]', 'w[cs/R0] = ', 'g[cs/R0] = '

    # form s and t intervals
    n_intervals = np.shape(s_intervals)[0]
    ids_s_intervals, _, lines_s = \
        mix.get_interval(s, s_intervals, 's', '0.3f')
    ids_t_intervals, t_final_intervals, lines_t = \
        mix.get_interval(t * coef_norm, t_intervals, 't', '0.1e')
    del s

    # poloidal angle
    id_chi, chi1 = mix.find(chi, chi1)
    line_chi1 = '\chi = {:0.1f}'.format(chi1)

    # signal description
    line_Phi = '\Phi({:s})'.format(line_chi1)

    # --- ESTIMATION ---
    if not flag_est:
        return

    wg_est, Phis_init, Phis_filt, Phis_fft_init, Phis_fft_filt, w2_grids = \
        {}, {}, {}, {}, {}, {}
    filt, Phi, t_int, w_int, ids_s, ids_t = None, None, None, None, None, None
    for id_int in range(n_intervals):
        ids_s, ids_t = ids_s_intervals[id_int], ids_t_intervals[id_int]
        t_int = t_final_intervals[id_int]
        oo_filter = filters[id_int]

        # averaging along s axis at a particular angle chi
        Phi = np.mean(
            dd[name_potsc]['data'][:, ids_s[0]:ids_s[-1] + 1], axis=1)

        # filtering
        filt = ymath.filtering(t, Phi, oo_filter)

        Phis_init[str(id_int)] = Phi[ids_t[0]:ids_t[-1] + 1]
        Phis_filt[str(id_int)] = filt['filt'][ids_t[0]:ids_t[-1] + 1]
        Phis_fft_init[str(id_int)] = filt['fft_init_2']
        Phis_fft_filt[str(id_int)] = filt['fft_filt_2']
        w2_grids[str(id_int)] = filt['w2']

        # estimation of the instability spectrum
        Phi_work = Phis_filt[str(id_int)]
        if Phi_work is None:
            Phi_work = Phi
        wg_est[str(id_int)] = ymath.estimate_wg(t_int, Phi_work)
    del filt, Phi, t_int, w_int, ids_s, ids_t

    # plot FFT
    for id_int in range(n_intervals):
        curves_est = crv.Curves().xlab(line_t).ylab('FFT:\ \Phi')\
            .tit('FFT:\ ' + line_Phi + ':\ ' + lines_s[id_int])
        curves_est.new('init') \
            .XS(w2_grids[str(id_int)]) \
            .YS(Phis_fft_init[str(id_int)]) \
            .leg('init').col('grey')
        curves_est.new('init') \
            .XS(w2_grids[str(id_int)]) \
            .YS(Phis_fft_filt[str(id_int)]) \
            .leg('filt').col('blue').sty(':')
        cpr.plot_curves(curves_est)

    # plot time evolution
    for id_int in range(n_intervals):
        curves_est = crv.Curves().xlab(line_t).ylab('\Phi')\
            .tit(line_Phi + ':\ ' + lines_s[id_int])
        curves_est.flag_semilogy = True
        curves_est.new('init')\
            .XS(t_final_intervals[id_int])\
            .YS(Phis_init[str(id_int)])\
            .leg('init').col('grey')
        curves_est.new('init') \
            .XS(t_final_intervals[id_int]) \
            .YS(Phis_filt[str(id_int)]) \
            .leg('filt').col('blue').sty(':')
        curves_est.new('peaks')\
            .XS(wg_est[str(id_int)]['x_peaks'])\
            .YS(wg_est[str(id_int)]['y_peaks'])\
            .leg('peaks').sty('o').col('green')
        curves_est.new('fitting')\
            .XS(wg_est[str(id_int)]['x_fit'])\
            .YS(wg_est[str(id_int)]['y_fit'])\
            .leg('fitting').col('red').sty('--')
        cpr.plot_curves(curves_est)

    print('--- Estimation ---')
    for id_int in range(n_intervals):
        print('E -> *** ' + lines_s[id_int] + ':\ ' + lines_t[id_int] + ' ***')
        print('E -> ' + line_w + '{:0.3e}'.format(wg_est[str(id_int)]['w']))
        print('E -> ' + line_g + '{:0.3e}'.format(wg_est[str(id_int)]['g']))


def calc_wg(dd, oo):
    # initial signal and grids
    rd.potsc_grids(dd)
    t = dd['potsc_grids']['t']
    s = dd['potsc_grids']['s']
    chi = dd['potsc_grids']['chi']

    # parameters to treat the results
    sel_norm = oo.get('sel_norm', 'wc')
    oo_filter_g = oo.get('filter_g', {'sel_filt': None})
    oo_filter_w = oo.get('filter_w', {'sel_filt': None})
    chi1 = oo.get('chi1', 0.0)

    oo_w_fft = oo.get('w_fft', None)
    n_max_w_fft = None
    if oo_w_fft is None:
        flag_w_fft = False
    else:
        flag_w_fft = True
        n_max_w_fft = oo_w_fft.get('n_max', 1)

    # full potential at a chosen poloidal angle
    oo_phi = {'chi_s': [chi1]}
    name_potsc = rd.potsc_chi(dd, oo_phi)[0]

    # time normalization
    coef_norm_t, coef_norm_w = None, None
    line_norm_t, line_norm_w, line_norm_g = None, None, None
    label_norm_t, label_norm_w, label_norm_g = None, None, None
    if sel_norm == 'khz':
        coef_norm_t = 1.e3 / dd['wc']
        coef_norm_w = 1
        line_norm_t,  line_norm_w,  line_norm_g  = 't(ms)', 'w(kHz)', 'g(1e3/s)'
        label_norm_t, label_norm_w, label_norm_g = \
            't, ms', '\omega(kHz)', '\gamma(1e3/s)'
    if sel_norm == 'wc':
        coef_norm_t = 1
        coef_norm_w = 2 * np.pi
        line_norm_t, line_norm_w,   line_norm_g  = 't[1/wc]', 'w[wci]', 'g[wci]'
        label_norm_t, label_norm_w, label_norm_g = \
            't[\omega_c^{-1}]', '\omega[\omega_c]', '\gamma[\omega]'
    if sel_norm == 'csa':
        coef_norm_t = (dd['cs'] / dd['a0']) / dd['wc']
        coef_norm_w = 2 * np.pi
        line_norm_t, line_norm_w, line_norm_g = 't[a0/cs]', 'w[cs/a0]', 'g[cs/a0 ]'
        label_norm_t, label_norm_w, label_norm_g = \
            't[a_0/c_s]', '\omega[c_a/a_0]', '\gamma[c_s/a_0]'
    if sel_norm == 'csr':
        coef_norm_t = (dd['cs'] / dd['R0']) / dd['wc']
        coef_norm_w = 2 * np.pi
        line_norm_t, line_norm_w, line_norm_g = 't[R0/cs]', 'w[cs/R0]', 'g[cs/R0]'
        label_norm_t, label_norm_w, label_norm_g = \
            't[R_0/c_s]', '\omega[c_a/R_0]', '\gamma[c_s/R_0]'

    oo_filter_g['norm_w'] = coef_norm_w
    oo_filter_w['norm_w'] = coef_norm_w

    # radial domain
    s, ids_s = mix.get_array_oo(oo, s, 's')
    line_s = 's = [{:0.3f}, {:0.3f}]'.format(s[0], s[-1])

    # time renormalization
    t_full = t * coef_norm_t
    del t

    # poloidal angle
    id_chi, chi1 = mix.find(chi, chi1)
    line_chi1 = '\chi = {:0.1f}'.format(chi1)

    # signal description
    line_Phi = '\Phi({:s})'.format(line_chi1)
    line_Phi_w = '\Phi({:s}) \cdot \exp(-g t)'.format(line_chi1)

    # averaging along s axis at a particular angle chi
    Phi_init_g = np.mean(dd[name_potsc]['data'][:, ids_s[0]:ids_s[-1] + 1], axis=1)

    # filtering of the growing signal
    filt_g = ymath.filtering(t_full, Phi_init_g, oo_filter_g)

    # information about the time domain, where the filtering of the growing signal
    # has been performed
    line_filt_tg = label_norm_t + ' = [{:0.3e}, {:0.3e}]' \
        .format(filt_g['x'][0], filt_g['x'][-1])

    # chose a time interval inside of the time domain
    # where the filtering has been performed
    tg, ids_tg = mix.get_array_oo(oo, filt_g['x'], 't')

    # information about the time domain where the growth rate is estimated:
    line_tg = line_norm_t + ' = [{:0.3e}, {:0.3e}]'.format(tg[0], tg[-1])
    if filt_g['filt'] is None:
        Phi_filt_g = None
        Phi_work_g = Phi_init_g[ids_tg[0]:ids_tg[-1] + 1]
    else:
        Phi_filt_g = filt_g['filt']
        Phi_work_g = Phi_filt_g[ids_tg[0]:ids_tg[-1] + 1]

    # estimation of the instability growth rate
    g_est = ymath.estimate_g(tg, Phi_work_g)

    # get rid of the growth:
    Phi_init_w = Phi_init_g * np.exp(- g_est['g'] * t_full)

    # filtering of the oscillating signal
    filt_w = ymath.filtering(t_full, Phi_init_w, oo_filter_w)

    # information about the time domain, where the filtering of the oscillating signal
    # has been performed
    line_filt_tw = label_norm_t + ' = [{:0.3e}, {:0.3e}]' \
        .format(filt_w['x'][0], filt_w['x'][-1])

    # chose a time interval inside of the time domain
    # where the filtering has been performed
    tw, ids_tw = mix.get_array_oo(oo, filt_w['x'], 't')

    # information about the time domain where the frequency is estimated:
    line_tw = line_norm_t + ' = [{:0.3e}, {:0.3e}]'.format(tw[0], tw[-1])
    if filt_w['filt'] is None:
        Phi_filt_w = None
        Phi_work_w = Phi_init_w[ids_tw[0]:ids_tw[-1] + 1]
    else:
        Phi_filt_w = filt_w['filt']
        Phi_work_w = Phi_filt_w[ids_tw[0]:ids_tw[-1] + 1]

    # estimation of the frequency:
    if flag_w_fft:
        w_est_fft = ymath.estimate_w_max_fft(filt_w['w'],
                                             filt_w['fft_filt'], {'n_max': n_max_w_fft})
    w_est = ymath.estimate_w(tw, Phi_work_w)

    # plot FFT
    curves = crv.Curves().xlab(label_norm_w).ylab('FFT:\ \Phi')\
        .tit('FFT:\ ' + line_Phi)\
        .titn(line_s + ', ' + line_filt_tg)
    curves.new() \
        .XS(filt_g['w2']) \
        .YS(filt_g['fft_init_2']) \
        .leg('initial').col('grey')
    curves.new() \
        .XS(filt_g['w2']) \
        .YS(filt_g['fft_filt_2']) \
        .leg('filtered').col('blue').sty(':')
    cpr.plot_curves(curves)

    curves = crv.Curves().xlab(label_norm_w).ylab('FFT:\ \Phi') \
        .tit('FFT:\ ' + line_Phi_w)\
        .titn(line_s + ', ' + line_filt_tw)
    curves.new() \
        .XS(filt_w['w2']) \
        .YS(filt_w['fft_init_2']) \
        .leg('initial').col('grey')
    curves.new() \
        .XS(filt_w['w2']) \
        .YS(filt_w['fft_filt_2']) \
        .leg('filtered').col('blue').sty(':')
    cpr.plot_curves(curves)

    if flag_w_fft:
        curves = crv.Curves().xlab(label_norm_w).ylab('FFT:\ \Phi') \
            .tit('FFT:\ ' + line_Phi_w) \
            .titn(line_s + ', ' + line_filt_tw)
        curves.new() \
            .XS(filt_w['w']) \
            .YS(filt_w['fft_filt']) \
            .leg('filtered').col('blue').sty('-')
        curves.new() \
            .XS(w_est_fft['w_max']) \
            .YS(w_est_fft['f_max']) \
            .leg('chosen\ fft\ peaks').col('orange').sty('o')
        cpr.plot_curves(curves)

    # plot filtered signals
    curves = crv.Curves().xlab(label_norm_t).ylab(line_Phi)\
        .tit(line_Phi + '\ (-> G)\ VS\ ' + line_Phi_w + '\ (-> W)')\
        .titn(line_s)
    curves.flag_semilogy = True
    curves.new() \
        .XS(t_full) \
        .YS(Phi_init_g) \
        .leg('G:\ init.').col('grey')
    curves.new() \
        .XS(filt_g['x']) \
        .YS(Phi_filt_g) \
        .leg('G:\ filt.').col('blue').sty(':')
    curves.new() \
        .XS(t_full) \
        .YS(Phi_init_w) \
        .leg('W:\ init.').col('red')
    curves.new() \
        .XS(filt_w['x']) \
        .YS(Phi_filt_w) \
        .leg('W:\ filt.').col('green').sty(':')
    cpr.plot_curves(curves)

    # plot fitted signals
    curves = crv.Curves().xlab(label_norm_t).ylab(line_Phi) \
        .tit(line_Phi + ':\ ' + line_s)\
        .titn(label_norm_g + ' = {:0.3e}'.format(g_est['g']))
    curves.flag_semilogy = True
    curves.new() \
        .XS(tg) \
        .YS(Phi_work_g) \
        .leg('signal').col('blue')
    curves.new() \
        .XS(g_est['x_peaks']) \
        .YS(g_est['y_peaks']) \
        .leg('peaks').sty('o').col('green')
    curves.new() \
        .XS(g_est['x_fit']) \
        .YS(g_est['y_fit']) \
        .leg('fitting').col('red').sty('--')
    cpr.plot_curves(curves)

    if w_est is not None:
        curves = crv.Curves().xlab(label_norm_t).ylab(line_Phi) \
            .tit(line_Phi_w + ':\ ' + line_s) \
            .titn(label_norm_w + ' = {:0.3e}'.format(w_est['w']))
    else:
        curves = crv.Curves().xlab(label_norm_t).ylab(line_Phi) \
            .tit(line_Phi_w + ':\ ' + line_s)
    curves.flag_semilogy = True
    curves.new() \
        .XS(tw) \
        .YS(Phi_work_w) \
        .leg('signal').col('blue')
    if w_est is not None:
        curves.new('peaks') \
            .XS(w_est['x_peaks']) \
            .YS(w_est['y_peaks']) \
            .leg('peaks').sty('o').col('green')
        curves.new('fitting') \
            .XS(w_est['x_fit']) \
            .YS(w_est['y_fit']) \
            .leg('fitting').col('red').sty('--')
    cpr.plot_curves(curves)

    print('--- Estimation ---')
    print('*** Growth rate: ' + line_s + ':\ ' + line_tg + ' ***')
    print('E -> ' + line_norm_g + ' = {:0.3e}'.format(g_est['g']))
    if w_est is not None:
        print('*** Frequency: '   + line_s + ':\ ' + line_tw + ' ***')
        print('E -> ' + line_norm_w + ' = {:0.3e}'.format(w_est['w']))
    if flag_w_fft:
        for i_w_fft in range(n_max_w_fft):
            print(
                    'w_fft(id_max = {:d}) = {:0.3e},   max(id_max = {:d}) = {:0.3e}'
                        .format(i_w_fft, w_est_fft['w_max'][i_w_fft],
                                i_w_fft, w_est_fft['f_max'][i_w_fft])
                  )


def calc_gamma_chimax(dd, s1, s2, oo={}):
    rd.potsc(dd)
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    _, ids_s = mix.get_array(s, s1, s2)
    t, ids_t = mix.get_array_oo(oo, t, 't')
    chi, ids_chi = mix.get_array_oo(oo, chi, 'chi')

    # non-zonal Phi in chosen intervals
    pot_chi = mix.get_slice(dd['potsc']['data'], ids_t, ids_chi, ids_s)

    # averaging on s
    pot_aver_s = np.mean(pot_chi, axis=2)

    # maximums along chi:
    pot_max = np.amax(np.abs(pot_aver_s), axis=1)

    # gamma calculation
    wg_est = ymath.estimate_g(t, pot_max)
    wg_adv = ymath.advanced_g(t, pot_max, {'est': wg_est})

    # plotting: estimation: time evolution and peaks
    curves = crv.Curves().xlab('t').ylab('\Phi')\
        .tit('Estimation:\                g = {:0.3e}'.format(wg_est['g']))\
        .tit('\\\\Full\ signal\ fitting:\ g = {:0.3e}'.format(wg_adv['g']))
    curves.flag_semilogy = True
    curves.new('max') \
        .XS(t).YS(pot_max).leg('max.\ along\ \chi')
    if 'x_peaks' in wg_est:
        curves.new('peaks') \
            .XS(wg_est['x_peaks']).YS(wg_est['y_peaks']) \
            .leg('peaks').sty('o').col('grey')
    curves.new('fitting') \
        .XS(t).YS(wg_est['y_fit']) \
        .leg('estimation').col('red').sty('--')
    curves.new('fitting') \
        .XS(wg_adv['x_fit']).YS(wg_adv['y_fit']) \
        .leg('full\ signal\ fitting').col('green').sty(':')
    cpr.plot_curves(curves)

    # print results (keep it to have from screen)
    print('--- Estimation ---')
    print('g = {:0.3e}'.format(wg_est['g']))
    print('--- Advanced ---')
    print('g = {:0.3e}'.format(wg_adv['g']))


def plot_schi_max_along_chi(dd, t1, oo={}):
    # plot max of phi on (s, chi) at a particular time step
    rd.potsc(dd)
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    s, ids_s     = mix.get_array_oo(oo, s, 's')
    chi, ids_chi = mix.get_array_oo(oo, chi, 'chi')
    id_t1, _ = mix.find(t, t1)

    # non-zonal Phi in chosen intervals
    pot_chi = mix.get_slice(dd['potsc']['data'], id_t1, ids_chi, ids_s)

    # maximums along chi:
    pot_max = np.amax(np.abs(pot_chi), axis=0)

    # build array for max:
    chi_max = np.zeros(np.size(s))
    for i in range(np.size(pot_max)):
        max1 = pot_max[i]
        pot_s1 = pot_chi[:, i]
        id_chi1 = np.where(pot_s1 == max1)
        if len(id_chi1[0]) is 0:
            id_chi1 = np.where(pot_s1 == -max1)
        id_chi1 = id_chi1[0][0]
        chi_max[i] = chi[id_chi1]

    # plotting:
    curves = crv.Curves().xlab('s').ylab('\chi').tit('\Phi')
    curves.new('schi') \
        .XS(s).YS(chi).ZS(pot_chi)\
        .leg('\Phi').cmp('hot')
    curves.new('max').XS(s).YS(chi_max).leg('max')
    cpr.plot_curves_3d(curves)


def plot_rz_max_along_chi(dd, t1, oo={}):
    rd.potsc(dd)
    t = dd['potsc']['t']  # create a new reference
    r = dd['potsc']['r']  # create a new reference
    z = dd['potsc']['z']  # create a new reference

    # intervals
    id_t1, _ = mix.find(t, t1)

    # signal in the chosen intervals
    pot_nz = np.array(dd['potsc']['data'][id_t1, :, :])  # actually copy data

    # maximums along chi:
    pot_max = np.amax(np.abs(pot_nz), axis=0)

    # build array for max:
    # r_max = np.zeros(np.shape(pot_nz))
    # z_max = np.zeros(np.shape(pot_nz))
    r_max = np.zeros(np.size(pot_max))
    z_max = np.zeros(np.size(pot_max))
    for id_j in range(np.size(pot_max)):
        max1 = pot_max[id_j]
        pot_s1 = pot_nz[:, id_j]
        id_i = np.where(pot_s1 == max1)
        if len(id_i[0]) is 0:
            id_i = np.where(pot_s1 == -max1)
        id_i = id_i[0][0]
        r_max[id_j] = r[id_i, id_j]
        z_max[id_j] = z[id_i, id_j]

    # plotting:
    curves = crv.Curves().xlab('r').ylab('z').tit('\Phi')
    curves.new('rz') \
        .XS(r).YS(z).ZS(pot_nz)\
        .leg('\Phi').cmp('hot')
    curves.new('max').XS(r_max).YS(z_max).leg('max')
    cpr.plot_curves_3d(curves)


def anim_st_chi0(dd, chi0, oo={}):
    rd.potsc(dd)
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    s, ids_s = mix.get_array_oo(oo, s, 's')
    t, ids_t = mix.get_array_oo(oo, t, 't')
    id_chi1, _ = mix.find(chi, chi0)

    # signal in the chosen intervals
    pot_nz = mix.get_slice(dd['potsc']['data'], ids_t, id_chi1, ids_s)

    # form 2d curve:
    curves = crv.Curves().xlab('s').wlab('t').zlab('\Phi').tit('\Phi')
    curves.new('anim_Phi_st').XS(s).WS(t).ZS(pot_nz).leg('\Phi')
    curves.set_limits()
    curves.flag_norm = True

    # animation:
    cpr.animation_curves_2d(curves)





