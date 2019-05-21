import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import zf_gam as zf
import general as gn
import write_data as wr
import transport
from scipy.fftpack import next_fast_len
from polycoherence import _plot_signal, polycoherence, plot_polycoherence
import numpy as np
from scipy import interpolate
import h5py as h5
from scipy.signal import correlate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(zf)
    mix.reload_module(gn)
    mix.reload_module(wr)
    mix.reload_module(transport)


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


# phibar along potsc grids:
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


# non-zonal potential at some poloidal angles
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


# non-zonal potential at a time point
def phinz_t(dd, oo):
    names_potsc = rd.potsc_t(dd, oo)
    phibar_interp(dd)

    t_points = oo.get('t_points', [0.0])
    nt_points = np.size(t_points)

    names = []
    for count_t in range(nt_points):
        one_t = t_points[count_t]
        var_name = 'phinz-t-' + '{:0.3e}'.format(one_t)
        names.append(var_name)
        if var_name in dd:
            continue

        data = {}
        one_name_potsc = names_potsc[count_t]
        data['id_t_point'], data['t_point'] = mix.find(dd['potsc_grids']['t'], one_t)
        phibar_t1   = dd['phibar_interp']['data'][data['id_t_point'], :]
        data['data'] = dd[one_name_potsc]['data'] - phibar_t1[None, :]
        dd[var_name] = data
    return names


#  radial derivative of the non-zonal potential at a particular chi:
def ernz_r(dd, oo):
    chi_s = oo.get('chi_s', [0.0])

    names_phinz = phinz(dd, oo)

    nchi = np.size(chi_s)
    names_ernz_chi = []

    for count_chi in range(nchi):
        one_chi       = chi_s[count_chi]
        name_ernz_chi = 'ernz_r-chi-{:0.3f}'.format(one_chi)
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


# poloidal derivative of the non-zonal potential at a particular chi:
def ernz_chi(dd, oo):
    phibar_interp(dd)

    chi_s = oo.get('chi_s', [0.0])
    nchi_points = np.size(chi_s)
    t, s, chi = dd['potsc_grids']['t'], dd['potsc_grids']['s'], dd['potsc_grids']['chi']
    nt, ns, nchi = np.size(t), np.size(s), np.size(chi)

    # number of radial points to read at once
    max_allowed_size = dd['max_size_Gb'] * 1024 ** 3  # in bytes
    float_size = 8  # in bytes
    id_s_step = int( round( max_allowed_size / (float_size * nt * nchi) ) )

    names = []
    for count_chi in range(nchi_points):
        one_chi       = chi_s[count_chi]
        var_name = 'ernz_chi-chi-{:0.3f}'.format(one_chi)
        names.append(var_name)
        if var_name in dd:
            continue

        data = rd.read_signal(dd['path_ext'], var_name)
        if data is None:
            data = {}

            # --- calculate data ---
            f = h5.File(dd['path_orb'], 'r')
            data['data'] = np.zeros([ns, nt])  # transposed w.r.t final matrix
            data['id_chi_1'], data['chi_1'] = mix.find(chi, one_chi)
            line_read = '/data/var2d/generic/potsc/data'
            id_s_end = 0
            while id_s_end < ns:
                # read Phi at several points along s: Phi(t,chi,ids_current)
                id_s_begin = id_s_end
                id_s_end  += id_s_step
                if id_s_end >= ns:
                    id_s_end = ns
                ids_current = [i for i in range(id_s_begin, id_s_end)]
                potsc_s1 = np.array(f[line_read][:, :, ids_current])

                # find non-zonal Phi(t,chi,ids_current)
                phibar_s1 = dd['phibar_interp']['data'][:, ids_current]
                phinz_s1 = potsc_s1 - phibar_s1[:, None, :]
                del potsc_s1

                # chi-derivative of the non-zonal Phi(t,chi,ids_current)
                loc = np.gradient(phinz_s1, chi, axis=1)
                del phinz_s1

                # save non-zonal poloidal electric field E_chi(t,chi1,ids_current)
                count_id_s = -1
                for id_s_loc in ids_current:
                    count_id_s += 1
                    data['data'][id_s_loc] = \
                        - loc[:, data['id_chi_1'], count_id_s] / s[id_s_loc]
                del loc
            f.close()

            # for the sake of generality, transpose data (ns,nt) -> (nt,ns)
            data['data'] = data['data'].T

            # save the data to an external file
            desc = 'non-zonal poloidal electric field at chi = {:0.3f}'.format(data['chi_1'])
            wr.save_data_adv(dd['path_ext'], data, {'name': var_name, 'desc': desc})

        # --- save data to the structure ---
        dd[var_name] = data
    return names


def choose_var(dd, oo):
    chi_point = oo.get('chi_point', 0.0)
    oo_var = {'chi_s': [chi_point]}
    opt_var = oo.get('opt_var', 'ernz_r')
    var_name, tit_var = '', ''
    if opt_var == 'phinz':
        var_name = phinz(dd, oo_var)[0]
        tit_var = '(\Phi - \overline{\Phi}):\ '
    if opt_var == 'ernz_r':
        var_name = ernz_r(dd, oo_var)[0]
        tit_var = '-\partial_s(\Phi - \overline{\Phi}):\ '
    if opt_var == 'ernz_chi':
        var_name = ernz_chi(dd, oo_var)[0]
        tit_var = '-s^{-1}\partial_{\chi}(\Phi - \overline{\Phi}):\ '
    vvar = dd[var_name]['data']
    s = dd['potsc_grids']['s']
    t = dd['potsc_grids']['t']
    line_chi = '\chi = {:0.3f}'.format(dd[var_name]['chi_1'])
    tit_var = tit_var + line_chi

    res = {
        'var': vvar,
        's': s,
        't': t,
        'tit': tit_var
    }

    return res


def choose_var_t(dd, oo):
    t_point = oo.get('t_point', 0.0)
    oo_var = {'t_points': [t_point]}
    opt_var = oo.get('opt_var', 'potsc_t')
    var_name, tit_var = '', ''
    if opt_var == 'potsc':
        var_name = rd.potsc_t(dd, oo_var)[0]
        tit_var = '\Phi:\ '
    if opt_var == 'phinz':
        var_name = phinz_t(dd, oo_var)[0]
        tit_var = '(\Phi - \overline{\Phi}):\ '
    vvar = dd[var_name]['data']
    s   = dd['potsc_grids']['s']
    chi = dd['potsc_grids']['chi']
    line_chi = 't[\omega_{ci}^{-1}]' + ' = {:0.3e}'.format(dd[var_name]['t_point'])
    tit_var = tit_var + line_chi

    res = {
        'var': vvar,
        's': s,
        'chi': chi,
        'tit': tit_var
    }

    return res


def plot_st(dd, oo):
    out = choose_var(dd, oo)

    oo_st = dict(oo)
    oo_st.update(out)
    gn.plot_st(dd, oo_st)


def plot_aver_st(dd, oo):
    out = choose_var(dd, oo)
    vvar = out['var']
    s = out['s']
    t = out['t']
    tit_var = out['tit']

    # # --- averaging in time ---
    # oo_avt = dict(oo)
    # oo_avt.update({
    #     'vars': [vvar, vvar],
    #     'ts': [t, t], 'ss': [s, s],
    #     'opts_av': ['rms', 'mean'],
    #     'tit': tit_var, 'vars_names': ['', '']
    # })
    # gn.plot_avt(dd, oo_avt)

    # --- averaging in space ---
    oo_avs = oo
    oo_avs.update({
        'vars': [vvar], 'ts': [t], 'ss': [s],
        'opts_av': ['rms'],
        'tit': tit_var, 'vars_names': ['']
    })
    gn.plot_avs(dd, oo_avs)


def plot_fft(dd, oo):
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'

    out = choose_var(dd, oo)
    vvar = out['var']
    r = out['s']
    t = out['t']
    tit_var = out['tit']

    # radial coordinate normalization
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
    vvar = out['var']
    s = out['s']
    t = out['t']
    tit_var = out['tit']

    # plotting
    oo_fft = dict(oo)
    oo_fft.update({
        'vars': [vvar], 'ts_wci': [t], 'ss': [s],
        'tit': tit_var,
        'vars_names': ['']
    })
    gn.plot_fft_1d(dd, oo_fft)


def compare_nz_var_aver_st(dd, oo):
    # compare NZ signal with some other signal

    # --- non-zonal signal ---
    oo_nz = dict(oo)
    oo_nz['opt_var'] = oo['opt_var_nz']
    out = choose_var(dd, oo_nz)
    varnz, s_nz, t_nz, tit_varnz = out['var'], out['s'], out['t'], out['tit']

    # --- signal to compare with ---
    var_type = oo.get('var_type', 'zonal')
    oo_var = dict(oo)
    oo_var['opt_var'] = oo['opt_var']

    if var_type == 'zonal':
        out = zf.choose_var(dd, oo_var)
    if var_type == 'transport':
        out = transport.choose_var(dd, oo_var)
    if var_type == 'nonzonal':
        out = choose_var(dd, oo_var)
    vvar, s, t, tit_var = out['var'], out['s'], out['t'], out['tit']

    # preliminary parameters
    oo_av = dict(oo)
    oo_av.update({
        'vars': [vvar, varnz],
        'ts': [t, t_nz],
        'ss': [s, s_nz],
        'vars_names': [tit_var, tit_varnz]
    })

    # --- averaging in time ---
    oo_av.update({
        'opts_av': ['rms', 'rms'],
        'tit': 'averaging\ in\ time',
    })
    gn.plot_avt(dd, oo_av)

    # --- averaging in space ---
    oo_avs = dict(oo)
    oo_avs.update({
        'opts_av': ['mean', 'mean'],
        'tit': 'averaging\ in\ space',
    })
    gn.plot_avs(dd, oo_av)


def compare_nz_var_fft_1d(dd, oo):
    # compare NZ signal with some other signal

    # --- non-zonal signal ---
    oo_nz = dict(oo)
    oo_nz['opt_var'] = oo['opt_var_nz']
    out = choose_var(dd, oo_nz)
    varnz, s_nz, t_nz, tit_varnz = out['var'], out['s'], out['t'], out['tit']

    # --- signal to compare with ---
    var_type = oo.get('var_type', 'zonal')
    oo_var = dict(oo)
    oo_var['opt_var'] = oo['opt_var']

    if var_type == 'zonal':
        out = zf.choose_var(dd, oo_var)
    if var_type == 'transport':
        out = transport.choose_var(dd, oo_var)
    if var_type == 'nonzonal':
        out = choose_var(dd, oo_var)
    vvar, s, t, tit_var = out['var'], out['s'], out['t'], out['tit']

    # plotting
    oo_fft = dict(oo)
    oo_fft.update({
        'vars': [varnz, vvar],
        'ts_wci': [t_nz, t],
        'ss': [s_nz, s],
        'tit': 'FFT',
        'vars_names': [tit_varnz, tit_var]
    })
    gn.plot_fft_1d(dd, oo_fft)

    return


def find_gamma(dd, oo):
    out = choose_var(dd, oo)
    vvar, s, t, tit_var = out['var'], out['s'], out['t'], out['tit']

    oo_gamma = oo  # oo must have some 'opt_av'
    oo_gamma.update({
        'var': vvar, 't': t, 's': s,
        'tit': tit_var, 'var_name': ''
    })
    gn.find_gamma(dd, oo_gamma)

    return


def plot_t(dd, oo):
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


def plot_schi(dd, t1, oo):
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


def plot_rz(dd, oo):
    res = choose_var_t(dd, oo)
    r = dd['potsc_grids']['r']  # create a new reference
    z = dd['potsc_grids']['z']  # create a new reference

    # form 3d curve:
    curves = crv.Curves().xlab('r').ylab('z').tit(res['tit'])
    curves.new('Phi_schi')\
        .XS(r)\
        .YS(z)\
        .ZS(res['var'].T)\
        .cmp('jet')
    cpr.plot_curves_3d(curves)


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


def test_correlation(dd):
    Nt = 401
    t_min = 1.0
    t_max = 4.3
    t = np.linspace(t_min, t_max, Nt)

    delta_t = 0.05
    T1, T2 = 0.25, 0.25
    w1, w2 = 2*np.pi/T1, 2*np.pi/T2
    y1, y2 = np.sin(w1 * t), np.sin(w1 * (t + delta_t))

    oo_dt = {
        'var1': y1, 'var2': y2,
        'grid_t1': t, 'grid_t2': t,
        'vars_names': ['y1', 'y2']
    }
    gn.find_time_delay(dd, oo_dt)


def test1_bicoherence():
    N = 5001
    t = np.linspace(0, 100, N)
    fs = 1 / (t[1] - t[0])
    s1 = np.cos(2 * np.pi * 4 * t + 0.2)
    s2 = 3 * np.cos(2 * np.pi * 5 * t + 0.5)
    np.random.seed(0)
    noise = 5 * np.random.normal(0, 1, N)
    signal = s1 + s2 + 0.5 * s1 * s2 + noise
    _plot_signal(t, signal)

    kw = dict(nperseg=N // 10, noverlap=N // 20, nfft=next_fast_len(N // 2))
    freq1, freq2, bicoh = polycoherence(signal, fs, **kw)
    plot_polycoherence(freq1, freq2, bicoh)

    freq1, fre2, bispec = polycoherence(signal, fs, norm=None, **kw)
    plot_polycoherence(freq1, fre2, bispec)

    # freq1, freq2, bicoh = polycoherence(signal, fs, flim1=(0, 30), flim2=(0, 30), **kw)
    # plot_polycoherence(freq1, freq2, bicoh)







