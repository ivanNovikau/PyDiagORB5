import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import general as gn
import numpy as np
from scipy import interpolate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(gn)


def chi(species_name, dd):
    if 'chi' in dd:
        return

    # read heat flux
    rd.radial_heat_flux(dd)
    rd.nT_evol(dd, species_name)
    efluxw_rad = dd[species_name].efluxw_rad
    Lx = dd['Lx']
    rho_star_inv = Lx / 2.

    # signals:
    t = efluxw_rad['t']  # the same for all 1d signals
    s_flux  = efluxw_rad['s']
    s = dd[species_name].nT_evol['s']
    rad_flux = efluxw_rad['data']
    n = dd[species_name].nT_evol['n']
    gradT = dd[species_name].nT_evol['gradT']

    func_rad_flux_interp = \
        interpolate.interp1d(s_flux, rad_flux, axis=1)
    rad_flux = func_rad_flux_interp(s)

    chi_var = - rad_flux / (n * gradT)
    chi_norm = chi_var * rho_star_inv**2

    dd['chi'] = {
        'data': chi_var,
        'data_norm': chi_norm,
        's': s,
        't': t
    }


def choose_var(dd, oo):
    species_name = oo.get('species_name', 'deuterium')

    opt_var = oo.get('opt_var', 'chi_norm')
    vvar, t, s, tit_var = [], [], [], []
    if opt_var == 'chi_norm':
        chi(species_name, dd)
        vvar = dd['chi']['data_norm']
        tit_var = species_name + ':\ \chi/\chi_B'
        t = dd['chi']['t']
        s = dd['chi']['s']
    if opt_var == 'efluxw_rad':
        efluxw_rad = dd[species_name].efluxw_rad
        vvar = efluxw_rad['data']
        tit_var = species_name + ':\ efluxw\_rad'
        t = efluxw_rad['t']
        s = efluxw_rad['s']

    res = {
        'var': vvar,
        's': s,
        't': t,
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
    vvar, s, t, tit_var = out['var'], out['s'], out['t'], out['tit']

    # --- averaging in time ---
    oo_avt = dict(oo)
    # oo_avt.update({
    #     'vars': [vvar, vvar],
    #     'ts': [t, t], 'ss': [s, s],
    #     'opts_av': ['rms', 'mean'],
    #     'tit': tit_var, 'vars_names': ['', '']
    # })
    oo_avt.update({
        'vars': [vvar],
        'ts': [t], 'ss': [s],
        'opts_av': ['rms'],
        'tit': tit_var, 'vars_names': ['']
    })
    gn.plot_avt(dd, oo_avt)

    # --- averaging in space ---
    oo_avs = dict(oo)
    oo_avs.update({
        'vars': [vvar], 'ts': [t], 'ss': [s],
        'opts_av': ['mean'],
        'tit': tit_var, 'vars_names': ['']
    })
    gn.plot_avs(dd, oo_avs)


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


def T(species_name, dd, oo={}):
    rd.nT_evol(dd, species_name)

    # signals:
    t = dd[species_name].nT_evol['t']
    s = dd[species_name].nT_evol['s']
    T = dd[species_name].nT_evol['T']

    # intervals
    t_wci, ids_t = mix.get_array_oo(oo, t, 't')
    s, ids_s = mix.get_array_oo(oo, s, 's')
    T     = mix.get_slice(T, ids_t, ids_s)

    # *** plotting ***
    curves = crv.Curves().xlab('t[a/s]').ylab('s')\
        .tit(species_name + ':\ temperature')
    curves.new('f').XS(t_wci).YS(s).ZS(T).lev(40)
    cpr.plot_curves_3d(curves)


def grad_T(species_name, dd, oo={}):
    rd.nT_evol(dd, species_name)

    # signals:
    t = dd[species_name].nT_evol['t']
    s = dd[species_name].nT_evol['s']
    gradT = dd[species_name].nT_evol['gradT']

    # intervals
    t_wci, ids_t = mix.get_array_oo(oo, t, 't')
    s, ids_s = mix.get_array_oo(oo, s, 's')
    gradT = mix.get_slice(gradT, ids_t, ids_s)

    # *** plotting ***
    # grad T(t)
    curves = crv.Curves().xlab('t[a/s]').ylab('s')\
        .tit(species_name + ':\ \\nabla T')
    curves.new('f').XS(t_wci).YS(s).ZS(gradT).lev(40)
    cpr.plot_curves_3d(curves)


def n(species_name, dd, oo={}):
    rd.nT_evol(dd, species_name)

    # signals:
    t = dd[species_name].nT_evol['t']
    s = dd[species_name].nT_evol['s']
    n = dd[species_name].nT_evol['n']

    # intervals
    t_wci, ids_t = mix.get_array_oo(oo, t, 't')
    s, ids_s = mix.get_array_oo(oo, s, 's')
    n = mix.get_slice(n, ids_t, ids_s)

    # *** plotting ***
    # n(t)
    curves = crv.Curves().xlab('t[a/s]').ylab('s')\
        .tit(species_name + ':\ density')
    curves.new('f').XS(t_wci).YS(s).ZS(n).lev(40)
    cpr.plot_curves_3d(curves)


def grad_n(species_name, dd, oo={}):
    rd.nT_evol(dd, species_name)

    # signals:
    t = dd[species_name].nT_evol['t']
    s = dd[species_name].nT_evol['s']
    gradn = dd[species_name].nT_evol['gradn']

    # intervals
    t_wci, ids_t = mix.get_array_oo(oo, t, 't')
    s, ids_s = mix.get_array_oo(oo, s, 's')
    gradn = mix.get_slice(gradn, ids_t, ids_s)

    # *** plotting ***
    # grad n(t)
    curves = crv.Curves().xlab('t[a/s]').ylab('s')\
        .tit(species_name + ':\ \\nabla n')
    curves.new('f').XS(t_wci).YS(s).ZS(gradn).lev(40)
    cpr.plot_curves_3d(curves)


def efluxw_rad_2d(species_name, dd, oo={}):
    rd.radial_heat_flux(dd)
    efluxw_rad = dd[species_name].efluxw_rad

    # intervals
    t_wci, ids_t = mix.get_array_oo(oo, efluxw_rad['t'], 't')
    s, ids_s = mix.get_array_oo(oo, efluxw_rad['s'], 's')
    data = mix.get_slice(efluxw_rad['data'], ids_t, ids_s)

    line_t_wci = mix.test_array(t_wci, 't[wci^{-1}]', ':0.2e')

    # change normalization of the time grid
    cs = dd['cs']
    wc = dd['wc']
    a0 = dd['a0']
    norm_cs = cs / a0
    t_cs = t_wci * norm_cs / wc

    # find FFT of the signal:
    w, _, nw = ymath.w_for_fft(t_cs)
    w = w * 2 * np.pi
    ff, _ = ymath.fft_y(data, nw, 0)

    # intervals for the frequency:
    w, ids_w = mix.get_array_oo(oo, w, 'w')
    ff = mix.get_slice(ff, ids_w)

    # time evolution in different domains:
    ns_av = oo.get('ns_av', 1)

    data_s_av = {}
    lines_s_av = {}
    for is_av in range(ns_av):
        s_av, ids_s_av = mix.get_array_oo(oo, efluxw_rad['s'],
                                          's_av{:d}'.format(is_av + 1))
        temp = mix.get_slice(efluxw_rad['data'], ids_t, ids_s_av)
        data_s_av[is_av] = np.mean(temp, axis=1)
        lines_s_av[is_av] = mix.test_array(s_av, 's', ':0.3f')

    curves = crv.Curves().xlab('t[wci^{-1}]')\
        .tit(species_name + '\ <heat\ flux>_s')
    for is_av in range(ns_av):
        curves.new('y{:d}'.format(is_av)) \
            .XS(t_wci).YS(data_s_av[is_av]).leg(lines_s_av[is_av])
    cpr.plot_curves(curves)

    # average in time
    t_av_wci, ids_t_av = mix.get_array_oo(oo, efluxw_rad['t'], 't_av')
    data_av = mix.get_slice(efluxw_rad['data'], ids_t_av, ids_s)
    data_av = np.mean(data_av, axis=0)

    line_t_av_wci = mix.test_array(t_av_wci, 't[wci^{-1}]', ':0.2e')

    curves = crv.Curves().xlab('s') \
        .tit('norm.\ <heat\ flux>_t:\ ') \
        .tit(species_name + ':\ ' + line_t_av_wci)
    curves.flag_norm = True
    curves.new('y').XS(s).YS(data_av)
    cpr.plot_curves(curves)

    # fourier spectrum
    curves = crv.Curves().xlab('s').ylab('\omega[cs/a]') \
        .tit(species_name + ':\ FFT:\ heat\ flux:\ ') \
        .tit(line_t_wci)
    curves.new('f').XS(s).YS(w).ZS(ff.T).lev(40)
    cpr.plot_curves_3d(curves)