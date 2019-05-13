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


def chi(species_name, dd, oo={}):
    # read heat flux
    rd.radial_heat_flux(dd)
    rd.nT_evol(dd, species_name)
    efluxw_rad = dd[species_name].efluxw_rad
    Lx = dd['Lx']
    rho_star_inv = Lx / 2.

    # signals:
    t = efluxw_rad['t'] # the same for all 1d signals
    s_flux  = efluxw_rad['s']
    s = dd[species_name].nT_evol['s']
    rad_flux = efluxw_rad['data']
    n = dd[species_name].nT_evol['n']
    gradT = dd[species_name].nT_evol['gradT']

    func_rad_flux_interp = \
        interpolate.interp1d(s_flux, rad_flux, axis=1)
    rad_flux = func_rad_flux_interp(s)

    chi = - rad_flux / (n * gradT)
    chi_norm = chi * rho_star_inv**2

    # intervals
    t_wci, ids_t = mix.get_array_oo(oo, t, 't')
    s, ids_s = mix.get_array_oo(oo, s, 's')
    rad_flux = mix.get_slice(rad_flux, ids_t, ids_s)
    chi_norm = mix.get_slice(chi_norm, ids_t, ids_s)

    line_t_wci = mix.test_array(t_wci, 't[wci^{-1}]', ':0.2e')

    # --- PLOTTING ---
    # efluxw_rad
    curves = crv.Curves().xlab('t[a/s]').ylab('s')\
        .tit(species_name + ':\ efluxw\_rad')
    curves.new('f').XS(t_wci).YS(s).ZS(rad_flux).lev(40)
    cpr.plot_curves_3d(curves)

    # chi/chi_B
    curves = crv.Curves().xlab('t[a/s]').ylab('s')\
        .tit(species_name + ':\ \chi/\chi_B')
    curves.new('f').XS(t_wci).YS(s).ZS(chi_norm).lev(40)
    cpr.plot_curves_3d(curves)


def chi_aver_st(species_name, dd, oo={}):
    project_name = oo.get('project_name', '')
    flag_save = oo.get('flag_save', False)
    save_name = oo.get('save_name', None)
    list_curves_load_avs = oo.get('list_curves_load_avs', {None})
    list_curves_load_avt = oo.get('list_curves_load_avt', {None})

    # read heat flux
    rd.radial_heat_flux(dd)
    rd.nT_evol(dd, species_name)
    efluxw_rad = dd[species_name].efluxw_rad
    Lx = dd['Lx']
    rho_star_inv = Lx / 2.

    # signals:
    t = efluxw_rad['t']  # the same for all 1d signals
    s_flux = efluxw_rad['s']
    rad_flux = efluxw_rad['data']

    s = dd[species_name].nT_evol['s']
    n = dd[species_name].nT_evol['n']
    gradT = dd[species_name].nT_evol['gradT']

    func_rad_flux_interp = \
        interpolate.interp1d(s_flux, rad_flux, axis=1)
    rad_flux = func_rad_flux_interp(s)

    chi = - rad_flux / (n * gradT)
    chi_norm = chi * rho_star_inv ** 2

    # intervals
    t_int, ids_t_int = mix.get_array_oo(oo, t, 't')
    s_int, ids_s_int = mix.get_array_oo(oo, s, 's')

    # time evolution in different radial domains:
    ns_av = oo.get('ns_av', 1)

    data_s_av, lines_s_av = {}, {}
    for is_av in range(ns_av):
        s_av, ids_s_av = mix.get_array_oo(oo, s, 's_av{:d}'.format(is_av + 1))
        temp = mix.get_slice(chi_norm, ids_t_int, ids_s_av)
        data_s_av[is_av]  = np.mean(temp, axis=1)
        lines_s_av[is_av] = mix.test_array(s_av, 's', ':0.3f')

    curves_avs = crv.Curves().xlab('t[wci^{-1}]').tit(species_name + '\ <chi/chi_B>_s')
    for is_av in range(ns_av):
        curves_avs.new() \
            .XS(t_int)\
            .YS(data_s_av[is_av])\
            .leg(project_name + ':\ ' + lines_s_av[is_av])
    for curves_load_avs in list_curves_load_avs:
        curves_avs.load(curves_load_avs)
    curves_avs.set_colors_styles()
    cpr.plot_curves(curves_avs)

    # average in time
    t_av_wci, ids_t_av = mix.get_array_oo(oo, t, 't_av')
    data_av = mix.get_slice(chi_norm, ids_t_av, ids_s_int)
    data_av = np.mean(data_av, axis=0)

    ## FOR NORMALIZATION
    # data_av[np.isnan(data_av)] = -np.inf
    # data_av[np.isinf(data_av)] = -np.inf
    # data_av = data_av / np.max(data_av)
    # data_av[np.isinf(data_av)] = np.nan

    line_t_av_wci = mix.test_array(t_av_wci, 't[wci^{-1}]', ':0.2e')

    curves_avt = crv.Curves().xlab('s') \
        .tit(species_name + ':\ <chi/chi_B>_t:\ ')
    curves_avt.new()\
        .XS(s_int)\
        .YS(data_av)\
        .leg(project_name + ':\ ' + line_t_av_wci)
    for curves_load_avt in list_curves_load_avt:
        curves_avt.load(curves_load_avt)
    curves_avt.set_colors_styles()
    cpr.plot_curves(curves_avt)

    # save data
    if flag_save:
        if 'saved_data' not in dd:
            dd['saved_data'] = {}
        dd['saved_data'][save_name + '-avs'] = curves_avs
        dd['saved_data'][save_name + '-avt'] = curves_avt


def fft_chi(species_name, dd, oo={}):
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

    chi = - rad_flux / (n * gradT)
    chi_norm = chi * rho_star_inv**2

    # change normalization of the time grid
    cs = dd['cs']
    wc = dd['wc']
    a0 = dd['a0']
    norm_cs = cs / a0
    t = t * norm_cs / wc

    # intervals
    t_wci, ids_t = mix.get_array_oo(oo, t, 't')
    s, ids_s = mix.get_array_oo(oo, s, 's')
    rad_flux = mix.get_slice(rad_flux, ids_t, ids_s)
    chi_norm = mix.get_slice(chi_norm, ids_t, ids_s)

    line_t_wci = mix.test_array(t_wci, 't[wci^{-1}]', ':0.2e')

    # FFT
    w, _, nw = ymath.w_for_fft(t)
    w = w * 2 * np.pi
    fft_rad_flux, _ = ymath.fft_y(rad_flux, nw, 0)
    fft_chi_norm, _ = ymath.fft_y(chi_norm, nw, 0)

    w, ids_w = mix.get_array_oo(oo, w, 'w')
    fft_rad_flux = mix.get_slice(fft_rad_flux, ids_w)
    fft_chi_norm = mix.get_slice(fft_chi_norm, ids_w)

    # --- PLOTTING ---
    # efluxw_rad
    curves = crv.Curves().xlab('s').ylab('\omega[cs/a]')\
        .tit(species_name + ':\ FFT:\ efluxw\_rad:')\
        .titn(line_t_wci)
    curves.new('f').XS(s).YS(w).ZS(fft_rad_flux.T).lev(40)
    cpr.plot_curves_3d(curves)

    # chi/chi_B
    curves = crv.Curves().xlab('s').ylab('\omega[cs/a]')\
        .tit(species_name + ':\ FFT:\ \chi/\chi_B:')\
        .titn(line_t_wci)
    curves.new('f').XS(s).YS(w).ZS(fft_chi_norm.T).lev(40)
    cpr.plot_curves_3d(curves)


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