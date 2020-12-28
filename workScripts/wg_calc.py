import Mix as mix
import common as cm
import Global_variables as GLO
import numpy as np


def reload():
    mix.reload_module(mix)
    mix.reload_module(cm)
    mix.reload_module(GLO)


def wg_calc_n1_s1_chi1(dd, **oo):
    # *************************************
    # *** CALCULATE (w,g) of one n-mode ***
    # *************************************
    reload()

    tmin, tmax  = oo['tmin'], oo['tmax']
    n_mode      = oo['n_mode']
    chi_point   = oo['chi_point']
    s1          = oo['s1']
    filt_global = oo['filt_global']
    filt_freq   = oo['filt_freq']
    min_n_peaks = oo['min_n_peaks']
    threshold_w = oo['threshold_w']
    threshold_g = oo['threshold_g']
    n_samples   = oo['n_samples']

    variables = ['n1']
    if 'variables' in oo:
        variables = oo['variables']

    flag_stat = True
    if 'flag_stat' in oo:
        flag_stat = oo['flag_stat']

    sel_norm_wg = GLO.DEF_NORM_WG
    if 'sel_norm_wg' in oo:
        sel_norm_wg = oo['sel_norm_wg']

    # signal
    ch_signal = GLO.create_signals_dds(
        GLO.def_fields3d_n1,
        [dd],
        planes=['ts'],
        operations=['point-s'],
        domains=[s1],
        variables=variables
    )[0]
    ch_signal['n1'] = n_mode
    ch_signal['chi-point'] = chi_point

    # For FreqGrowth calculation
    oo_wg = {
        # BASIC
        't_work': [tmin, tmax],
        'flag_two_stages': True,
        'sel_norm_wg': sel_norm_wg,
        # FILTERING
        'filt_global': filt_global,
        'filt_freq': filt_freq,
        # STATISTICS
        'flag_stat': flag_stat,
        'n_samples': n_samples,
        'min_n_peaks': min_n_peaks,
        'threshold_w': threshold_w,
        'threshold_g': threshold_g,
    }

    # styling
    ff = dict(GLO.DEF_PLOT_FORMAT)
    ff.update({
        'flag_norm': True,
        'flag_semilogy': True
    })

    # plotting:
    oo = {
        'signal': ch_signal,
        'ff': ff,
        'flag_subplots': True,
    }
    cm.calc_wg(oo, oo_wg)


def fft_1d_n1_s1_chi1(oo):
    reload()

    # get input parameters
    dds = oo['dds']
    ss = oo['ss']
    n_modes = oo['n_modes']
    chi_points = oo.get('chi_points', None)

    n_signals = len(dds)

    variables = oo.get('variables', ['n1'] * n_signals)
    legends = oo.get('legends', [])
    w_start = oo.get('w_start', None)
    w_end = oo.get('w_end', None)
    styles = oo.get('styles', ['-', ':'])
    flag_norm = oo.get('flag_norm', False)
    x_fft_domain = oo.get('x_fft_domain', None)
    sel_norm_w = oo.get('sel_norm_w', 'frequency-wc')
    title = oo.get('title', None)

    # signals
    ch_signals = GLO.create_signals_dds(
        GLO.def_fields3d_n1,
        dds=dds,
        planes=['ts'] * n_signals,
        operations=['point-s'] * n_signals,
        domains=ss,
        variables=variables
    )
    for id1, one_signal in enumerate(ch_signals):
        one_signal['n1'] = n_modes[id1]
        one_signal['chi-point'] = chi_points[id1] if chi_points is not None else 0.0

    # ylabel
    ylabel = 'FFT\ [\Phi_{n}(s_1)]'
    if flag_norm:
        ylabel = 'norm.\ ' + ylabel

    # styling
    ff = dict(GLO.DEF_PLOT_FORMAT)
    ff.update({
        'xlabel': '\omega',
        'ylabel': ylabel,
        'title': title,
        'legends': legends,
        'styles': styles,
        'flag_norm': flag_norm,
    })

    # post-processing: 1d FFT
    oo_fft_one = dict(GLO.DEF_OPERATION_FFT_1D)

    if x_fft_domain is not None:
        oo_fft_one['domain'] = x_fft_domain

    oo_fft_one['oo_fft'].update({
        'flag_f2': False,  # plot two-sided spectrum or one sided
    })

    post_var_one = [oo_fft_one]
    oo_postprocessing = [post_var_one] * n_signals

    # plotting
    oo = {
        'signals': ch_signals,
        'ff': ff,
        'oo_postprocessing': oo_postprocessing,
        'sel_norm_x': sel_norm_w,
    }
    if w_start is not None:
        oo['x_start'] = w_start
    if w_end is not None:
        oo['x_end'] = w_end

    cm.plot_vars_1d(oo)


