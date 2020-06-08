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

    flag_stat = True
    if 'flag_stat' in oo:
        flag_stat = oo['flag_stat']

    # signal
    ch_signal = GLO.create_signals_dds(
        GLO.def_fields3d_n1,
        [dd],
        planes=['ts'],
        operations=['point-s'],
        domains=[s1],
    )[0]
    ch_signal['n1'] = n_mode
    ch_signal['chi-point'] = chi_point

    # For FreqGrowth calculation
    oo_wg = {
        # BASIC
        't_work': [tmin, tmax],
        'flag_two_stages': True,
        'sel_norm_wg': 'wc',
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

