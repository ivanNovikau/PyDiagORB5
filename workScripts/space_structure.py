import Mix as mix
import common as cm
import Global_variables as GLO
import numpy as np


def reload():
    mix.reload_module(mix)
    mix.reload_module(cm)
    mix.reload_module(GLO)


def poloidal_n1_t1(dd, **oo):
    n_mode = oo['n_mode']
    t_point = oo['t_point']
    plane = oo['plane']  # 'rz' or 'schi'
    slims = oo['slims']  # [s_left, s_right], always along s-grid, [0,1]
    chilims = oo['chilims']  # in radians, counterclockwise, always along chi-grid
        # [np.pi/2, None] - from pi/2 to 2pi
        # [-np.pi/2, np.pi/2] - low-field-side
        # etc.

    # signal
    signal = GLO.create_signals_dds(
        GLO.def_fields3d_n1_schi,
        [dd],
        planes=[plane],
    )[0]
    signal['t-point'] = t_point
    signal['n1'] = n_mode

    # styling
    ff = dict(GLO.DEF_PLOT_FORMAT)
    ff.update({
        'xlabel': 's' if plane == 'schi' else 'R',
        'ylabel': '\chi' if plane == 'schi' else 'Z'
    })

    # plotting
    oo = {
        'signal': signal,
        'ff': ff,
    }
    if slims[0] is not None:
        oo['s_start'] = slims[0]
    if slims[-1] is not None:
        oo['s_end'] = slims[-1]
    if chilims[0] is not None:
        oo['chi_start'] = chilims[0]
    if chilims[-1] is not None:
        oo['chi_end'] = chilims[-1]

    cm.plot_vars_2d(oo)