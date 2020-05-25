import Mix as mix
import common as cm
import Global_variables as GLO
import numpy as np


def reload():
    mix.reload_module(mix)
    mix.reload_module(cm)
    mix.reload_module(GLO)


def plot_safety_factor(dds, legends, flag_ivis=False):
    s_domain = [0.0, 1.0]
    styles = ['-', ':']

    q_signal_res = GLO.create_signals_dds(
        GLO.def_safety_factor,
        dds=dds,
    )
    ff = dict(GLO.DEF_PLOT_FORMAT)
    ff.update({
        'xlabel': 's',
        'ylabel': 'safety factor',
        'styles': styles,
        'legends': legends,
        'flag_ivis': flag_ivis,
    })
    oo_plot = {
        'signals': q_signal_res,
        'ff': ff,
        'x_start': s_domain[0],
        'x_end': s_domain[1],
    }
    cm.plot_vars_1d(oo_plot)


def plot_magnetic_field_configuration(dd, flag_ivis=False):
    oo_format = {
        's_domain': [0.0, 1.0],
        'flag_subplots': False,
        'q_fluxes': [1.0, 3.0],
        'flag_graphic': True,
        'flag_plot_q': False,
        'flag_ivis': flag_ivis,
    }
    cm.Bq_equil(dd, oo_format=oo_format)



