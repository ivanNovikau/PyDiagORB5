import Mix as mix
import common as cm
import Global_variables as GLO
import numpy as np
import matplotlib.pyplot as mpl
from matplotlib import animation


def reload():
    mix.reload_module(mix)
    mix.reload_module(cm)
    mix.reload_module(GLO)


def get_pot_rz(dd, t, variable, n_mode, title):
    import fields3d
    mix.reload_module(fields3d)

    get_f = lambda t1: fields3d.choose_one_var_rz({
        'dd': dd,
        'variable': variable,
        'n1': n_mode,
        't-point': t1
    })

    init = get_f(t[0])
    res = {
        't': t,
        'title': title,
        'r': init['r'],
        'z': init['z'],
        's': init['s'],
        'chi': init['chi'],
    }

    data = init['data']
    res['data'] = np.zeros(
        [len(res['t']), np.size(data, 0), np.size(data, 1)]
    )
    res['data'][0] = data

    for it_d, t_point in enumerate(t[1:]):
        res['data'][it_d] = get_f(t_point)['data']
    return res


def anim_2d(two_signals, file_to_save):
    ncols = 1 if two_signals[1] is None else 2
    fig, axes = mpl.subplots(ncols=ncols, nrows=1)
    if ncols == 1:
        axes = [axes]

    mpl.xlabel('R')
    mpl.ylabel('Z')

    oplot = [None] * ncols

    # for id_col in range(ncols):
    #     ss = two_signals[id_col]
    #     ax = axes[id_col]
    #
    #     ax.set_title(
    #         ss['title'] + ':\ t = {:0.3e}'.format(ss['t'][0])
    #     )
    #     oplot[id_col] = ax.contourf(
    #         ss['r'], ss['z'], ss['data'][0]
    #     )

    def set_plot(id_t):
        for id_col in range(ncols):
            ss = two_signals[id_col]
            ax = axes[id_col]

            ax.set_title(
                ss['title'] + ': t = {:0.3e}'.format(ss['t'][id_t])
            )
            oplot[id_col] = axes[id_col].contourf(
                ss['r'],
                ss['z'],
                ss['data'][id_t].T
            )
        return oplot,

    anim = animation.FuncAnimation(
        fig,
        set_plot,
        frames=100, interval=100, blit=True,
    )
    anim.save(file_to_save)


