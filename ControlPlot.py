import Mix as mix
import curve as crv
import matplotlib.pyplot as mpl
from matplotlib import rc as mrc
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(crv)


def plot_x1(x, y, oo={}):
    fig, ax = mpl.subplots(figsize=(5, 3))
#    ax.plot(x, y)
    ax.semilogy(x, abs(y))


def plot_x1x2(x, y, z, oo={}):
    fig, ax = mpl.subplots(figsize=(5, 3))

    if x.ndim < 2:
        X, Y = np.meshgrid(x, y)
    else:
        X = x
        Y = y
    cs = ax.contour(X, Y, z)
    fig.colorbar(cs, shrink=0.8, extend='both')
    mpl.xlabel(oo.get('xlabel', ''))
    mpl.ylabel(oo.get('ylabel', ''))


def plot_curves(curves):
    # -> curves - class crv.Curves

    # # change font size
    # map_font = {'size': curves.fontS}
    # mrc('font', **map_font)

    # number of curves
    ncurves = curves.n()

    # Build plots
    fig, ax = mpl.subplots(figsize=(10, 6))
    for icrv in range(ncurves):
        curve = curves.list(icrv)
        if curves.flag_semilogy:
            ref_lines, = ax.semilogy(curve.xs, abs(curve.ys), curve.style)
        else:
            ref_lines, = ax.plot(curve.xs, curve.ys, curve.style)

        # set legend
        ref_lines.set_label(r'$' + curve.legend + '$')

        # set format for every line
        mpl.setp(ref_lines, linewidth=curve.width,
                 color=curve.color, markersize=curve.markersize)

    # set labels:
    mpl.xlabel(r'$' + curves.xlabel + '$', fontsize=curves.fontS)
    mpl.ylabel(r'$' + curves.ylabel + '$', fontsize=curves.fontS)

    # fontsize of axes ticks
    ax.xaxis.set_tick_params(labelsize=curves.fontS)
    ax.yaxis.set_tick_params(labelsize=curves.fontS)

    ax.legend(fontsize = curves.fontS)
    mpl.grid(True)




