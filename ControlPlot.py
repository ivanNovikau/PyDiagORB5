import Mix as mix
import matplotlib.pyplot as mpl
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)


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
