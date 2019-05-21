import Mix as mix
import ymath
import curve as crv
import matplotlib.pyplot as mpl
from matplotlib import animation
import matplotlib.ticker as ticker
import numpy as np
import plotly.offline as py
import plotly.graph_objs as go

FIG_SIZE_W = 14
FIG_SIZE_H = 9.5

# FIG_SIZE_W = 10
# FIG_SIZE_H = 6


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)


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


def plot_curves(curves, method=None):
    if method is 'mat':
        plot_curves_mat(curves)

    if method is 'plotly':
        plot_curves_plotly(curves)

    if method is None:
        # plot_curves_plotly(curves)
        plot_curves_mat(curves)


def plot_curves_mat(curves):
    # -> curves - class crv.Curves

    # number of curves
    ncurves = curves.n()

    # Build plots
    fig, ax = mpl.subplots(figsize=(FIG_SIZE_W, FIG_SIZE_H))

    for icrv in range(ncurves):
        curve = curves.list(icrv)

        y_res = curve.ys
        if y_res is None:
            continue

        if curves.flag_norm:
            y_res = ymath.find_norm(y_res)
        if curves.flag_semilogy:
            y_res = abs(y_res)

        if curves.flag_semilogy:
            ref_lines, = ax.semilogy(curve.xs, abs(y_res), curve.style)
        else:
            ref_lines, = ax.plot(curve.xs, y_res, curve.style)

        # set legend
        ref_lines.set_label(r'$' + curve.legend + '$')

        # set format for every line
        mpl.setp(ref_lines, linewidth=curve.width,
                 color=curve.color,
                 markersize=curve.markersize,
                 markerfacecolor=curve.markerfacecolor,
                 markeredgewidth=curve.width/2)

    # set labels:
    if len(curves.xlabel) is not 0:
        mpl.xlabel(r'$' + curves.xlabel + '$', fontsize=curves.fontS)
    if len(curves.ylabel) is not 0:
        mpl.ylabel(r'$' + curves.ylabel + '$', fontsize=curves.fontS)

    # axes ticks:
    if curves.xticks_labels is np.nan:
        mpl.xticks(curves.xticks) if curves.xticks is not np.nan else 0
    else:
        mpl.xticks(curves.xticks, curves.xticks_labels) if curves.xticks is not np.nan else 0

    if curves.yticks_labels is np.nan:
        mpl.yticks(curves.yticks) if curves.yticks is not np.nan else 0
    else:
        mpl.yticks(curves.yticks, curves.yticks_labels) if curves.yticks is not np.nan else 0

    # fontsize of axes ticks
    ax.xaxis.set_tick_params(labelsize=curves.fontS)
    ax.yaxis.set_tick_params(labelsize=curves.fontS)
    ax.xaxis.get_offset_text().set_fontsize(curves.fontS)
    ax.yaxis.get_offset_text().set_fontsize(curves.fontS)

    # format of axis labels
    mpl.ticklabel_format(axis='x', style=curves.x_style, scilimits=(-2, 2))

    if curves.flag_semilogy is False:
        mpl.ticklabel_format(axis='y', style=curves.y_style, scilimits=(-2, 2))

    # set limits:
    if curves.xlimits is not None:
        ax.set_xlim(curves.xlimits[0], curves.xlimits[-1])
    if curves.ylimits is not None:
        ax.set_ylim(curves.ylimits[0], curves.ylimits[-1])

    # set legend
    if curves.flag_legend:
        ax.legend(fontsize=curves.fontS, loc=curves.legend_position, facecolor=curves.legend_fcol)

    # set title
    mpl.title(r'$' + curves.title + '$', fontsize=curves.fontS)

    mpl.grid(True)
    mpl.show()


def plot_curves_plotly(curves):
    # -> curves - class crv.Curves

    # number of curves
    ncurves = curves.n()

    # Build plots
    data = []
    for icrv in range(ncurves):
        curve = curves.list(icrv)

        y_res = curve.ys
        if curves.flag_norm:
            y_res = y_res / np.max(np.abs(y_res))
        if curves.flag_semilogy:
            y_res = abs(y_res)

        # lines or marker
        modeL = {
            'o': 'markers',
            '-': 'lines',
            '-o': 'lines + markers',
            '--': 'lines',
            '--o': 'lines + markers',
            ':': 'lines',
            ':o': 'lines + markers',
            '-.': 'lines',
            '-.o': 'lines + markers',
        }.get(curve.style, 'lines')

        # style of the line
        dashL = {
            ':': 'dot',
            ':o': 'dot',
            '--': 'dash',
            '--o': 'dash',
            '-.': 'dashdot',
            '-.o': 'dashdot'
        }.get(curve.style, None)

        # data for the current curve
        trace = go.Scatter(
            x=curve.xs,
            y=y_res,
            name=r'$' + curve.legend + '$',  # legend
            mode=modeL,
            line=dict(
                width=curve.width,
                dash=dashL,
                color=curve.color
            ),
            marker=dict(
                color=curve.color,
                size=curve.markersize
            )
        )
        data.append(trace)

    if curves.flag_semilogy:
        type_y = 'log'
    else:
        type_y = 'linear'

    layout = go.Layout(
        title=r'$' + curves.title + '$',
        titlefont=dict(
            size=curves.fontS
        ),
        xaxis=dict(
            title=r'$' + curves.xlabel + '$',
            titlefont=dict(
                size=curves.fontS
            ),
            tickfont=dict(
                size=curves.axisFS,
            )
        ),
        yaxis=dict(
            type=type_y,
            title=r'$' + curves.ylabel + '$',
            titlefont=dict(
                size=curves.fontS
            ),
            tickfont=dict(
                size=curves.axisFS,
            )
        )
    )

    # build plot
    fig = go.Figure(data=data, layout=layout)
    py.iplot(fig)


def plot_curves_3d(curves, method=None):
    if method is 'mat':
        plot_curves_3d_mat(curves)

    if method is 'plotly':
        plot_curves_3d_plotly(curves)

    if method is None:
        plot_curves_3d_mat(curves)


def plot_curves_3d_plotly(curves):
    # -> curves - class crv.Curves

    # number of curves
    ncurves = curves.n()

    # data for every curve
    data = []
    for icrv in range(ncurves):
        curve = curves.list(icrv)
        trace = go.Contour(
            x=curve.xs,
            y=curve.ys,
            z=curve.zs,
            name=r'$' + curve.legend + '$',  # legend
        )
        data.append(trace)

    layout = go.Layout(
        title=curves.title,
        titlefont=dict(
            size=curves.fontS
        ),
        xaxis=dict(
            title=r'$' + curves.xlabel + '$',
            titlefont=dict(
                size=curves.fontS
            ),
            tickfont=dict(
                size=curves.axisFS,
            )
        ),
        yaxis=dict(
            title=r'$' + curves.ylabel + '$',
            titlefont=dict(
                size=curves.fontS
            ),
            tickfont=dict(
                size=curves.axisFS,
            )
        )
    )

    # build plot
    fig = go.Figure(data=data, layout=layout)
    py.iplot(fig)


def plot_curves_3d_mat(curves):
    # -> curves - class crv.Curves

    # number of curves
    ncurves = curves.n()

    # inititialization of the figure
    fig, ax = mpl.subplots(figsize=(FIG_SIZE_W, FIG_SIZE_H))

    # data from the first curve, that has to be 3d plot
    curve_one = curves.list(0)
    ZZ = curve_one.zs
    if curve_one.xs.ndim < 2:
        XX, YY = np.meshgrid(curve_one.xs, curve_one.ys)
    else:
        XX = curve_one.xs
        YY = curve_one.ys

    # build 3d plot
    # -> colormaps: hot, jet, pink, bone, bwr, RdYlBu, RdYlGn, RdBu, ocean, seismic, RdGy
    cs = ax.contourf(XX, YY, ZZ.T, levels=curve_one.levels, cmap=curve_one.colormap)

    def fmt(x, pos):
        if x == 0:
            return r'${:0.2f}$'.format(x)
        if 2 < np.log10(np.abs(x)) or np.log10(np.abs(x)) < -1:
            return r'${:.2e}$'.format(x)
        else:
            return r'${:0.2f}$'.format(x)
    cb = fig.colorbar(cs, shrink=0.8, extend='both', format=ticker.FuncFormatter(fmt))
    cb.ax.tick_params(labelsize=curves.fontS)

    # data from other curves, that should be lines
    for icrv in range(1, ncurves):
        curve = curves.list(icrv)

        # ref_lines, = ax.plot(curve.xs, curve.ys, curve.style)
        ref_lines = ax.errorbar(curve.xs, curve.ys,
                yerr=curve.ys_err, xerr=curve.xs_err, fmt=curve.style,
                elinewidth=curve.width/3, ecolor=curve.color)

        # set legend
        ref_lines.set_label(r'$' + curve.legend + '$')

        # set format for every line
        mpl.setp(ref_lines[0], linewidth=curve.width/2,
                 color=curve.color, markersize=curve.markersize/3,
                 markerfacecolor=curve.markerfacecolor,
                 markeredgewidth=curve.width/2)

    # set labels
    if len(curves.xlabel) is not 0:
        mpl.xlabel(r'$' + curves.xlabel + '$', fontsize=curves.fontS)
    if len(curves.ylabel) is not 0:
        mpl.ylabel(r'$' + curves.ylabel + '$', fontsize=curves.fontS)

    # font size of axes ticks
    ax.xaxis.set_tick_params(labelsize=curves.fontS)
    ax.yaxis.set_tick_params(labelsize=curves.fontS)
    ax.xaxis.get_offset_text().set_fontsize(curves.fontS)
    ax.yaxis.get_offset_text().set_fontsize(curves.fontS)

    # format of axis labels
    mpl.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))

    # set limits:
    if curves.xlimits is not None:
        ax.set_xlim(curves.xlimits[0], curves.xlimits[-1])
    if curves.ylimits is not None:
        ax.set_ylim(curves.ylimits[0], curves.ylimits[-1])

    # legend
    if ncurves > 1:
        # ax.legend(fontsize=curves.fontS, loc=curves.legend_position,
        #           framealpha=1, facecolor='grey')
        ax.legend(fontsize=curves.fontS, loc='upper left',
                  framealpha=1, facecolor='grey')

    # set title
    mpl.title(r'$' + curves.title + '$', fontsize=curves.fontS)


def animation_curves_2d(curves):
    # -> curves - class crv.Curves

    # number of curves
    ncurves = curves.n()

    fig, (ax) = mpl.subplots(1, 1, figsize=(FIG_SIZE_W, FIG_SIZE_H))
    ax.set_xlim(curves.xlimits[0], curves.xlimits[-1])
    if curves.flag_norm:
        ax.set_ylim(0, 1)
    line, = ax.plot([], [])
    line.set_data([], [])

    curve_one = curves.list(0)

    def animate(i, curves, line):
        if curves.flag_norm:
            ff = curve_one.zs[i, :]
            ff = ff / np.max(np.abs(ff))
        else:
            ff = curve_one.zs[i, :]
        line.set_data(curve_one.xs, ff)
        return line,

    # anim = animation.FuncAnimation(fig, animate, init_func=init,
    #                                frames=np.size(curve_one.ws), interval=20, blit=True)
    # anim = animation.FuncAnimation(fig, animate, init_func=init,
    #                                frames=10, interval=20, blit=True)
    anim = animation.FuncAnimation(fig, animate, 25, fargs=(curves, line),
                                   interval=100, blit=True)
    # anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    mpl.show()



