import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import matplotlib.pyplot as mpl
import numpy as np
import types
from matplotlib import animation, ticker
from IPython.display import HTML
import matplotlib.colors


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)
    mix.reload_module(GLO)


def plot_curves(curves):
    if curves.is_empty():
        return

    # number of curves
    ncurves = curves.n()

    # Build plots
    fig, ax = mpl.subplots(figsize=(GLO.FIG_SIZE_W, GLO.FIG_SIZE_H))
    axes = mpl.gca()

    # set curves
    for icrv in range(ncurves):
        curve = curves.list(icrv)

        y_res = curve.ys
        if y_res is None:
            continue

        if curve.flag_hist:
            ax.hist(y_res, curve.xs, alpha=curve.pr_alpha,
                    label=r'\boldmath $' + curve.legend  + '$',
                    density=True)
        else:
            if curves.flag_norm:
                y_res = ymath.find_norm(y_res, curve.data_norm_to)
            if curves.flag_semilogy:
                y_res = abs(y_res)

            if not curve.flag_errorbar:
                if curves.flag_semilogy:
                    ref_lines, = ax.semilogy(curve.xs, abs(y_res), curve.style)
                else:
                    ref_lines, = ax.plot(curve.xs, y_res, curve.style)

                if curve.style == ':':
                    ref_lines.set_dashes([0.5, 0.4])
                mpl.setp(ref_lines, linewidth=curve.width,
                         color=curve.color,
                         markersize=curve.markersize,
                         markerfacecolor=curve.markerfacecolor,
                         markeredgewidth=curve.width / 2)
            else:
                ref_lines = ax.errorbar(curve.xs, curve.ys,
                                        yerr=curve.ys_err, xerr=curve.xs_err, fmt=curve.style,
                                        elinewidth=curve.width / 3, ecolor=curve.color)
                mpl.setp(ref_lines[0], linewidth=curve.width,
                         color=curve.color,
                         markersize=curve.markersize,
                         markerfacecolor=curve.markerfacecolor,
                         markeredgewidth=curve.width / 2)

            # set legend
            if curve.legend == "_":
                ref_lines.set_label("_")
            else:
                ref_lines.set_label(r'\boldmath $' + curve.legend + '$')

    # format the plot
    format_plot(fig, ax, axes, curves)


def plot_curves_3d(curves):
    if curves.is_empty():
        return

    # number of curves
    ncurves = curves.n()

    # initialization of the figure
    fig, ax = mpl.subplots(figsize=(GLO.FIG_SIZE_W, GLO.FIG_SIZE_H))
    axes = mpl.gca()

    # data from the first curve, that has to be 3d plot
    curve_one = curves.list(0)
    ZZ = curve_one.zs
    if curve_one.xs.ndim < 2:
        XX, YY = np.meshgrid(curve_one.xs, curve_one.ys)
    else:
        XX = curve_one.xs
        YY = curve_one.ys

    # --- contour plot ---
    divnorm = None
    if curve_one.colormap_center is not None:
        divnorm = matplotlib.colors.DivergingNorm(vcenter=curve_one.colormap_center)
    cs = ax.contourf(XX, YY, ZZ.T,
                     levels=curve_one.levels, cmap=curve_one.colormap,
                     norm=divnorm)

    # color bar
    if curves.flag_colorbar:
        cb = fig.colorbar(cs, shrink=0.8, extend='both')
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((0, 0))
        cb.ax.tick_params(labelsize=GLO.FONT_SIZE_TICKS)
        cb.ax.yaxis.get_offset_text().set_fontsize(GLO.FONT_SIZE_ORDER)

        register_offset(cb.ax.yaxis, bottom_offset)
        cb.update_ticks()

    # --- other curves ---
    for icrv in range(1, ncurves):
        curve = curves.list(icrv)

        if not curve.flag_errorbar:
            ref_lines, = ax.plot(curve.xs, curve.ys, curve.style)
            if curve.style == ':':
                ref_lines.set_dashes([0.6, 0.6])
            mpl.setp(ref_lines, linewidth=curve.width/2,
                     color=curve.color,
                     markersize=curve.markersize / 2,
                     markerfacecolor=curve.markerfacecolor,
                     markeredgewidth=curve.width / 2)
        else:
            ref_lines = ax.errorbar(curve.xs, curve.ys,
                    yerr=curve.ys_err, xerr=curve.xs_err, fmt=curve.style,
                    elinewidth=curve.width/3, ecolor=curve.color)
            mpl.setp(ref_lines[0], linewidth=curve.width / 2,
                     color=curve.color, markersize=curve.markersize / 3,
                     markerfacecolor=curve.markerfacecolor,
                     markeredgewidth=curve.width / 2)

        # set legend
        ref_lines.set_label(r'\boldmath $' + curve.legend + '$')

    # format the plot
    format_plot(fig, ax, axes, curves, flag_2d=True)


def animation_1d(curves):
    # WORKS, but in development

    # PLOTTING along Y, ANIMATION along X, DATA is Z

    # number of curves
    ncurves = curves.n()
    ngeoms = curves.n_geoms

    # Build plots
    fig, ax = mpl.subplots(figsize=(GLO.FIG_SIZE_W, GLO.FIG_SIZE_H))
    axes = mpl.gca()

    # set limits:
    if curves.xlimits is not None:
        ax.set_xlim(curves.xlimits[0], curves.xlimits[-1])
    if curves.ylimits is not None:
        ax.set_ylim(curves.ylimits[0], curves.ylimits[-1])

    # set labels:
    if curves.xlabel is not None:
        mpl.xlabel(r'\boldmath $' + curves.xlabel + '$', fontsize=GLO.FONT_SIZE * 1.7)
    if curves.ylabel is not None:
        mpl.ylabel(r'\boldmath $' + curves.ylabel + '$', fontsize=GLO.FONT_SIZE * 1.7)

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
    ax.xaxis.set_tick_params(labelsize=GLO.FONT_SIZE)
    ax.yaxis.set_tick_params(labelsize=GLO.FONT_SIZE)
    ax.xaxis.get_offset_text().set_fontsize(GLO.FONT_SIZE)
    ax.yaxis.get_offset_text().set_fontsize(GLO.FONT_SIZE)

    # format of axis labels
    mpl.ticklabel_format(axis='x', style=curves.x_style, scilimits=(-2, 2))
    if curves.flag_maxlocator:
        ax.xaxis.set_major_locator(mpl.MaxNLocator(curves.maxlocator))

    if curves.flag_semilogy is False:
        mpl.ticklabel_format(axis='y', style=curves.y_style, scilimits=(-2, 2))

    # set legend
    if curves.flag_legend:
        ax.legend(fontsize=GLO.FONT_SIZE * GLO.LEG_SCALE, loc=curves.legend_position,
                  facecolor=curves.legend_fcol)

    # set title
    if curves.title is not None:
        mpl.title(r'\boldmath $' + curves.title + '$', fontsize=GLO.FONT_SIZE * 1.5, pad='18', usetex=True)
    if GLO.FLAG_LATEX:
        mpl.rc('text', usetex=True)
        mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
        mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

    # set grid
    mpl.grid(True)

    # set empty curves
    ref_lines = [None] * ncurves
    for icrv in range(ncurves):
        curve = curves.list(icrv)
        if not curve.flag_errorbar:
            if curves.flag_semilogy:
                ref_lines[icrv], = ax.semilogy([], [], curve.style)
            else:
                ref_lines[icrv], = ax.plot([], [], curve.style)

            if curve.style == ':':
                ref_lines[icrv].set_dashes([0.5, 0.4])
            mpl.setp(ref_lines[icrv], linewidth=curve.width,
                     color=curve.color,
                     markersize=curve.markersize,
                     markerfacecolor=curve.markerfacecolor,
                     markeredgewidth=curve.width / 2)
        else:
            ref_lines[icrv] = ax.errorbar([], [],
                                    yerr=curve.ys_err, xerr=curve.xs_err, fmt=curve.style,
                                    elinewidth=curve.width / 3, ecolor=curve.color)
            mpl.setp(ref_lines[icrv][0], linewidth=curve.width,
                     color=curve.color,
                     markersize=curve.markersize,
                     markerfacecolor=curve.markerfacecolor,
                     markeredgewidth=curve.width / 2)

        # set legend
        if curve.legend == "_":
            ref_lines[icrv].set_label("_")
        else:
            ref_lines[icrv].set_label(r'\boldmath $' + curve.legend + '$')

    # # draw geometrical figures:
    # for igeom in range(ngeoms):
    #     one_geom = curves.list_geoms[igeom]
    #     one_geom.draw(mpl, ax, axes, {})

    # if FLAG_LATEX:
    #     fig.tight_layout()

    # initialization function: plot the background of each frame
    def init():
        for icrvL in range(ncurves):
            ref_lines[icrvL].set_data([], [])
        return ref_lines,  # !!!

    # animation function. This is called sequentially
    def animate(i, Y_res, Z_res):
        for icrvL in range(ncurves):
            if i < np.shape(z_res)[0]:
                ref_lines[icrvL].set_data(Y_res[icrvL][:], Z_res[icrvL][i, :])
        return ref_lines,  # !!!

    nx_max = np.zeros(ncurves)
    for icrv in range(ncurves):
        nx_max[icrv] = np.size(curves.list(icrv).xs)
    nx_max = int(np.max(nx_max))

    Y_res, Z_res = [], []
    for icrv in range(ncurves):
        curve = curves.list(icrv)

        z_res = curve.zs
        if z_res is None:
            continue

        if curves.flag_norm:
            z_res = ymath.find_norm(z_res, curve.data_norm_to)
        if curves.flag_semilogy:
            z_res = abs(z_res)
        Z_res.append(z_res)
        Y_res.append(curve.ys)

    anim = animation.FuncAnimation(
        fig,
        animate, fargs=(Y_res, Z_res),
        init_func=init,
        frames=nx_max, interval=20, blit=True
    )

    HTML(anim.to_html5_video())


def format_plot(fig, ax, axes, curves, flag_2d=False):
    ncurves = curves.n()
    ngeoms = curves.n_geoms
    ntexts = len(curves.list_text)

    # set labels:
    if curves.xlabel is not None:
        mpl.xlabel(r'\boldmath $' + curves.xlabel + '$', fontsize=GLO.FONT_SIZE_LABELS)
    if curves.ylabel is not None:
        mpl.ylabel(r'\boldmath $' + curves.ylabel + '$', fontsize=GLO.FONT_SIZE_LABELS)

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
    ax.xaxis.set_tick_params(labelsize=GLO.FONT_SIZE_TICKS)
    ax.yaxis.set_tick_params(labelsize=GLO.FONT_SIZE_TICKS)
    ax.xaxis.get_offset_text().set_fontsize(GLO.FONT_SIZE_ORDER)
    ax.yaxis.get_offset_text().set_fontsize(GLO.FONT_SIZE_ORDER)
    register_offset(ax.yaxis, top_offset)

    # format of axis labels
    mpl.ticklabel_format(axis='x', style=curves.x_style, scilimits=(-2, 2))
    if curves.flag_maxlocator:
        ax.xaxis.set_major_locator(mpl.MaxNLocator(curves.maxlocator))

    if curves.flag_semilogy is False:
        mpl.ticklabel_format(axis='y', style=curves.y_style, scilimits=(-2, 2))

    # set limits:
    if curves.xlimits is not None:
        ax.set_xlim(curves.xlimits[0], curves.xlimits[-1])
    if curves.ylimits is not None:
        ax.set_ylim(curves.ylimits[0], curves.ylimits[-1])

    # set legend
    if ncurves > 1 and curves.flag_legend:
        ax.legend(fontsize=GLO.FONT_SIZE * GLO.LEG_SCALE, loc=curves.legend_position,
                  facecolor=curves.legend_fcol,
                  labelspacing=0.1, handlelength=1, handletextpad=0.4)

    # set title
    if curves.title is not None:
        mpl.title(r'\boldmath $' + curves.title + '$', fontsize=GLO.FONT_SIZE_TITLE, pad='18', usetex=True)
    if GLO.FLAG_LATEX:
        mpl.rc('text', usetex=True)
        mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
        mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

    # draw geometrical figures:
    for igeom in range(ngeoms):
        one_geom = curves.list_geoms[igeom]
        if one_geom is not None:
            one_geom.draw(mpl, ax, axes, {})

    # add text:
    for itext in range(ntexts):
        loc_text = curves.list_text[itext]
        mpl.text(
            loc_text.x,
            loc_text.y,
            r'\boldmath $' + loc_text.line + '$',
            fontsize=GLO.FONT_SIZE * loc_text.coef_width,
            color=loc_text.color
        )

    # set grid
    if not flag_2d:
        mpl.grid(True)

    if GLO.FLAG_LATEX:
        fig.tight_layout()


def bottom_offset(self, bboxes, bboxes2):
    bottom = self.axes.bbox.ymin
    self.offsetText.set(va="top", ha="left")
    self.offsetText.set_position(
        (0, bottom - self.OFFSETTEXTPAD * 8 * self.figure.dpi / 72.0))


def top_offset(self, bboxes, bboxes2):
    top = self.axes.bbox.ymax
    self.offsetText.set(va="bottom", ha="left")
    self.offsetText.set_position(
        (0, top + self.OFFSETTEXTPAD * 6 * self.figure.dpi / 72.0))


def register_offset(axis, func):
    axis._update_offset_text_position = types.MethodType(func, axis)


