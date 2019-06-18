import Mix as mix
import Constants as cst
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cst)


class Geom:
    name = ''
    geom_type = 'NONE'

    def draw(self, mpl, ax, axes, oo):
        return


class Fill(Geom):
    geom_type = 'FILL'

    xs = None
    ys = None
    color = 'b'
    alpha = 0.2

    def draw(self, mpl, ax, axes, oo):
        xlims = axes.get_xlim()
        ylims = axes.get_ylim()

        xs_plot = np.zeros(np.size(self.xs))
        ys_plot = np.zeros(np.size(self.xs))
        for ix in range(np.size(self.xs)):
            if self.xs[ix] == 'liml':
                xs_plot[ix] = xlims[0]
            elif self.xs[ix] == 'limr':
                xs_plot[ix] = xlims[-1]
            else:
                xs_plot[ix] = self.xs[ix]

            if self.ys[ix] == 'limb':
                ys_plot[ix] = ylims[0]
            elif self.ys[ix] == 'limu':
                ys_plot[ix] = ylims[-1]
            else:
                ys_plot[ix] = self.ys[ix]

        mpl.fill(xs_plot, ys_plot, self.color, alpha=self.alpha)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
