import Mix as mix
import Constants as cst
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cst)


class Curve:
    name = ''
    xs = None
    ys = None
    zs = None
    xlabel = ''
    ylabel = ''
    legend = ''
    style = '-'
    width = None
    color = 'blue'
    markersize = None

    def XS(seld, v):
        seld.xs = v
        return seld

    def YS(self, v):
        self.ys = v
        return self

    def ZS(self, v):
        self.zs = v
        return self

    def xlab(self, v):
        self.xlabel = v
        return self

    def ylab(self, v):
        self.ylabel = v
        return self

    def leg(self, v):
        self.legend = v
        return self

    def sty(self, v):
        self.style = v
        return self

    def w(self, v):
        self.width = v
        return self

    def col(self, v):
        self.color = v
        return self

    def ms(self, v):
        self.markersize = v
        return self

    def name(self, v):
        self.name = v
        return self


class Curves:
    list_curves = None
    map_curves = None
    n_curves = 0

    xlabel = ''
    ylabel = ''

    flag_semilogy = False

    axisFS = None
    lineW = None
    markerS = None
    fontS = None

    def __init__(self):
        self.list_curves = []
        self.map_curves = {}
        self.set_work()

    def new(self, name_curve=None):
        new_curve = Curve()
        self.list_curves.append(new_curve)
        self.n_curves += 1
        if name_curve is None:
            name_curve = 'curve_' + str(self.n_curves-1)
        self.map_curves[name_curve] = new_curve.name(name_curve)

        new_curve.w(self.lineW).ms(self.markerS)

        return new_curve

    def xlab(self, v):
        self.xlabel = v
        return self

    def ylab(self, v):
        self.ylabel = v
        return self

    def n(self):
        return self.n_curves

    def list(self, i_curve):
        return self.list_curves[i_curve]

    def map(self, name_curve):
        return self.map_curves.get(name_curve, None)

    def set_print(self):
        self.axisFS = 48
        self.lineW = 10
        self.markerS = 32
        self.fontS = 48

        for curve in self.list_curves:
            curve.w(self.lineW).ms(self.markerS)

    def set_work(self):
        self.axisFS = 10
        self.lineW = 4
        self.markerS = 12
        self.fontS = 10

        for curve in self.list_curves:
            curve.w(self.lineW).ms(self.markerS)