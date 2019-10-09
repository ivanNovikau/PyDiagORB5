import Mix as mix
import Constants as cst
import Geom as geom
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cst)


class Curve:
    name = ''
    xs = None
    xs_err = None
    ys_err = None
    ys = None
    zs = None
    ws = None
    legend = "_"
    style = '-'
    width = None
    color = 'blue'
    markersize = None
    markerfacecolor = "None"
    colormap = 'jet'  # hot, jet, pink, hot_r, jet_r etc.
    levels = 10  # for contour plot
    pr_alpha = 1
    flag_errorbar = False

    flag_hist = False

    data_norm_to = None

    def XS(self, v):
        self.xs = v
        return self

    def XS_ERR(self, v):
        self.xs_err = v
        return self

    def YS(self, v):
        self.ys = v
        return self

    def YS_ERR(self, v):
        self.ys_err = v
        return self

    def ZS(self, v):
        self.zs = v
        return self

    def WS(self, v):
        self.ws = v
        return self

    def leg(self, v):
        if v is not None:
            if self.legend == "_":
                self.legend = ''

            if isinstance(v, list):
                self.legend += v[0]
                for i_line in range(1, np.shape(v)[0]):
                    self.legn(v[i_line])
            else:
                self.legend += v
        return self

    def legn(self, v):
        if v is not None:
            if self.legend == "_":
                self.legend = ''
            self.legend += '$\n \\boldmath $' + v
        return self

    def sty(self, v):
        self.style = v
        return self

    def new_sty(self, id0):
        def_styles = ['-', ':', '--', '-.']
        if id0>= len(def_styles):
            id0 = np.mod(id0, len(def_styles))
        self.style = def_styles[id0]
        return self

    def w(self, v):
        self.width = v
        return self

    def col(self, v):
        if v is not None:
            self.color = v
        return self

    def ms(self, v):
        self.markersize = v
        return self

    def mfc(self, v):
        self.markerfacecolor = v
        return self

    def name(self, v):
        self.name = v
        return self

    def cmp(self, v):
        self.colormap = v
        return self

    def lev(self, v):
        self.levels = v
        return self

    def set_hist(self):
        self.flag_hist = True
        return self

    def alpha(self, v):
        self.pr_alpha = v
        return self

    def norm_to(self, v):
        self.data_norm_to = v
        return self

    def set_errorbar(self, v, ys=None, xs=None):
        self.flag_errorbar = v
        self.ys_err = ys
        self.xs_err = xs
        return self


class Curves:
    list_curves = None
    map_curves = None
    n_curves = 0

    list_geoms = None
    n_geoms = 0

    list_text = None

    xlabel = None
    ylabel = None
    zlabel = None
    wlabel = None
    title = None

    flag_semilogy = False
    flag_norm = False

    axisFS = None
    lineW = None
    markerS = None
    fontS = None

    xlimits = None
    ylimits = None
    zlimits = None

    xticks = np.nan
    yticks = np.nan
    xticks_labels = np.nan
    yticks_labels = np.nan

    flag_legend = True
    legend_position = 'best'  # 'upper right', 'center left'
    legend_fcol = 'lightgray'

    def_colors = ['b', 'r', 'g', 'black', 'm', 'c',
                  'y', 'k', 'cyan', 'Purple', 'gray'
                  'lightcoral']
    flag_diff_styles = False
    def_styles = ['-', ':', '-.', ':']

    x_style = 'sci'  # 'sci', 'plain'
    y_style = 'sci'  # 'sci', 'plain'

    flag_maxlocator = False
    maxlocator = 6

    def __init__(self):
        self.list_curves = []
        self.map_curves = {}
        self.set_work()

        self.list_geoms = []
        self.list_text = []

    def new(self, name_curve=None):
        new_curve = Curve()
        self.list_curves.append(new_curve)
        self.n_curves += 1
        if name_curve is None:
            name_curve = 'curve_' + str(self.n_curves-1)
        self.map_curves[name_curve] = new_curve.name(name_curve)

        new_curve.w(self.lineW).ms(self.markerS)
        new_curve.col(self.new_color())

        if self.flag_diff_styles:
            new_curve.sty(self.new_style())

        return new_curve

    def newg(self, one_geom):
        self.list_geoms.append(one_geom)
        self.n_geoms += 1

    def newt(self, one_text):
        self.list_text.append(one_text)

    def load(self, curves_to_load):
        if curves_to_load is None:
            return
        for one_curve in curves_to_load.list_curves:
            self.n_curves += 1
            self.list_curves.append(one_curve)
            self.map_curves[one_curve.name] = one_curve

    def set_colors_styles(self):
        count_curve = -1
        for one_curve in self.list_curves:
            count_curve = count_curve + 1

            # color
            if count_curve + 1 <= len(self.def_colors):
                one_color = self.def_colors[count_curve]
            else:
                one_color = ','.join('{}'.format(*k) for k in enumerate(np.random.rand(3, 1)))
                one_color = 'rgb({})'.format(one_color)

            # style
            if count_curve + 1 <= len(self.def_styles):
                one_style = self.def_styles[count_curve]
            else:
                one_style = '-.'

            one_curve.col(one_color).sty(one_style)

    def new_color(self):
        if self.n_curves <= len(self.def_colors):
            one_color = self.def_colors[self.n_curves - 1]
        else:
            one_color = ','.join('{}'.format(*k) for k in enumerate(np.random.rand(3, 1)))
            one_color = 'rgb({})'.format(one_color)
        return one_color

    def new_style(self):
        if self.n_curves <= len(self.def_colors):
            one_style = self.def_styles[self.n_curves - 1]
        else:
            one_style = ':'
        return one_style

    def xlab(self, v):
        self.xlabel = v
        return self

    def ylab(self, v):
        self.ylabel = v
        return self

    def zlab(self, v):
        self.zlabel = v
        return self

    def wlab(self, v):
        self.wlabel = v
        return self

    def tit(self, v):
        if v is not None:
            if self.title is None:
                self.title = ''

            if isinstance(v, list):
                self.title += v[0]
                for i_line in range(1, np.shape(v)[0]):
                    self.titn(v[i_line])
            else:
                self.title += v

        return self

    def titn(self, v):
        if v is not None:
            if self.title is None:
                self.title = ''
            self.title += '$\n \\boldmath $' + v
        return self

    def new_tit(self, v):
        self.title = v
        return self

    def n(self):
        return self.n_curves

    def list(self, i_curve):
        return self.list_curves[i_curve]

    def map(self, name_curve):
        return self.map_curves.get(name_curve, None)

    def set_print(self):
        self.axisFS = 24
        self.lineW = 10
        self.markerS = 32
        self.fontS = 38

        for curve in self.list_curves:
            curve.w(self.lineW).ms(self.markerS)
        return self

    def set_work(self):
        self.axisFS = 22
        self.lineW = 6
        self.markerS = 14
        self.fontS = 22

        for curve in self.list_curves:
            curve.w(self.lineW).ms(self.markerS)
        return self

    def xlim(self, v):
        if v is not None:
            self.xlimits = v
        return self

    def ylim(self, v):
        if v is not None:
            self.ylimits = v
        return self

    def zlim(self, v):
        self.zlimits = v
        return self

    def set_limits(self):
        cr1 = self.list_curves[0]
        if cr1.xs is not None:
            self.xlimits = [cr1.xs[0], cr1.xs[-1]]
        if cr1.ys is not None:
            self.ylimits = [cr1.ys[0], cr1.ys[-1]]
        if cr1.zs is not None:
            self.zlimits = [cr1.zs[0], cr1.zs[-1]]

    def set_diff_styles(self):
        self.flag_diff_styles = True
        return self

    def leg_pos(self, v):
        self.legend_position = v
        return self

    def xt(self, v, lab=np.nan):
        self.xticks = v
        self.xticks_labels = lab
        return self

    def yt(self, v, lab=np.nan):
        self.yticks = v
        self.yticks_labels = lab
        return self

    def xsty(self, v):  # 'sci', 'plain'
        self.x_style = v
        return self

    def ysty(self, v):  # 'sci', 'plain'
        self.y_style = v
        return self

    def sort_legends(self):
        count_curve = -1
        new_list_curves = []
        list_curves_emp_leg = []
        for one_curve in self.list_curves:
            count_curve = count_curve + 1
            leg = one_curve.legend

            if leg == '\ ' or leg == ' ' or len(leg) == 0:
                list_curves_emp_leg.append(one_curve)
            else:
                new_list_curves.append(one_curve)
        new_list_curves.extend(list_curves_emp_leg)
        self.list_curves = new_list_curves

    def is_empty(self):
        if len(self.list_curves) == 0:
            return True
        else:
            return False


class PlText:
    x = None
    y = None
    line = ''
    color = 'black'

    def __init__(self, oo):
        self.init_from_oo(oo)

    def init_from_oo(self, oo):
        self.x = oo.get('x', None)
        self.y = oo.get('y', None)
        self.line = oo.get('line', '')
        self.color = oo.get('color', 'black')




