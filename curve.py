import Mix as mix
import Constants as cst
import Geom as geom
import Global_variables as GLO


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cst)
    mix.reload_module(GLO)


class Curve:
    CurveName = ''
    ff = None  # format of the curve
    xs = None
    xs_err = None
    ys_err = None
    ys = None
    zs = None
    ws = None
    data_norm_to = None

    def __init__(self):
        self.ff = dict(GLO.DEF_CURVE_FORMAT)

    def name(self, v):
        self.CurveName = v
        return self

    def XS(self, v):
        self.xs = v
        return self

    def YS(self, v):
        self.ys = v
        return self

    def ZS(self, v):
        self.zs = v
        return self

    def WS(self, v):
        self.ws = v
        return self

    def set_ff(self, v):
        self.ff = dict(v)
        return self

    def norm_to(self, v):
        self.data_norm_to = v
        return self

    def set_errorbar(self, v, ys=None, xs=None):
        self.ff['flag_errorbar'] = v
        self.ys_err = ys
        self.xs_err = xs
        return self


class Curves:
    list_curves = None
    map_curves = None
    n_curves = 0
    list_text = None
    list_geoms = None
    n_geoms = 0
    ff = None  # format

    def __init__(self):
        self.list_curves = []
        self.map_curves = {}

        self.list_geoms = []
        self.list_text = []

        # create default format:
        self.ff = GLO.DEF_PLOT_FORMAT

    def new(self, name_curve=None):
        new_curve = Curve()
        self.list_curves.append(new_curve)
        self.n_curves += 1
        if name_curve is None:
            name_curve = 'curve_' + str(self.n_curves-1)
        self.map_curves[name_curve] = new_curve.name(name_curve)

        return new_curve

    def load(self, curves_to_load):
        if curves_to_load is None:
            return
        for one_curve in curves_to_load.list_curves:
            self.n_curves += 1
            self.list_curves.append(one_curve)
            self.map_curves[one_curve.name] = one_curve

    def set_ff(self, v):
        self.ff = dict(v)
        return self

    def newg(self, one_geom):
        self.list_geoms.append(one_geom)
        self.n_geoms += 1

    def newt(self, one_text):
        self.list_text.append(one_text)

    def n(self):
        return self.n_curves

    def set_fixed_limits(self):
        cr1 = self.list_curves[0]
        if cr1.xs is not None:
            self.ff['xlimits'] = [cr1.xs[0], cr1.xs[-1]]
        if cr1.ys is not None:
            self.ff['ylimits'] = [cr1.ys[0], cr1.ys[-1]]
        if cr1.zs is not None:
            self.ff['zlimits'] = [cr1.zs[0], cr1.zs[-1]]

    def sort_legends(self):
        count_curve = -1
        new_list_curves = []
        list_curves_emp_leg = []
        for one_curve in self.list_curves:
            count_curve = count_curve + 1
            leg = one_curve.ff['legend']

            if leg is None:
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
    coef_width = 1

    def __init__(self, oo):
        self.init_from_oo(oo)

    def init_from_oo(self, oo):
        self.x = oo.get('x', None)  # in units of x-axis
        self.y = oo.get('y', None)  # in units of y-axis
        self.color = oo.get('color', 'black')
        self.coef_width = oo.get('coef_width', 1)
        self.line = mix.create_line_from_list(oo.get('line', ''))

        # self.line = line_res if line_res != '' else '\ '




