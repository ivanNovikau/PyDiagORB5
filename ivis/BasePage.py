import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr
import ivis.IVButtons as ivb
import ivis.IVMenus as ivm


def reload():
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)
    mix.reload_module(GLO)
    mix.reload_module(cpr)
    mix.reload_module(ivb)
    mix.reload_module(ivm)


class BasePage:
    frame = None
    mw = None
    elements = {}

    def __init__(self, frame, mw):
        self.frame = frame
        self.mw = mw

    def update_elements(self):
        curves = self.mw.curves
        pass
