import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr
import ivis.IVButtons as ivb
import ivis.IVMenus as ivm
import ivis.BasePage as ivis_base_page

import numpy as np
import re
import matplotlib.lines as mlines
import matplotlib.patches as patches

import tkinter as tk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

from matplotlib.backend_bases import key_press_handler, button_press_handler
from matplotlib.ticker import AutoLocator


def reload():
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)
    mix.reload_module(GLO)
    mix.reload_module(cpr)
    mix.reload_module(ivb)
    mix.reload_module(ivm)
    mix.reload_module(ivis_base_page)


class PageAxes(ivis_base_page.BasePage):
    id_current_text = None

    def __init__(self, **kwargs):
        super(PageAxes, self).__init__(**kwargs)

        self.create_elements()

    def create_elements(self):
        # row counter
        cnt = mix.Counter()

        # Buttons
        self.elements['bDefault'] = ivb.BButton(
            master=self.frame,
            text="Default axis",
            command=self.default_axis
        )
        self.elements['bDefault'].grid(row=cnt.next(), column=0)

    def default_axis(self):
        # --- create a figure and plot data ---
        mw = self.mw
        ax = mw.get_ax()

        mw.curves = crv.copy_curves(mw.curves_default, ax)

        cpr.plot(
            mw.curves_default, mw.fig, ax,
            FIG_W=mw.curves.ff['figure_width'] / 2,
            FIG_H=mw.curves.ff['figure_height'] / 2,
        )
        mw.draw()

        # update elements:
        mw.fLeft.sections['ax']['pages']['text'].update_default_elements()
