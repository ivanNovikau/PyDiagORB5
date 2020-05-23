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


class PageCurves(ivis_base_page.BasePage):

    def __init__(self, **kwargs):
        super(PageCurves, self).__init__(**kwargs)

        self.create_elements()

    def create_elements(self):
        curves = self.mw.curves

        # row counter
        cnt = mix.Counter()

        # # labels:
        # ivb.BLabel(master=cf, text="List of curves: ").grid(row=0, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        # ivb.BLabel(master=cf, text="Color: ").grid(row=1, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        # ivb.BLabel(master=cf, text="Colormap: ").grid(row=2, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        # ivb.BLabel(master=cf, text="Styles: ").grid(row=3, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        # ivb.BLabel(master=cf, text="Width: ").grid(row=4, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        # ivb.BLabel(master=cf, text="Marker size: ").grid(row=5, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        #
        # cf.columnconfigure(0, weight=1)
        # n_rows = 6
        # for id_row in range(n_rows):
        #     cf.rowconfigure(id_row, weight=1)