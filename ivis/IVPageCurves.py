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
    id_selected_legend = None

    def __init__(self, **kwargs):
        super(PageCurves, self).__init__(**kwargs)

        self.create_elements()

    def create_elements(self):
        curves = self.mw.curves

        # row counter
        cnt = mix.Counter()

        # create curves' dictionaries
        self.elements['curve'] = {}

        # OptionMenu with available curves
        self.elements['curve']['om'] = ivb.LabelledOptionMenu(
            self.frame,
            "Curves: ",
            [cnt.next(), 0],
            curves.get_legends(),
            self.get_selected_curve
        )

        # *** set relative sizes ***
        self.frame.columnconfigure(0, weight=1)
        self.frame.columnconfigure(1, weight=2)
        for id_row in range(cnt.n_elements):
            self.frame.rowconfigure(id_row, weight=0)

        # --- Frame with properties of a selected additional text ---
        self.create_frame_curve_properties(cnt.next())

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

    def update_elements(self):
        self.update_curves_om()

    def update_curves_om(self):
        curves = self.mw.curves
        self.elements['curve']['om'].update_options(
            curves.get_legends()
        )

    def create_frame_curve_properties(self, id_row):
        fcrv = self.elements['curve']

        # create a frame
        fcrv['frame'] = tk.Frame(
            self.frame,
            bg=GLO.IVIS_selected_element_inf_frame
        )
        fcrv['frame'].grid(
            row=id_row, column=0, columnspan=2,
            sticky=tk.N + tk.S + tk.E + tk.W,
        )

        # row counter
        cnt = mix.Counter()

        # legend
        fcrv['legend'] = ivb.LabelledEntry(
            fcrv['frame'], "Legend: ", [cnt.next(), 0], ''
        ).var
        fcrv['legend'].trace(
            'w',
            self.curve_write_legend
        )

    def get_selected_curve(self, opt_selected, id_selected):
        self.id_selected_legend = id_selected
        fcrv = self.elements['curve']
        if self.id_selected_legend is not None:
            ocurve = self.mw.curves.list_curves[self.id_selected_legend]

            fcrv['legend'].set(ocurve.ff['legend'])

            fcrv['frame'].tkraise()
        else:
            pass

    def curve_write_legend(self, *args):
        if self.id_selected_legend is not None:
            res = self.elements['curve']['legend'].get()
            ocurve = self.mw.curves.list_curves[self.id_selected_legend]
            ocurve.ff['legend'] = res




