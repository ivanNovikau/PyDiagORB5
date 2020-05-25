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
        for id_column in range(2):
            self.frame.columnconfigure(id_column, weight=1)
            for id_row in range(cnt.n_elements):
                self.frame.rowconfigure(id_row, weight=0)

        # --- Frame with properties of a selected additional text ---
        self.create_frame_curve_properties(cnt.next())

        # --- Create comparison table ---
        rows_list = ['X', 'Y', 'Z']
        columns_list = self.mw.curves.get_legends()

        self.elements['table-comparison'] = ivb.Table(
            'Curves Comparison Table', self.frame, (cnt.next(), 0, 1, 2),
            rows_list, columns_list,
            default_rows_columns=(0, 0, 3, 0)
        )

        # ivb.BLabel(master=cf, text="Colormap: ").grid(row=2, column=0, sticky=tk.N + tk.S + tk.E + tk.W)

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

        # XY data
        fcrv['x'] = ivb.LabelledEntry(
            fcrv['frame'], "X: ", [cnt.next(), 0], ''
        ).var
        fcrv['y'] = ivb.LabelledEntry(
            fcrv['frame'], "Y: ", [cnt.next(), 0], ''
        ).var
        fcrv['z'] = ivb.LabelledEntry(
            fcrv['frame'], "Z: ", [cnt.next(), 0], ''
        ).var

        # add an checkbutton "on fly":
        fcrv['flag_on_fly'] = tk.IntVar()
        tk.Checkbutton(
            fcrv['frame'],
            text="See XYZ on fly",
            variable=fcrv['flag_on_fly']
        ).grid(row=cnt.next(), column=0)

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

    def follow_curve_xyz_data(self, id_curve, xcurve, ycurve, zcurve=None):
        fcrv = self.elements['curve']

        # follow XYZ data of a selected curve
        if self.id_selected_legend == id_curve:
            if fcrv['flag_on_fly'].get():
                fcrv['x'].set(xcurve)
                fcrv['y'].set(ycurve)
                fcrv['z'].set(zcurve)

        # follow XYZ data in the table:
        if fcrv['flag_on_fly'].get():
            self.elements['table-comparison'].set_cell(0, id_curve, xcurve)  # set x
            self.elements['table-comparison'].set_cell(1, id_curve, ycurve)  # set y
            self.elements['table-comparison'].set_cell(2, id_curve, zcurve)  # set z

    def set_curve_xyz_point(self, xdata, ydata):
        if xdata is None or ydata is None:
            return

        fcrv = self.elements['curve']

        # set XYZ data of a selected curve:
        if self.id_selected_legend is not None:
            if not fcrv['flag_on_fly'].get():
                one_curve = self.mw.curves.list_curves[self.id_selected_legend]

                id_x, xcurve, _ = mix.get_ids(one_curve.xs, xdata)
                if not one_curve.get_flag_2d():
                    ycurve = one_curve.ys[id_x]
                else:
                    id_y, ycurve, _ = mix.get_ids(one_curve.ys, ydata)
                    zcurve = one_curve.zs[id_x, id_y]
                    fcrv['z'].set(zcurve)

                fcrv['x'].set(xcurve)
                fcrv['y'].set(ycurve)

        # set XYZ data in the table:
        if not fcrv['flag_on_fly'].get():
            tcomp = self.elements['table-comparison']
            for id_curve in tcomp.map_column_ids.keys():
                one_curve = self.mw.curves.list_curves[id_curve]
                id_x, xcurve, _ = mix.get_ids(one_curve.xs, xdata)
                if not one_curve.get_flag_2d():
                    ycurve = one_curve.ys[id_x]
                else:
                    id_y, ycurve, _ = mix.get_ids(one_curve.ys, ydata)
                    zcurve = one_curve.zs[id_x, id_y]
                    tcomp.set_cell(2, id_curve, zcurve)  # set z
                tcomp.set_cell(0, id_curve, xcurve)  # set x
                tcomp.set_cell(1, id_curve, ycurve)  # set y
