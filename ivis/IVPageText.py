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


class PageText(ivis_base_page.BasePage):
    id_current_text = None

    def __init__(self, **kwargs):
        super(PageText, self).__init__(**kwargs)

        self.create_elements()

    def update_elements(self):
        curves = self.mw.curves
        self.elements['title'].set(curves.ff['title'] if curves.ff['title'] is not None else "")
        self.elements['xlabel'].set(curves.ff['xlabel'] if curves.ff['xlabel'] is not None else "")
        self.elements['ylabel'].set(curves.ff['ylabel'] if curves.ff['ylabel'] is not None else "")
        self.elements['xaticks'].set("")
        self.elements['yaticks'].set("")
        self.update_atext_om()

    def update_atext_om(self):
        curves = self.mw.curves
        self.elements['atext']['om'].update_options(
            mix.get_atext_from_curves(curves)
        )

    def create_elements(self):
        curves = self.mw.curves

        # row counter
        cnt = mix.Counter()

        # --- entries ---
        self.elements['title'] = ivb.LabelledEntry(
            self.frame, "Title: ", [cnt.next(), 0], curves.ff['title']
        ).var
        self.elements['title'].trace('w', self.write_title)

        self.elements['xlabel'] = ivb.LabelledEntry(
            self.frame, "X label: ", [cnt.next(), 0], curves.ff['xlabel']
        ).var
        self.elements['xlabel'].trace('w', self.write_xlabel)

        self.elements['ylabel'] = ivb.LabelledEntry(
            self.frame, "Y label: ", [cnt.next(), 0], curves.ff['ylabel']
        ).var
        self.elements['ylabel'].trace('w', self.write_ylabel)

        self.elements['xaticks'] = ivb.LabelledEntry(
            self.frame, "X additional ticks: ", [cnt.next(), 0], ""
        ).var
        self.elements['xaticks'].trace('w', self.write_xaticks)

        self.elements['yaticks'] = ivb.LabelledEntry(
            self.frame, "Y additional ticks: ", [cnt.next(), 0], ""
        ).var
        self.elements['yaticks'].trace('w', self.write_yaticks)

        # --- Create a dictionary to save elements to describe additional text ---
        self.elements['atext'] = {}

        # --- OptionMenu with additional texts ---
        self.elements['atext']['om'] = ivb.LabelledOptionMenu(
            self.frame,
            "Additional text: ",
            [cnt.next(), 0],
            mix.get_atext_from_curves(curves),
            self.get_selected_text
        )

        # *** set relative sizes ***
        self.frame.columnconfigure(0, weight=1)
        self.frame.columnconfigure(1, weight=2)
        for id_row in range(cnt.n_elements):
            self.frame.rowconfigure(id_row, weight=0)

        # --- Frame with properties of a selected additional text ---
        self.create_frame_additional_text(cnt.next())

    def create_frame_additional_text(self, id_row):
        # create a frame
        self.elements['atext']['frame'] = tk.Frame(
            self.frame,
            bg=GLO.IVIS_selected_element_inf_frame
        )
        self.elements['atext']['frame'].grid(
            row=id_row, column=0, columnspan=2,
            sticky=tk.N + tk.S + tk.E + tk.W,
        )

        # add an entry with text
        self.elements['atext']['text'] = ivb.LabelledEntry(
            self.elements['atext']['frame'], "Text: ", [0, 0], ''
        ).var
        self.elements['atext']['text'].trace(
            'w',
            self.atext_write_text
        )

        # add an entry with text x coodinate
        self.elements['atext']['x'] = ivb.LabelledEntry(
            self.elements['atext']['frame'], "X: ", [1, 0], ''
        ).var
        self.elements['atext']['x'].trace(
            'w',
            self.atext_write_x
        )

        # add an entry with text y coodinate
        self.elements['atext']['y'] = ivb.LabelledEntry(
            self.elements['atext']['frame'], "Y: ", [2, 0], ''
        ).var
        self.elements['atext']['y'].trace(
            'w',
            self.atext_write_y
        )

        # add an entry with text color
        self.elements['atext']['color'] = ivb.LabelledOptionMenu(
            self.elements['atext']['frame'],
            "Color: ",
            [3, 0],
            ["black", "red", "green", "blue", "grey"],
            None
        )
        self.elements['atext']['color'].var.trace(
            'w',
            self.atext_write_color
        )

        # add an checkbutton invisible:
        self.elements['atext']['flag_invisible'] = tk.IntVar()
        tk.Checkbutton(
            self.elements['atext']['frame'],
            text="Invisible",
            variable=self.elements['atext']['flag_invisible']
        ).grid(row=4, column=0)
        self.elements['atext']['flag_invisible'].trace(
            'w',
            self.atext_set_invisible
        )

    def get_selected_text(self, opt_selected, id_selected):
        self.id_current_text = id_selected
        if self.id_current_text is not None:
            el_text = self.mw.curves.list_text[id_selected]

            self.elements['atext']['text'].set(
                mix.delete_bold_keyword(el_text.line)
            )
            self.elements['atext']['x'].set(el_text.x)
            self.elements['atext']['y'].set(el_text.y)
            self.elements['atext']['color'].var.set(el_text.color)
            self.elements['atext']['flag_invisible'].set(el_text.flag_invisible)

            self.elements['atext']['frame'].tkraise()
        else:
            self.elements['atext']['text'].set('')
            self.elements['atext']['x'].set('')
            self.elements['atext']['y'].set('')
            self.elements['atext']['color'].var.set('black')
            self.elements['atext']['flag_invisible'].set(0)

    def atext_write_text(self, *args):
        if self.id_current_text is not None:
            res_line = self.elements['atext']['text'].get()
            res_line = mix.create_line_from_list(res_line)
            self.mw.curves.list_text[self.id_current_text].line = res_line

    def atext_write_x(self, *args):
        if self.id_current_text is not None:
            res_x = float(self.elements['atext']['x'].get())
            self.mw.curves.list_text[self.id_current_text].x = res_x

    def atext_write_y(self, *args):
        if self.id_current_text is not None:
            res_y = float(self.elements['atext']['y'].get())
            self.mw.curves.list_text[self.id_current_text].y = res_y

    def atext_write_color(self, *args):
        if self.id_current_text is not None:
            res_color = self.elements['atext']['color'].var.get()
            self.mw.curves.list_text[self.id_current_text].color = res_color

    def atext_set_invisible(self, *args):
        if self.id_current_text is not None:
            self.mw.curves.list_text[self.id_current_text].flag_invisible = \
                self.elements['atext']['flag_invisible'].get()

    def write_title(self, *args):
        self.mw.curves.ff['title'] = self.elements['title'].get()

    def write_xlabel(self, *args):
        self.mw.curves.ff['xlabel'] = self.elements['xlabel'].get()

    def write_ylabel(self, *args):
        self.mw.curves.ff['ylabel'] = self.elements['ylabel'].get()

    def write_xaticks(self, *args):
        ivis_add_xticks = list(
            mix.array_from_str(self.elements['xaticks'].get())
        )
        self.mw.curves.ff['ivis_add_xticks'] = ivis_add_xticks
        if len(ivis_add_xticks) > 0:
            self.mw.curves.ff['xticks'] = \
                self.mw.curves_default.ff['xticks'] + ivis_add_xticks

    def write_yaticks(self, *args):
        ivis_add_yticks = list(
            mix.array_from_str(self.elements['yaticks'].get())
        )
        self.mw.curves.ff['ivis_add_yticks'] = ivis_add_yticks
        if len(ivis_add_yticks) > 0:
            self.mw.curves.ff['yticks'] = \
                self.mw.curves_default.ff['yticks'] + ivis_add_yticks




