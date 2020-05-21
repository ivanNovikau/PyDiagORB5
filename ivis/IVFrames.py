import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr
import ivis.IVButtons as ivb
import ivis.IVMenus as ivm

import numpy as np
import re

import tkinter as tk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

from matplotlib.backend_bases import key_press_handler, button_press_handler


def reload():
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)
    mix.reload_module(GLO)
    mix.reload_module(cpr)
    mix.reload_module(ivb)
    mix.reload_module(ivm)


# *** Basic Frame ***
class BFrame(tk.Frame):
    mw = None  # main window

    def __init__(self, mw, **kwargs):
        if 'bg' not in kwargs:
            kwargs['bg'] = GLO.IVIS_frame_color
        super(BFrame, self).__init__(**kwargs)
        self.mw = mw


# *** Figure Frame (where figures are plotted) ***
class FigureFrame(BFrame):
    # frames
    fUpper = None
    fBottom = None

    # popup that appear then right button is pressed
    popup_menu_canvas = None

    # canvas where figure is plotted
    fig_canvas = None

    def __init__(self, mw, **kwargs):
        super(FigureFrame, self).__init__(mw, **kwargs)

        # upper figure frame
        self.fUpper = tk.Frame(
            master=self,
            height=GLO.IVIS_height_tab_frame,
            bg=GLO.IVIS_color_tabs_frame,
        )

        # create canvas
        self.fig_canvas = FigureCanvasTkAgg(self.mw.fig, master=self)
        self.fig_canvas.draw()

        # bottom figure frame
        self.fBottom = BottomFigureFrame(
            mw=self.mw,
            figure_frame=self,
            master=self,
            height=40
        )

        # figure style
        self.mw.fig.patch.set_facecolor(GLO.IVIS_canvas_color)

        # press events
        self.fig_canvas.mpl_connect("key_press_event", self.on_key_press)
        self.fig_canvas.mpl_connect("button_press_event", self.on_button_press)

        # create the popup menu of a canvas:
        self.popup_menu_canvas = ivm.PopupCanvasMenu(
            mw=self.mw,
            master=self.fig_canvas.get_tk_widget(),
            tearoff=0
        )

        # set position of frame elemens
        self.arrange_elements()

    def arrange_elements(self):
        self.fUpper.grid(row=0, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.fig_canvas.get_tk_widget().grid(
            row=1, column=0, sticky=tk.N + tk.S + tk.E + tk.W
        )
        self.fBottom.grid(row=2, column=0, sticky=tk.N + tk.S + tk.E + tk.W)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=0)

    def on_key_press(self, event):
        print("key: you pressed {}".format(event.key))
        key_press_handler(event, self.fig_canvas)

    def on_button_press(self, event):
        if event.button == 3:
            self.popup_menu_canvas.call(event)

        button_press_handler(event, self.fig_canvas)


# *** Bottom Figure Frame (with some properties) ***
class BottomFigureFrame(BFrame):
    fFigure = None  # figure frame

    # label with axes coordinates
    lbXYdata = None

    def __init__(self, mw, figure_frame, **kwargs):
        super(BottomFigureFrame, self).__init__(mw, **kwargs)
        self.fFigure = figure_frame

        # Add labels
        self.lbXYdata = tk.Label(
            master=self,
            text='XY data',
            bg=GLO.IVIS_label_color
        )

        # set position of frame elements
        self.lbXYdata.grid(row=0, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)


# *** Left Frame (with Figure, Axes properties and post-processing) ***
class LeftFrame(BFrame):
    # frames
    fFigProp = None
    fAxProp = None
    fProc = None

    def __init__(self, mw, **kwargs):
        super(LeftFrame, self).__init__(mw, **kwargs)

        # Add frames:
        self.fFigProp = FigPropFrame(self.mw, master=self)
        self.fAxProp = AxPropFrame(self.mw, master=self)
        self.fProc = ProcessingFrame(self.mw, master=self)

        # set position of elements
        self.fFigProp.grid(row=0, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.fAxProp.grid(row=1, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.fProc.grid(row=2, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)


# *** Frame with Figure properties  ***
class FigPropFrame(BFrame):
    def __init__(self, mw, **kwargs):
        super(FigPropFrame, self).__init__(mw, **kwargs)

        # style
        self.config(
            highlightbackground=GLO.IVIS_border_color,
            highlightthickness=2,
        )


# *** Frame with Axes properties  ***
class AxPropFrame(BFrame):
    tabController = None
    n_columns = 5
    omAText = None

    # frames with information about different elements
    fElementInf = None  # root element frame
    feiText = None
    feiText_elements = {}
    id_current_text = None

    def __init__(self, mw, **kwargs):
        super(AxPropFrame, self).__init__(mw, **kwargs)

        # style
        self.config(
            highlightbackground=GLO.IVIS_border_color,
            highlightthickness=2,
        )

        # Create Tab Controller
        self.tabController = ivb.TabController(self)

        # Create pages
        self.create_frame_text()
        self.create_frame_curves()

        # Activate Text tab:
        self.tabController.call_page('Text')

        # self.fText.tkraise()

    def create_frame_text(self):
        # --- FUNCTIONS ------------------------------------------------------------------
        def rebuild_plot():
            ivis_add_xticks = list(mix.array_from_str(eXAticks.get()))
            ivis_add_yticks = list(mix.array_from_str(eYAticks.get()))
            self.mw.curves.ff.update({
                'title': eTitle.get(),
                'xlabel': eXlabel.get(),
                'ylabel': eYlabel.get(),
                'xticks': self.mw.curves_default.ff['xticks'] + ivis_add_xticks,
                'yticks': self.mw.curves_default.ff['yticks'] + ivis_add_yticks,
                'ivis_add_xticks': ivis_add_xticks,
                'ivis_add_yticks': ivis_add_yticks,
            })

            # update the plot format
            ax = self.mw.fig.axes[0]
            ax.texts = []
            cpr.format_plot(
                self.mw.fig, ax, self.mw.curves, self.mw.flag_2d
            )
            self.mw.fig.canvas.draw()

        def default_plot():
            ax = self.mw.fig.axes[0]
            self.mw.curves = crv.copy_curves(self.mw.curves_default, ax)
            curves = self.mw.curves

            eTitle.set(curves.ff['title'] if curves.ff['title'] is not None else "")
            eXlabel.set(curves.ff['xlabel'] if curves.ff['xlabel'] is not None else "")
            eYlabel.set(curves.ff['ylabel'] if curves.ff['ylabel'] is not None else "")
            eXAticks.set("")
            eYAticks.set("")
            self.omAText.update_options(
                mix.get_atext_from_curves(curves)
            )

            ax.texts = []
            cpr.format_plot(
                self.mw.fig, ax, curves, self.mw.flag_2d
            )
            self.mw.fig.canvas.draw()

        def text_selected(opt_selected, id_selected):
            if id_selected is not None:
                self.id_current_text = id_selected

                el_text = self.mw.curves.list_text[id_selected]

                self.feiText_elements['text'].set(
                    mix.delete_bold_keyword(el_text.line)
                )
                self.feiText_elements['x'].set(el_text.x)
                self.feiText_elements['y'].set(el_text.y)
                self.feiText_elements['color'].var.set(el_text.color)
                self.feiText_elements['flag_invisible'].set(el_text.flag_invisible)

                self.feiText.tkraise()

        def write_text(*args):
            if self.id_current_text is not None:
                res_line = self.feiText_elements['text'].get()
                res_line = mix.create_line_from_list(res_line)
                self.mw.curves.list_text[self.id_current_text].line = res_line

        def write_text_x(*args):
            if self.id_current_text is not None:
                res_x = float(self.feiText_elements['x'].get())
                self.mw.curves.list_text[self.id_current_text].x = res_x

        def write_text_y(*args):
            if self.id_current_text is not None:
                res_y = float(self.feiText_elements['y'].get())
                self.mw.curves.list_text[self.id_current_text].y = res_y

        def write_text_color(*args):
            if self.id_current_text is not None:
                res_color = self.feiText_elements['color'].var.get()
                self.mw.curves.list_text[self.id_current_text].color = res_color

        def write_text_invisible(*args):
            if self.id_current_text is not None:
                self.mw.curves.list_text[self.id_current_text].flag_invisible = \
                    self.feiText_elements['flag_invisible'].get()

        # ---------------------------------------------------------------
        cf = self.tabController.create_page('Text')
        curves = self.mw.curves

        # Entries:
        eTitle = ivb.LabelledEntry(cf, "Title: ", [0, 0], curves.ff['title']).var
        eXlabel = ivb.LabelledEntry(cf, "X label: ", [1, 0], curves.ff['xlabel']).var
        eYlabel = ivb.LabelledEntry(cf, "Y label: ", [2, 0], curves.ff['ylabel']).var
        eXAticks = ivb.LabelledEntry(cf, "X additional ticks: ", [3, 0], "").var
        eYAticks = ivb.LabelledEntry(cf, "Y additional ticks: ", [4, 0], "").var

        self.omAText = ivb.LabelledOptionMenu(
            cf,
            "Additional text: ",
            [5, 0],
            mix.get_atext_from_curves(curves),
            text_selected
        )

        cf.columnconfigure(0, weight=1)
        cf.columnconfigure(1, weight=2)
        n_rows = 6
        for id_row in range(n_rows):
            cf.rowconfigure(id_row, weight=0)

        # --- FRAME with information from a selected element ---
        # - Root frame -
        self.fElementInf = tk.Frame(
            cf,
            bg=GLO.IVIS_selected_element_inf_frame
        )

        n_rows += 1
        self.fElementInf.grid(
            row=n_rows-1,
            column=0, columnspan=2,
            sticky=tk.N + tk.S + tk.E + tk.W,
        )

        # - Text frame -
        self.feiText = tk.Frame(
            self.fElementInf,
            bg=GLO.IVIS_selected_element_inf_frame
        )
        self.feiText.grid(
            row=0, column=0, sticky=tk.N + tk.S + tk.E + tk.W,
        )

        # entry with text
        self.feiText_elements['text'] = ivb.LabelledEntry(
            self.feiText, "Text: ", [0, 0], ''
        ).var
        self.feiText_elements['text'].trace(
            'w',
            write_text
        )

        # entry with text x coodinate
        self.feiText_elements['x'] = ivb.LabelledEntry(
            self.feiText, "X: ", [1, 0], ''
        ).var
        self.feiText_elements['x'].trace(
            'w',
            write_text_x
        )

        # entry with text y coodinate
        self.feiText_elements['y'] = ivb.LabelledEntry(
            self.feiText, "Y: ", [2, 0], ''
        ).var
        self.feiText_elements['y'].trace(
            'w',
            write_text_y
        )

        # entry with text color
        self.feiText_elements['color'] = ivb.LabelledOptionMenu(
            self.feiText,
            "Color: ",
            [3, 0],
            ["black", "red", "green", "blue", "grey"],
            None
        )
        self.feiText_elements['color'].var.trace(
            'w',
            write_text_color
        )

        # checkbutton invisible:
        self.feiText_elements['flag_invisible'] = tk.IntVar()
        tk.Checkbutton(
            self.feiText,
            text="Invisible",
            variable=self.feiText_elements['flag_invisible']
        ).grid(row=4, column=0)
        self.feiText_elements['flag_invisible'].trace(
            'w',
            write_text_invisible
        )

        # --- Buttons ---
        n_rows += 1
        ivb.BButton(
            master=cf,
            text='Rebuild plot',
            command=rebuild_plot
        ).grid(row=n_rows+1, column=0)
        ivb.BButton(
            master=cf,
            text='Default',
            command=default_plot
        ).grid(row=n_rows + 1, column=1)

        cf.rowconfigure(6, weight=1)

    def create_frame_curves(self):
        cf = self.tabController.create_page('Curves')

        # labels:
        ivb.BLabel(master=cf, text="List of curves: ").grid(row=0, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        ivb.BLabel(master=cf, text="Color: ").grid(row=1, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        ivb.BLabel(master=cf, text="Colormap: ").grid(row=2, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        ivb.BLabel(master=cf, text="Styles: ").grid(row=3, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        ivb.BLabel(master=cf, text="Width: ").grid(row=4, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        ivb.BLabel(master=cf, text="Marker size: ").grid(row=5, column=0, sticky=tk.N + tk.S + tk.E + tk.W)

        cf.columnconfigure(0, weight=1)
        n_rows = 6
        for id_row in range(n_rows):
            cf.rowconfigure(id_row, weight=1)


# *** Frame to work with post-processing  ***
class ProcessingFrame(BFrame):
    def __init__(self, mw, **kwargs):
        super(ProcessingFrame, self).__init__(mw, **kwargs)

        # style
        self.config(
            highlightbackground=GLO.IVIS_border_color,
            highlightthickness=2,
        )






