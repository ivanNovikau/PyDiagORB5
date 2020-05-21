import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr
import ivis.IVButtons as ivb
import ivis.IVMenus as ivm

import numpy as np

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
        cf = self.tabController.create_page('Text')
        curves = self.mw.curves

        # Entries:
        eTitle = ivb.LabelledEntry(cf, "Title: ", [0, 0], curves.ff['title']).var
        eXlabel = ivb.LabelledEntry(cf, "X label: ", [1, 0], curves.ff['xlabel']).var
        eYlabel = ivb.LabelledEntry(cf, "Y label: ", [2, 0], curves.ff['ylabel']).var
        eXAticks = ivb.LabelledEntry(cf, "X additional ticks: ", [3, 0], "").var
        eYAticks = ivb.LabelledEntry(cf, "Y additional ticks: ", [4, 0], "").var
        # ivb.BLabel(master=cf, text="Additional texts: ").grid(row=5, column=0, sticky=tk.N + tk.S + tk.E + tk.W)

        cf.columnconfigure(0, weight=1)
        cf.columnconfigure(1, weight=2)
        n_rows = 6
        for id_row in range(n_rows):
            cf.rowconfigure(id_row, weight=0)

        # Buttons to rebuild or reset the plot:
        def rebuild_plot():
            ivis_add_xticks = list(mix.array_from_str(eXAticks.get()))
            ivis_add_yticks = list(mix.array_from_str(eYAticks.get()))
            curves.ff.update({
                'title': eTitle.get(),
                'xlabel': eXlabel.get(),
                'ylabel': eYlabel.get(),
                'xticks': self.mw.ff_default['xticks'] + ivis_add_xticks,
                'yticks': self.mw.ff_default['yticks'] + ivis_add_yticks,
                'ivis_add_xticks': ivis_add_xticks,
                'ivis_add_yticks': ivis_add_yticks,
            })

            ax = self.mw.fig.axes[0]

            cpr.format_plot(
                self.mw.fig, ax, curves, self.mw.flag_2d
            )
            self.mw.fig.canvas.draw()

        def default_plot():
            curves.ff = dict(self.mw.ff_default)

            eTitle.set(curves.ff['title'] if curves.ff['title'] is not None else "")
            eXlabel.set(curves.ff['xlabel'] if curves.ff['xlabel'] is not None else "")
            eYlabel.set(curves.ff['ylabel'] if curves.ff['ylabel'] is not None else "")
            eXAticks.set("")
            eYAticks.set("")

            ax = self.mw.fig.axes[0]

            cpr.format_plot(
                self.mw.fig, ax, curves, self.mw.flag_2d
            )
            self.mw.fig.canvas.draw()

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






