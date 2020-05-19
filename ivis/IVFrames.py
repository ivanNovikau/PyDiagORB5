import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr
import ivis.IVButtons as ivb
import ivis.IVMenus as ivm

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

    # style:
    bg = None  # background color

    def __init__(self, mw, **kwargs):
        super(BFrame, self).__init__(**kwargs)
        self.mw = mw


# *** Upper Frame (with File, Options etc.) ***
class UpperFrame(BFrame):

    def __init__(self, mw, **kwargs):
        super(UpperFrame, self).__init__(mw, **kwargs)

        # styles:
        self.bg = mix.to_rgb((180, 180, 180))
        self.config(
            bg=self.bg
        )


# *** Figure Frame (where figures are plotted) ***
class FigureFrame(BFrame):
    # style
    bg_canvas = mix.to_rgb((140, 140, 140))

    # frames
    fUpper = None
    fBottom = None

    # popup that appear then right button is pressed
    popup_menu_canvas = None

    # canvas where figure is plotted
    fig_canvas = None

    def __init__(self, mw, **kwargs):
        super(FigureFrame, self).__init__(mw, **kwargs)

        # styles:
        self.bg = mix.to_rgb((140, 140, 140))
        self.config(
            bg=self.bg
        )

        # upper figure frame
        self.fUpper = UpperFigureFrame(
            mw=self.mw,
            figure_frame=self,
            master=self,
            height=30
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
        self.mw.fig.patch.set_facecolor(self.bg_canvas)

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


# *** Upper Figure Frame (with figure tabs) ***
class UpperFigureFrame(BFrame):
    fFigure = None  # figure frame

    def __init__(self, mw, figure_frame, **kwargs):
        super(UpperFigureFrame, self).__init__(mw, **kwargs)
        self.fFigure = figure_frame

        # styles:
        self.bg = mix.to_rgb((160, 160, 160))
        self.config(
            bg=self.bg
        )


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
            bg=mix.to_rgb(GLO.IVIS_label_color)
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

        # styles:
        self.bg = mix.to_rgb((120, 120, 120))
        self.config(
            bg=self.bg
        )

        # Add frames:
        self.fFigProp = FigPropFrame(self.mw, master=self, bg=self.bg)
        self.fAxProp = AxPropFrame(self.mw, master=self, bg=self.bg)
        self.fProc = ProcessingFrame(self.mw, master=self, bg=self.bg)

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


# *** Frame with Axes properties  ***
class AxPropFrame(BFrame):
    def __init__(self, mw, **kwargs):
        super(AxPropFrame, self).__init__(mw, **kwargs)


# *** Frame to work with post-processing  ***
class ProcessingFrame(BFrame):
    def __init__(self, mw, **kwargs):
        super(ProcessingFrame, self).__init__(mw, **kwargs)






