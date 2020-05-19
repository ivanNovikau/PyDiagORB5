import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr

import tkinter as tk


def reload():
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)
    mix.reload_module(GLO)
    mix.reload_module(cpr)


# *** Basic Menu ***
class BMenu(tk.Menu):
    # Main window (where fig and curves can be found)
    mw = None

    def __init__(self, mw, **kwargs):
        super(BMenu, self).__init__(**kwargs)
        self.mw = mw


# *** Popup Canvas Menu ***
class PopupCanvasMenu(BMenu):

    # coordinates
    xdata, ydata = None, None  # read axes coordinates from figure
    xroot, yroot = None, None  # gui coordinates

    def __init__(self, mw, **kwargs):
        super(PopupCanvasMenu, self).__init__(mw, **kwargs)
        self.add_command(
            label='Add text',
            command=self.add_text_here
        )
        self.add_command(
            label='Delete previous text',
            command=self.delete_previous_text
        )

    def call(self, event):
        try:
            self.xdata = event.xdata
            self.ydata = event.ydata
            self.xroot = event.guiEvent.x_root
            self.yroot = event.guiEvent.y_root
            self.tk_popup(self.xroot, self.yroot, 0)
        finally:
            self.grab_release()

    def add_text_here(self):
        axes = self.mw.fig.axes
        ax = axes[0]

        ax.text(
            self.xdata,
            self.ydata,
            'Text here',
            fontsize=GLO.FONT_SIZE,
        )
        self.mw.fig.canvas.draw()

    def delete_previous_text(self):
        txt = self.mw.fig.axes[0].texts
        txt[-1].set_visible(False)
        self.mw.fig.canvas.draw()
