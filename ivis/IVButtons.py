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


# *** Basic Button ***
class BButton(tk.Button):

    mw = None  # Main Window

    # button style:
    bg = mix.to_rgb((160, 160, 160))
    padx = 10
    relief = tk.GROOVE
    activebackground = mix.to_rgb((120, 120, 120))

    def __init__(self, mw, **kwargs):
        super(BButton, self).__init__(
            bg=self.bg,
            padx=self.padx,
            relief=self.relief,
            activebackground=self.activebackground,
            **kwargs
        )
        self.mw = mw


# *** Quit button ***
class QuitButton(BButton):

    def __init__(self, mw, **kwargs):
        super(QuitButton, self).__init__(mw, **kwargs)
        self.grid(row=0, column=0)



