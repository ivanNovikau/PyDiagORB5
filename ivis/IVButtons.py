import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr

import tkinter as tk


# *** Buttons, Labels, Tab Controller ***


def reload():
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)
    mix.reload_module(GLO)
    mix.reload_module(cpr)


# *** Basic Button ***
class BButton(tk.Button):
    bg = GLO.IVIS_color_button
    padx = 10
    relief = tk.GROOVE
    activebackground = mix.to_rgb((120, 120, 120))

    def __init__(self, **kwargs):
        super(BButton, self).__init__(
            bg=self.bg,
            padx=self.padx,
            relief=self.relief,
            activebackground=self.activebackground,
            **kwargs
        )


# *** Basic Label ***
class BLabel(tk.Label):
    bg = GLO.IVIS_label_color

    def __init__(self, **kwargs):
        super(BLabel, self).__init__(
            bg=self.bg,
            **kwargs
        )


# *** Tab Controller ***
class TabController:
    master = None
    fTabs = None
    fMain = None

    pages = {}
    tab_buttons = {}

    n_pages = 0

    # name of active page
    name_active = None

    # color of active button
    bg_active = mix.to_rgb((203, 203, 255))

    def __init__(self, master):
        # parent frame
        self.master = master

        # tabs and main frames
        self.fTabs = tk.Frame(
            master=self.master, bg=GLO.IVIS_color_tabs_frame,
            height=GLO.IVIS_height_tab_frame
        )
        self.fMain = tk.Frame(
            master=self.master, bg=GLO.IVIS_frame_color
        )

        # arrange the frames
        self.fTabs.grid(row=0, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.fMain.grid(row=1, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.master.columnconfigure(0, weight=1)
        self.master.rowconfigure(0, weight=0)
        self.master.rowconfigure(1, weight=1)

    def create_page(self, name, bg=GLO.IVIS_frame_color):
        # create a new page:
        new_page = tk.Frame(master=self.fMain, bg=bg)
        self.pages[name] = new_page
        self.n_pages += 1

        new_page.grid(row=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.fMain.columnconfigure(0, weight=1)
        self.fMain.rowconfigure(0, weight=1)

        # create a corresponding button:
        new_button = BButton(
            master=self.fTabs,
            text=name,
            command=lambda: self.call_page(name),
        )
        self.tab_buttons[name] = new_button

        new_button.grid(column=self.n_pages - 1, row=0)

        # set the first create page as the selected one:
        if self.n_pages == 1:
            self.call_page(name)

        return new_page

    def call_page(self, page_name):
        # check is it a new active page:
        if page_name == self.name_active:
            return

        # reset color of previous active button
        if self.name_active is not None:
            self.tab_buttons[self.name_active].configure(
                bg=GLO.IVIS_color_button
            )

        # set new active page
        self.name_active = page_name
        self.pages[self.name_active].tkraise()

        self.tab_buttons[self.name_active].configure(
            bg=self.bg_active
        )







