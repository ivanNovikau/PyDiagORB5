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
    bg_active = GLO.IVIS_color_active_button

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


# *** Labelled entry ***
class LabelledEntry:
    master = None
    label = None
    entry = None
    var = None

    def __init__(self, master, text, pos_row_col, init_entry=''):
        self.master = master
        self.label = BLabel(master=self.master, text=text)

        self.var = tk.StringVar(value=init_entry)
        self.entry = tk.Entry(master=self.master, textvariable=self.var)

        # set position
        self.label.grid(
            row=pos_row_col[0], column=pos_row_col[-1],
            sticky=tk.N + tk.S + tk.E + tk.W,
        )
        self.entry.grid(
            row=pos_row_col[0], column=pos_row_col[-1]+1,
            sticky=tk.N + tk.S + tk.E + tk.W,
        )


# *** Labelled OptionMenu ***
class LabelledOptionMenu:
    master = None
    label = None
    entry = None
    var = None
    ids_elements = {}
    call_selected = None

    def __init__(
            self, master, text, pos_row_col, options, call_selected=None, init_var=None
    ):
        # Label
        self.master = master
        self.label = BLabel(master=self.master, text=text)
        self.call_selected = call_selected

        # OptionMenu
        res_options = self.form_options(options)

        self.var = tk.StringVar()
        if init_var is not None:
            self.var.set(init_var)

        self.entry = tk.OptionMenu(
            self.master,
            self.var,
            *res_options
        )
        self.var.trace(
            'w',
            self.option_selected
        )

        # set position
        self.label.grid(
            row=pos_row_col[0], column=pos_row_col[-1],
            sticky=tk.N + tk.S + tk.E + tk.W,
        )
        self.entry.grid(
            row=pos_row_col[0], column=pos_row_col[-1] + 1,
            sticky=tk.N + tk.S + tk.E + tk.W,
        )

    def form_options(self, options):
        res_options = ["---"]
        if len(options) > 0:
            res_options = []
            for id_opt, one_opt in enumerate(options):
                if self.call_selected is not None:
                    one_opt = '{:d}: '.format(id_opt) + one_opt
                    self.ids_elements[one_opt] = id_opt
                res_options.append(one_opt)
        return res_options

    def option_selected(self, *args):
        if self.call_selected is not None:
            opt_selected = self.var.get()
            if opt_selected == "---":
                id_selected = None
            else:
                id_selected = self.ids_elements[opt_selected]
            self.call_selected(opt_selected, id_selected)

    def update_options(self, options):
        menu = self.entry["menu"]
        menu.delete(0, "end")
        res_options = self.form_options(options)

        for one_text in res_options:
            menu.add_command(
                label=one_text,
                command=lambda value=one_text: self.var.set(value)
            )
        self.var.set(res_options[-1])






