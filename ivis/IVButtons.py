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

        rowspan, columnspan = 1, 1
        if len(pos_row_col) > 2:
            rowspan = pos_row_col[2]
            if len(pos_row_col) > 3:
                columnspan = pos_row_col[3]

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
            row=pos_row_col[0], column=pos_row_col[1],
            rowspan=rowspan,
            sticky=tk.N + tk.S + tk.E + tk.W,
        )
        self.entry.grid(
            row=pos_row_col[0], column=pos_row_col[1] + 1,
            rowspan=rowspan, columnspan=columnspan,
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


# *** Table ***
class Table:
    title = ''
    master = None
    frame = None
    bCreateHide = None
    bAddRow, bAddColumn = None, None
    mRow, mColumn = None, None
    rowNames, columnNames = None, None

    row_cnt    = None  # to count cells in a row
    column_cnt = None  # to count cells in a column

    cells = []  # variables[id_row, id_column]

    # find correspondence between a cell position and a variable
    map_row_ids, map_column_ids = {}, {}

    def __init__(self, title, master, pos_size, rows_list, columns_list,
                 default_rows_columns=(0, 0, 1, 1)):
        self.title = title
        self.master = master
        row, column, rowspan, columnspan = pos_size

        # --- start button ---
        self.bCreateHide = BButton(
            master=self.master,
            text='Create ' + title,
            command=self.create_table
        )
        self.bCreateHide.grid(
            row=row,        rowspan=1,
            column=column,  columnspan=columnspan,
            sticky=tk.N + tk.S + tk.E + tk.W
        )

        # --- Frame ---
        self.frame = tk.Frame(
            self.master,
            bg=GLO.IVIS_selected_element_inf_frame
        )
        self.frame.grid(
            row=row+1,      rowspan=rowspan,
            column=column,  columnspan=columnspan,
            sticky=tk.N + tk.S + tk.E + tk.W
        )

        # --- Buttons to add rows and columns ---
        self.mRow = tk.Menu(master=self.frame, tearoff=0)
        self.rowNames = list(rows_list)
        for id_row, one_row in enumerate(rows_list):
            self.mRow.add_command(
                label=one_row,
                command=lambda id=id_row: self.create_new_row(id)
            )

        self.mColumn = tk.Menu(master=self.frame, tearoff=0)
        self.columnNames = list(columns_list)
        for id_column, one_column in enumerate(columns_list):
            self.mColumn.add_command(
                label=one_column,
                command= lambda id=id_column: self.create_new_column(id)
            )

        self.bAddRow = BButton(
            master=self.frame,
            text='Add row',
            command=lambda flag_row=True, menu=self.mRow:
                self.call_menu(flag_row, menu)
        )
        self.bAddColumn = BButton(
            master=self.frame,
            text='Add column',
            command=lambda flag_row=False, menu=self.mColumn:
                self.call_menu(flag_row, menu)
        )

        # create default rows and columns
        self.create_default_table(default_rows_columns)

    def create_table(self):
        self.bCreateHide['state'] = tk.DISABLED
        self.bCreateHide['text'] = self.title

        self.bAddRow.grid(row=0, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.bAddColumn.grid(row=0, column=1, sticky=tk.N + tk.S + tk.E + tk.W)

        self.frame.rowconfigure(0, weight=1)
        self.frame.columnconfigure(0, weight=1)
        self.frame.columnconfigure(1, weight=1)

        self.row_cnt = mix.Counter()
        self.column_cnt = mix.Counter()

    def create_default_table(self, default_rows_columns):
        start_row, start_column, n_rows, n_columns = default_rows_columns

        self.create_table()

        for id_row in range(start_row, start_row + n_rows):
            self.create_new_row(id_row)

        for id_column in range(start_column, start_column + n_columns):
            self.create_new_column(id_column)

    def call_menu(self, flag_row, menu):
        try:
            button = self.bAddRow if flag_row else self.bAddColumn
            x, y = button.winfo_rootx(), button.winfo_rooty()
            menu.tk_popup(x+40, y+30, 0)
        finally:
            menu.grab_release()

    def create_new_row(self, id):
        # do not double a row:
        if id in self.map_row_ids.keys():
            return

        # increment the row counter
        self.row_cnt.next()

        # position of the row label and row cells
        row_pos = self.row_cnt.counter + 2
        start_column_pos = 0

        # label of the row
        row_label = BLabel(master=self.frame, text=self.rowNames[id])
        row_label.grid(
            row=row_pos, column=start_column_pos,
            sticky=tk.N + tk.S + tk.E + tk.W
        )

        # create cells in the row for every column
        cells_in_row = []
        for id_column in range(self.column_cnt.n_elements):
            cells_in_row.append(tk.StringVar())
            tk.Entry(
                master=self.frame, textvariable=cells_in_row[-1]
            ).grid(
                row=row_pos, column=start_column_pos + id_column+1,
                sticky=tk.N + tk.S + tk.E + tk.W
            )
        self.cells.append(cells_in_row)

        # relate cell position with actual variable
        self.map_row_ids[id] = self.row_cnt.counter

    def create_new_column(self, id):
        # do not double a column
        if id in self.map_column_ids.keys():
            return

        # increment the coumn counter
        self.column_cnt.next()

        # position of the column label and column cells
        column_pos = self.column_cnt.counter + 1
        start_row_pos = 1

        column_label = BLabel(master=self.frame, text=self.columnNames[id])
        column_label.grid(
            row=start_row_pos, column=column_pos,
            sticky=tk.N + tk.S + tk.E + tk.W
        )

        # create cells in the column for every row
        for id_row in range(self.row_cnt.n_elements):
            row_cell = self.cells[id_row]
            row_cell.append(tk.StringVar())
            tk.Entry(
                master=self.frame, textvariable=row_cell[-1]
            ).grid(
                row=start_row_pos + id_row +1, column=column_pos,
                sticky=tk.N + tk.S + tk.E + tk.W
            )

        for id_row in range(self.row_cnt.n_elements):
            id_row_pos = id_row + start_row_pos
            self.frame.rowconfigure(id_row_pos, weight=1)
            self.frame.columnconfigure(column_pos, weight=1)

        # relate cell position with actual variable
        self.map_column_ids[id] = self.column_cnt.counter

    def set_cell(self, id_row_var, id_column_var, v):
        # id_row_var corresponds to a position of a necessary variable (name) in rowNames
        # id_column_var corresponds to a position of a necessary variable (name)
        #   in columnNames

        if id_row_var not in self.map_row_ids.keys():
            return
        if id_column_var not in self.map_column_ids.keys():
            return

        id_row_cell = self.map_row_ids[id_row_var]
        id_column_cell =  self.map_column_ids[id_column_var]

        self.cells[id_row_cell][id_column_cell].set(v)
