import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr
import ivis.IVButtons as ivb
import ivis.IVMenus as ivm

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

    # to view curves XYZ data
    flag_curves_xyz_data = False

    # to zoom data
    flag_zoom = False
    xy_zoom = None

    n_temp_lines = 0

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
        self.fig_canvas.mpl_connect('motion_notify_event', self.on_motion)

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
        if event.button == 1:
            if self.flag_zoom:

                # start drawing a rectangle:
                if self.xy_zoom is None:
                    self.xy_zoom = {
                        'x1': event.xdata,
                        'y1': event.ydata
                    }
                # end drawing a rectangle
                else:
                    self.xy_zoom.update({
                        'x2': event.xdata,
                        'y2': event.ydata
                    })

                    # change axis limits
                    x1, x2, y1, y2 = \
                        self.xy_zoom['x1'], self.xy_zoom['x2'], \
                        self.xy_zoom['y1'], self.xy_zoom['y2']

                    x_lim = (x1, x2) if x1 < x2 else (x2, x1)
                    y_lim = (y1, y2) if y1 < y2 else (y2, y1)

                    ax = self.mw.get_ax()
                    ax.set_xlim(x_lim)
                    ax.set_ylim(y_lim)

                    # reset to None XY coordinates of the zoom box
                    self.xy_zoom = None

                    # update the axis
                    self.mw.draw()

        if event.button == 3:
            self.popup_menu_canvas.call(event)

        button_press_handler(event, self.fig_canvas)

    def on_motion(self, event):
        fB = self.mw.fFigure.fBottom
        fB.elements['xdata'].set(event.xdata)
        fB.elements['ydata'].set(event.ydata)
        fB.elements['xgui'].set(event.guiEvent.x_root)
        fB.elements['ygui'].set(event.guiEvent.y_root)

        if self.flag_curves_xyz_data:
            self.view_xyz(event)
        if self.flag_zoom:
            self.zoom_data(event)

    def view_xyz(self, event):
        if event.xdata is not None and event.ydata is not None:
            self.n_temp_lines = 0
            ax = self.mw.get_ax()

            x_lim = ax.get_xlim()
            y_lim = ax.get_ylim()

            # draw horizontal and vertical lines
            self.n_temp_lines += 1
            one_line = mlines.Line2D(
                [event.xdata, event.xdata], y_lim,
                color='grey',
                linestyle=':',
                linewidth=2
            )
            ax.add_line(one_line)

            # draw markers on curves
            for id_curve, one_curve in enumerate(self.mw.curves.list_curves):
                self.n_temp_lines += 1
                id_x, x_curve, _ = mix.get_ids(one_curve.xs, event.xdata)
                y_curve = one_curve.ys[id_x]

                color_plt = one_curve.ff['color'] \
                    if one_curve.ff['color'] is not None \
                    else GLO.new_color(id_curve)

                ax.plot(
                    x_curve, y_curve,
                    'o',
                    color=color_plt,
                    markersize=10,
                )

                self.n_temp_lines += 1
                one_line = mlines.Line2D(
                    [x_lim[0], x_curve], [y_curve, y_curve],
                    color=color_plt,
                    linestyle=':',
                    linewidth=2
                )
                ax.add_line(one_line)

            ax.set_xlim(x_lim)
            ax.set_ylim(y_lim)

            # update the axis
            self.mw.draw()

            # remove lines and markers
            for one_line in ax.lines[-self.n_temp_lines:-1]:
                one_line.remove()
            ax.lines[-1].remove()
        else:
            if self.n_temp_lines > 0:
                self.n_temp_lines = 0
                self.mw.draw()

    def zoom_data(self, event):
        if self.xy_zoom is None:
            return

        if event.xdata is not None and event.ydata is not None:
            ax = self.mw.get_ax()

            left_bottom_corner, w, h = mix.get_rectangular(
                self.xy_zoom['x1'], event.xdata,
                self.xy_zoom['y1'], event.ydata,
            )

            rect = patches.Rectangle(
                left_bottom_corner, w, h,
                ls=':',
                linewidth=2,
                edgecolor='grey',
                facecolor='none',
            )
            ax.add_patch(rect)

            # update the axis
            ax.xaxis.set_major_locator(AutoLocator())
            ax.yaxis.set_major_locator(AutoLocator())
            self.mw.draw()

            # delete the rectangular
            ax.patches[-1].remove()
        else:
            pass


# *** Bottom Figure Frame (with some properties) ***
class BottomFigureFrame(BFrame):
    fFigure = None  # figure frame

    # label with axes coordinates
    elements = {}
    lbXYdata = None

    def __init__(self, mw, figure_frame, **kwargs):
        super(BottomFigureFrame, self).__init__(mw, **kwargs)
        self.fFigure = figure_frame

        # Entries:
        self.elements['xdata'] = ivb.LabelledEntry(self, "X-data", [0, 0]).var
        self.elements['ydata'] = ivb.LabelledEntry(self, "Y-data", [0, 2]).var
        self.elements['xgui'] = ivb.LabelledEntry(self, "X-GUI", [0, 4]).var
        self.elements['ygui'] = ivb.LabelledEntry(self, "Y-GUI", [0, 6]).var

        # Buttons:
        self.elements['xyz'] = ivb.BButton(
            master=self,
            text='XYZ',
            command=self.view_xyz
        )
        self.elements['xyz'].grid(row=0, column=8)

        self.elements['zoom'] = ivb.BButton(
            master=self,
            text='zoom',
            command=self.zoom_data
        )
        self.elements['zoom'].grid(row=0, column=9)

    def view_xyz(self):
        self.mw.fFigure.flag_curves_xyz_data = \
            not self.mw.fFigure.flag_curves_xyz_data

        if self.mw.fFigure.flag_curves_xyz_data:
            self.elements['xyz'].configure(
                bg=GLO.IVIS_color_active_button
            )
        else:
            self.elements['xyz'].configure(
                bg=GLO.IVIS_color_button
            )

    def zoom_data(self):
        self.mw.fFigure.flag_zoom = not self.mw.fFigure.flag_zoom

        if self.mw.fFigure.flag_zoom:
            self.elements['zoom'].configure(
                bg=GLO.IVIS_color_active_button
            )
        else:
            self.elements['zoom'].configure(
                bg=GLO.IVIS_color_button
            )


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
    bUpdatePlot = None

    def __init__(self, mw, **kwargs):
        super(FigPropFrame, self).__init__(mw, **kwargs)

        # style
        self.config(
            highlightbackground=GLO.IVIS_border_color,
            highlightthickness=2,
        )

        # Buttons
        self.bUpdatePlot = ivb.BButton(
            master=self,
            text="Default plot",
            command=self.default_plot
        )

        # arrange elements
        self.bUpdatePlot.grid(row=0, column=0)

    def default_plot(self):
        # --- create a figure and plot data ---
        mw = self.mw
        ax = mw.get_ax()

        mw.curves = crv.copy_curves(mw.curves_default, ax)

        mw.fig, _, _ = cpr.plot(
            mw.curves_default, mw.fig, ax,
            FIG_W=mw.curves.ff['figure_width'] / 2,
            FIG_H=mw.curves.ff['figure_height'] / 2,
        )
        self.mw.draw()

        # update elemens in left panel:
        self.mw.fLeft.fAxProp.update_elements()


# *** Frame with Axes properties  ***
class AxPropFrame(BFrame):
    tabController = None
    n_columns = 5
    omAText = None

    # frames with information about different elements
    fElementInf = None  # root element frame
    feiText = None
    main_els = {}
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

    def rebuild_plot(self):
        ivis_add_xticks = list(mix.array_from_str(self.main_els['xaticks'].get()))
        ivis_add_yticks = list(mix.array_from_str(self.main_els['yaticks'].get()))
        self.mw.curves.ff.update({
            'title': self.main_els['title'].get(),
            'xlabel': self.main_els['xlabel'].get(),
            'ylabel': self.main_els['ylabel'].get(),
            'xticks': self.mw.curves_default.ff['xticks'] + ivis_add_xticks,
            'yticks': self.mw.curves_default.ff['yticks'] + ivis_add_yticks,
            'ivis_add_xticks': ivis_add_xticks,
            'ivis_add_yticks': ivis_add_yticks,
        })

        # update the plot format
        ax = self.mw.get_ax()
        ax.texts = []
        cpr.format_plot(
            self.mw.fig, ax, self.mw.curves, self.mw.flag_2d
        )
        self.mw.draw()

    def update_elements(self):
        ax = self.mw.get_ax()
        self.mw.curves = crv.copy_curves(self.mw.curves_default, ax)
        curves = self.mw.curves

        self.main_els['title'].set(curves.ff['title'] if curves.ff['title'] is not None else "")
        self.main_els['xlabel'].set(curves.ff['xlabel'] if curves.ff['xlabel'] is not None else "")
        self.main_els['ylabel'].set(curves.ff['ylabel'] if curves.ff['ylabel'] is not None else "")
        self.main_els['xaticks'].set("")
        self.main_els['yaticks'].set("")
        self.omAText.update_options(
            mix.get_atext_from_curves(curves)
        )

    def default_plot(self):
        self.update_elements()

        ax = self.mw.get_ax()
        ax.texts = []
        cpr.format_plot(
            self.mw.fig, ax, self.mw.curves, self.mw.flag_2d
        )
        self.mw.draw()

    def create_frame_text(self):
        # --- FUNCTIONS ------------------------------------------------------------------
        def text_selected(opt_selected, id_selected):
            self.id_current_text = id_selected
            if self.id_current_text is not None:
                el_text = self.mw.curves.list_text[id_selected]

                self.feiText_elements['text'].set(
                    mix.delete_bold_keyword(el_text.line)
                )
                self.feiText_elements['x'].set(el_text.x)
                self.feiText_elements['y'].set(el_text.y)
                self.feiText_elements['color'].var.set(el_text.color)
                self.feiText_elements['flag_invisible'].set(el_text.flag_invisible)

                self.feiText.tkraise()
            else:
                self.feiText_elements['text'].set('')
                self.feiText_elements['x'].set('')
                self.feiText_elements['y'].set('')
                self.feiText_elements['color'].var.set('black')
                self.feiText_elements['flag_invisible'].set(0)

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
        self.main_els['title'] = ivb.LabelledEntry(cf, "Title: ", [0, 0], curves.ff['title']).var
        self.main_els['xlabel'] = ivb.LabelledEntry(cf, "X label: ", [1, 0], curves.ff['xlabel']).var
        self.main_els['ylabel'] = ivb.LabelledEntry(cf, "Y label: ", [2, 0], curves.ff['ylabel']).var
        self.main_els['xaticks'] = ivb.LabelledEntry(cf, "X additional ticks: ", [3, 0], "").var
        self.main_els['yaticks'] = ivb.LabelledEntry(cf, "Y additional ticks: ", [4, 0], "").var

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
        self.fElementInf = tk.Frame(  # root frame
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
            command=self.rebuild_plot
        ).grid(row=n_rows+1, column=0)
        ivb.BButton(
            master=cf,
            text='Default',
            command=self.default_plot
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






