import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr
import ivis.IVButtons as ivb
import ivis.IVMenus as ivm
import ivis.IVPageText as ivp_text
import ivis.IVPageCurves as ivp_curves
import ivis.IVPageData as ivp_data
import ivis.IVPageFigures as ivp_figures
import ivis.IVPageAxes as ivp_axes

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
    mix.reload_module(ivp_text)
    mix.reload_module(ivp_curves)
    mix.reload_module(ivp_data)
    mix.reload_module(ivp_figures)
    mix.reload_module(ivp_axes)


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


# *** Left Frame (with Figure, Axes properties and processing) ***
class LeftFrame(BFrame):
    names_sections = ['fig', 'ax', 'proc']
    sections = {}

    def __init__(self, mw, **kwargs):
        super(LeftFrame, self).__init__(mw, **kwargs)

        # --- create sections' dictionaries ---
        self.sections['fig'] = {}
        self.sections['ax'] = {}
        self.sections['proc'] = {}

        # --- create frames ---
        cnt = mix.Counter()
        self.columnconfigure(0, weight=1)
        for oname in self.names_sections:
            sec = self.sections[oname]
            frame = self.sections[oname]['frame'] = tk.Frame(master=self)
            sec['pages'] = {}

            # style
            frame.config(
                highlightbackground=GLO.IVIS_border_color,
                highlightthickness=2,
            )

            # position and relative size
            frame.grid(
                row=cnt.next(), column=0, sticky=tk.N + tk.S + tk.E + tk.W
            )
            self.rowconfigure(cnt.counter, weight=1)

            # tab controller
            sec['tc'] = ivb.TabController(frame)

        # create elements in the sections:
        self.create_section_fig()
        self.create_section_ax()

    def create_section_fig(self):
        name_section = 'fig'
        sec = self.sections[name_section]
        tc, pages = sec['tc'], sec['pages']

        # Create Data page
        pages['data'] = ivp_data.PageData(
            frame=tc.create_page('data'),
            mw=self.mw,
        )

        # Create Figures page
        pages['figures'] = ivp_figures.PageFigures(
            frame=tc.create_page('figures'),
            mw=self.mw,
        )

        # Create Axes page
        pages['axes'] = ivp_axes.PageAxes(
            frame=tc.create_page('axes'),
            mw=self.mw,
        )

        # activate Data tab:
        tc.call_page('data')

    def create_section_ax(self):
        name_section = 'ax'
        sec = self.sections[name_section]
        tc, pages = sec['tc'], sec['pages']

        # Create Text page
        pages['text'] = ivp_text.PageText(
            frame=tc.create_page('text'),
            mw=self.mw,
        )

        # Create Curves page
        pages['curves'] = ivp_curves.PageCurves(
            frame=tc.create_page('curves'),
            mw=self.mw,
        )

        # activate Text tab:
        tc.call_page('text')




