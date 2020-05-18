import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO

import os
from matplotlib import rcParams
import re
import pandas as pd
import subprocess
import shlex

import matplotlib.pyplot as mpl
import numpy as np
import tkinter as tk
import tkinter.filedialog
import tkinter.messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler, button_press_handler


def reload():
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)
    mix.reload_module(GLO)


class Ivis:
    root = None  # main window
    canvas = None  # canvas with a plot
    curves = None

    fig = None  # figure

    ext_data = '.dat'
    ext_latex = '.tex'
    ext_png = '.png'
    ext_eps = '.eps'

    def __init__(self, fig, curves, **kwargs):
        # curves to plot:
        self.curves = curves

        # create the main window of the application (top widget)
        self.root = tk.Tk()
        self.root.wm_title("ivis")

        # Frame with a plot and basic toolbar
        figFrame = tk.Frame(master=self.root)

        # create plot
        self.fig = fig
        self.canvas = FigureCanvasTkAgg(fig, master=figFrame)  # A tk.DrawingArea.
        self.canvas.draw()

        # add the standard matpotlib toolbar
        toolbarFrame = tk.Frame(master=figFrame)
        toolbar = NavigationToolbar2Tk(self.canvas, toolbarFrame)
        toolbar.update()

        self.canvas.mpl_connect(
            "key_press_event",
            lambda event: self.on_key_press(event, self.canvas, toolbar)
        )
        self.canvas.mpl_connect(
            "button_press_event",
            lambda event: self.on_button_press(event, self.canvas, toolbar)
        )

        # create a quit button
        bQuit = tk.Button(master=self.root, text="Quit", command=self.root.destroy)
        bPgf = tk.Button(master=self.root, text="Save pgfplot",
                         command=self.on_press_pgfplot_save)

        # create the window grid
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky=tk.N+tk.S+tk.E+tk.W)
        toolbarFrame.grid(row=1, column=0)
        figFrame.columnconfigure(0, weight=1)
        figFrame.rowconfigure(0, weight=1)

        bQuit.grid(row=0, column=0)
        bPgf.grid(row=1, column=0)
        figFrame.grid(row=0, column=1, rowspan=2, sticky=tk.N+tk.S+tk.E+tk.W)
        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=1)

        tk.mainloop()

    def on_key_press(self, event, canvas, toolbar):
        print("you pressed {}".format(event.key))
        key_press_handler(event, canvas, toolbar)

    def on_button_press(self, event, canvas, toolbar):
        print("you pressed {}".format(event.xdata))
        button_press_handler(event, canvas, toolbar)

    def on_press_pgfplot_save(self):
        fname = self.get_pgfplot_filename()
        if fname in ["", ()]:
            return

        # save the file
        try:
            self.save_pgfplot(fname)
        except Exception as e:
            tk.messagebox.showerror("Error saving file", str(e))

    def get_pgfplot_filename(self):
        defaultextension = ''
        initialdir = os.path.expanduser(rcParams['savefig.directory'])
        initialfile = self.canvas.get_default_filename()

        fname = tk.filedialog.asksaveasfilename(
            master=self.root,
            title='Save the figure',
            defaultextension=defaultextension,
            initialdir=initialdir,
            initialfile=initialfile,
        )

        # Save dir for next time, unless empty str (i.e., use cwd).
        if initialdir != "":
            rcParams['savefig.directory'] = (
                os.path.dirname(str(fname)))

        return fname

    def save_pgfplot(self, fname):
        # save a bunch of files necessary to recreate the pgfplot

        # save 1d plot
        if self.curves.list_curves[0].zs is None:
            self.save_pgfplot_1d(fname)
        # save 2d plot
        else:
            self.save_pgfplot_2d(fname)

    def save_pgfplot_1d(self, fname):
        print('Save the 1d plot to files with a name: {}'.format(fname))

        # save data to a .dat file:
        file_data_name_curves = []
        for id_curve, one_curve in enumerate(self.curves.list_curves):
            file_name = fname + "{:d}".format(id_curve) + self.ext_data
            file_data_name_curves.append(file_name)
            output_df = pd.DataFrame({
                'X': one_curve.xs,
                'Y': one_curve.ys
            })
            output_df.to_csv(file_name, sep=" ", index=False)

        # read template to create a .tex file
        ff_template = open("ivis/template_plot_1d.txt", 'r')
        template_text = ff_template.read()
        ff_template.close()

        # --- save a corresponding .tex file ---
        file_name = fname + self.ext_latex
        ff_tex = open(file_name, 'w')

        result_text = template_text

        # set title
        result_text = self.template_set(
            result_text, 1, self.curves.ff['title']
        )

        # set XY labels
        result_text = self.template_set(
            result_text, 2, self.curves.ff['xlabel']
        )
        result_text = self.template_set(
            result_text, 3, self.curves.ff['ylabel']
        )

        # set file name where data are read from
        result_data_line = ''
        for id_curve, file_data_name in enumerate(file_data_name_curves):
            line_file_name = '\\addplot table [y=Y, x=X]{' \
                + '{}'.format(file_data_name) \
                + '};\n'
            line_legend = '\\addlegendentry{'\
                + '{}'.format(self.curves.list_curves[id_curve].ff['legend'])\
                + '};\n'
            result_data_line += line_file_name + line_legend
        result_data_line = result_data_line[:-1]

        result_text = self.template_set(
            result_text, 4, result_data_line
        )

        # set additional XY ticks:
        line_ticks = ''
        if not np.isnan(self.curves.ff['xticks']).any():
            line_ticks = '{}'.format(line_ticks)
            line_ticks = line_ticks[1:-1]
        result_text = self.template_set(
            result_text, 5, line_ticks
        )

        line_ticks = ''
        if not np.isnan(self.curves.ff['yticks']).any():
            line_ticks = '{}'.format(line_ticks)
            line_ticks = line_ticks[1:-1]
        result_text = self.template_set(
            result_text, 6, line_ticks
        )

        # save the resulting text into the .tex file
        ff_tex.write(result_text)
        ff_tex.close()

    def save_pgfplot_2d(self, fname):
        print('Save the 2d plot to files with a name: {}'.format(fname))

        # --- extract necessary data ---
        ref_curve = self.curves.list_curves[0]
        xmin, xmax = np.min(ref_curve.xs), np.max(ref_curve.xs)
        ymin, ymax = np.min(ref_curve.ys), np.max(ref_curve.ys)
        zmin, zmax = np.min(ref_curve.zs), np.max(ref_curve.zs)

        # --- create external figure .eps (using Image Magic) ---
        name_plot_png = fname + self.ext_png
        self.fig.savefig(name_plot_png)

        name_plot_eps = fname + self.ext_eps
        command_line = 'convert -trim ' + name_plot_png + ' ' + name_plot_eps
        os.system(command_line)

        # --- read template to create a .tex file ---
        ff_template = open("ivis/template_plot_2d.txt", 'r')
        template_text = ff_template.read()
        ff_template.close()

        # --- save the corresponding .tex file ---
        file_name = fname + self.ext_latex
        ff_tex = open(file_name, 'w')

        result_text = template_text

        # set the plot file name to read
        result_text = self.template_set(
            result_text, 1, name_plot_eps
        )

        # set title
        result_text = self.template_set(
            result_text, 2, self.curves.ff['title']
        )

        # set XY labels
        result_text = self.template_set(
            result_text, 3, self.curves.ff['xlabel']
        )
        result_text = self.template_set(
            result_text, 4, self.curves.ff['ylabel']
        )

        # set min and max values of axes:
        result_text = self.template_set(result_text, 5, xmin)
        result_text = self.template_set(result_text, 6, xmax)
        result_text = self.template_set(result_text, 7, ymin)
        result_text = self.template_set(result_text, 8, ymax)
        result_text = self.template_set(result_text, 9, zmin)
        result_text = self.template_set(result_text, 10, zmax)

        # set the colormap:
        resulting_colormap = ref_curve.ff['colormap'] \
            if ref_curve.ff['colormap'] is not None \
            else GLO.DEF_COLORMAP
        result_text = self.template_set(
            result_text, 11, resulting_colormap
        )

        # save the resulting text into the .tex file
        ff_tex.write(result_text)
        ff_tex.close()

    def template_set(self, template_text, id_to_change, value):
        # create identifiers
        id_left = "IVIS{:d}N".format(id_to_change)
        id_right =  "N{:d}IVIS".format(id_to_change)

        # find a line between the identifiers
        line_to_change = re.search(
            id_left + "(.*)" + id_right,
            template_text
        ).group(1)

        # change the line
        changed_line = line_to_change.format(value)

        # update the text created from the template
        result_text = template_text.replace(
            id_left + line_to_change + id_right, changed_line
        )
        return result_text





