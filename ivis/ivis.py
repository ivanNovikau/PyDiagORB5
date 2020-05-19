import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr
import ivis.IVButtons as ivb
import ivis.IVFrames as ivf

import os
from matplotlib import rcParams
import re
import pandas as pd

import numpy as np
import tkinter as tk
import tkinter.filedialog
import tkinter.messagebox

from matplotlib.backend_bases import key_press_handler, button_press_handler


def reload():
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)
    mix.reload_module(GLO)
    mix.reload_module(cpr)
    mix.reload_module(ivb)
    mix.reload_module(ivf)


def plot_data(curves, fig=None, ax=None):
    # --- modify font size if ivis is used ---
    if curves.ff['flag_ivis']:
        curves.ff['fontS'] /= GLO.IVIS_COEF_FONT_SIZE
        curves.ff['flag_graphic'] = False

    # --- create a figure and plot data ---
    fig, ax, css = plot(
        curves, fig, ax,
        FIG_W=curves.ff['figure_width']/2,
        FIG_H=curves.ff['figure_height']/2,
    )

    # --- ivis: interactive visualisation ---
    if curves.ff['flag_ivis']:
        oi = Ivis(fig=fig, curves=curves)

    return fig, ax, css


def plot(curves, fig=None, ax=None, FIG_W=None, FIG_H=None):
    flag_1d = False
    if curves.list_curves[0].zs is None:
        flag_1d = True

    if flag_1d:
        fig, ax, css = cpr.plot_curves(curves, fig, ax, FIG_W, FIG_H)
    else:
        fig, ax, css = cpr.plot_curves_3d(curves, fig, ax, FIG_W, FIG_H)

    return fig, ax, css


class Ivis:
    root = None  # main window
    curves = None
    fig = None

    # styles:
    colorbg_root = mix.to_rgb((120, 120, 120))
    colorbb_upper_frame = mix.to_rgb((140, 140, 140))

    ext_data = '.dat'
    ext_latex = '.tex'
    ext_png = '.png'
    ext_eps = '.eps'

    WINDOW_SIZE_W = 1900
    WINDOW_SIZE_H = 990
    WINDOW_POSITION_FROM_RIGHT = 2
    WINDOW_POSITION_FROM_DOWN = 2

    def __init__(self, fig, curves, **kwargs):
        # Set curves and figure to plot:
        self.curves = curves
        self.fig = fig

        # Main window (top widget)
        self.root = tk.Tk()
        self.root.geometry("{:d}x{:d}+{:d}+{:d}".format(
            int(self.WINDOW_SIZE_W), int(self.WINDOW_SIZE_H),
            int(self.WINDOW_POSITION_FROM_RIGHT), int(self.WINDOW_POSITION_FROM_DOWN)
        ))
        self.root.configure(
            bg=self.colorbg_root
        )
        self.root.wm_title("ivis")

        # Upper frame
        fUpper = ivf.UpperFrame(mw=self, master=self.root, height=30)

        # Figure frame
        figFrame = ivf.FigureFrame(mw=self, master=self.root)

        # Left frame
        fLeft = ivf.LeftFrame(mw=self, master=self.root)

        # # Button: quit
        # bQuit = ivb.QuitButton(
        #     mw=self,
        #     master=self.root,
        #     text="Quit",
        #     command=self.root.destroy,
        # )

        # # Button: save pgfplot
        # bPgf = tk.Button(master=self.root, text="Save pgfplot",
        #                  command=self.on_press_pgfplot_save)

        # --- create the window grid ---
        fUpper.grid(row=0, column=0, columnspan=2, sticky=tk.N + tk.S + tk.E + tk.W)
        fLeft.grid(row=1, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        figFrame.grid(row=1, column=1, sticky=tk.N + tk.S + tk.E + tk.W)
        self.root.rowconfigure(0, weight=0)
        self.root.rowconfigure(1, weight=1)
        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=1)

        # fUpper.grid(row=0, column=0, columnspan=2)
        # fRootProc.grid(row=1, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        # figFrame.grid(row=1, column=1, sticky=tk.N + tk.S + tk.E + tk.W)
        # self.root.rowconfigure(0, weight=0)
        # self.root.rowconfigure(1, weight=1)
        # self.root.columnconfigure(0, weight=1)
        # self.root.columnconfigure(1, weight=1)
        #
        # fFigProp.grid(row=0, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        # fAxProp.grid(row=1, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        # fProc.grid(row=2, column=0, sticky=tk.N + tk.S + tk.E + tk.W)
        # fRootProc.columnconfigure(0, weight=1)
        # fRootProc.rowconfigure(0, weight=1)
        # fRootProc.rowconfigure(1, weight=1)
        # fRootProc.rowconfigure(2, weight=1)

        # --- Launch tkinter loop ---
        tk.mainloop()

    def on_press_pgfplot_save(self):

        # get file name and path to this file from the filedialog
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
        # save files necessary to create a pgfplot

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

        # set file name where data are read from
        result_data_line = ''
        for id_curve, file_data_name in enumerate(file_data_name_curves):
            line_file_name = '\\addplot table [y=Y, x=X]{' \
                             + '{}'.format(file_data_name) \
                             + '};\n'
            line_legend = '\\addlegendentry{' \
                          + '{}'.format(self.curves.list_curves[id_curve].ff['legend']) \
                          + '};\n'
            result_data_line += line_file_name + line_legend
        result_data_line = result_data_line[:-1]

        result_text = self.template_set(
            result_text, 4, result_data_line
        )

        # set title
        resulting_title = '' \
            if self.curves.ff['title'] is None \
            else self.curves.ff['title']
        result_text = self.template_set(result_text, 1, resulting_title)

        # set XY labels
        result_text = self.template_set(
            result_text, 2, self.curves.ff['xlabel']
        )
        result_text = self.template_set(
            result_text, 3, self.curves.ff['ylabel']
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

        # set width of the figure:
        resulting_w = GLO.PGFPLOT_WIDTH
        if self.curves.ff['figure_width'] != GLO.FIG_SIZE_W:
            resulting_w = GLO.PGFPLOT_WIDTH * \
                 self.curves.ff['figure_width'] / self.curves.ff['figure_height']
        result_text = self.template_set(result_text, 7, resulting_w)

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

        # --- rebuild the figure without axes and enlarged ---
        self.curves.ff['flag_graphic'] = True
        self.curves.ff['flag_add_text_plot'] = False
        fig_new, _, _ = plot(self.curves)
        self.curves.ff['flag_graphic'] = False
        self.curves.ff['flag_add_text_plot'] = True

        # --- create external figure .eps (using Image Magic) ---
        name_plot_png = fname + self.ext_png
        fig_new.savefig(name_plot_png)

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
        resulting_title = '' \
            if self.curves.ff['title'] is None \
            else self.curves.ff['title']
        result_text = self.template_set(result_text, 2, resulting_title)

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

        # set width of the figure:
        resulting_w = GLO.PGFPLOT_WIDTH
        if self.curves.ff['figure_width'] != GLO.FIG_SIZE_W:
            resulting_w = GLO.PGFPLOT_WIDTH * \
                          self.curves.ff['figure_width'] / self.curves.ff['figure_height']
        result_text = self.template_set(result_text, 12, resulting_w)

        # add text:
        resulting_command = 'node[]{}'

        ntexts = len(self.curves.list_text)
        if ntexts != 0:
            resulting_command = ''

        for itext in range(ntexts):
            loc_text = self.curves.list_text[itext]
            resulting_command += '\\node[{}] at (axis cs: {:0.3e}, {:0.3e})'\
                .format(loc_text.color, loc_text.x, loc_text.y)
            resulting_command += ' {' + '{}'.format(loc_text.line) + '};\n'
        resulting_command = resulting_command[:-1]
        result_text = self.template_set(result_text, 400, resulting_command)

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





