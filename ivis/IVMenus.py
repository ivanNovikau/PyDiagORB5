import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr

import os
from matplotlib import rcParams
import re
import pandas as pd
import numpy as np
import tkinter as tk
import tkinter.filedialog
import tkinter.messagebox


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


# *** File Menu in main Menu bar ***
class FileMenu(BMenu):

    def __init__(self, mw, **kwargs):
        super(FileMenu, self).__init__(mw, **kwargs)
        self.add_command(label="Save pgfplot", command=self.on_press_pgfplot_save)
        self.add_separator()
        self.add_command(label="Exit", command=self.mw.root.destroy)

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
        # initialfile = self.canvas.get_default_filename()

        fname = tk.filedialog.asksaveasfilename(
            master=self.mw.root,
            title='Save the figure',
            defaultextension=defaultextension,
            initialdir=initialdir,
            # initialfile=initialfile,
        )

        # Save dir for next time, unless empty str (i.e., use cwd).
        if initialdir != "":
            rcParams['savefig.directory'] = (
                os.path.dirname(str(fname))
            )
        return fname

    def save_pgfplot(self, fname):
        # save 1d plot
        if self.mw.curves.list_curves[0].zs is None:
            self.save_pgfplot_1d(fname)
        # save 2d plot
        else:
            self.save_pgfplot_2d(fname)

    def save_pgfplot_1d(self, fname):
        print('Save the 1d plot to files with a name: {}'.format(fname))

        curves = self.mw.curves

        # --- save data to a .dat file ---
        file_data_name_curves = []
        for id_curve, one_curve in enumerate(curves.list_curves):
            file_name = fname + "{:d}".format(id_curve) + GLO.ext_data
            file_data_name_curves.append(file_name)
            output_df = pd.DataFrame({
                'X': one_curve.xs,
                'Y': one_curve.ys
            })
            output_df.to_csv(file_name, sep=" ", index=False)

        # --- read template to create a .tex file ---
        ff_template = open("ivis/template_plot_1d.txt", 'r')
        template_text = ff_template.read()
        ff_template.close()

        # --- save a corresponding .tex file ---
        file_name = fname + GLO.ext_latex
        ff_tex = open(file_name, 'w')
        result_text = template_text

        # set .dat file name
        result_data_line = ''
        for id_curve, file_data_name in enumerate(file_data_name_curves):
            line_file_name = '\\addplot table [y=Y, x=X]{' \
                             + '{}'.format(file_data_name) \
                             + '};\n'
            line_legend = '\\addlegendentry{' \
                          + '{}'.format(curves.list_curves[id_curve].ff['legend']) \
                          + '};\n'
            result_data_line += line_file_name + line_legend
        result_data_line = result_data_line[:-1]
        result_text = mix.template_set(result_text, 1, result_data_line)

        # set general figure properties:
        result_text = self.set_figure_general_properties(result_text)

        # save the resulting text into the .tex file
        ff_tex.write(result_text)
        ff_tex.close()

    def save_pgfplot_2d(self, fname):
        print('Save the 2d plot to files with a name: {}'.format(fname))

        # --- extract necessary data ---
        ref_curve = self.mw.curves.list_curves[0]
        xmin, xmax = np.min(ref_curve.xs), np.max(ref_curve.xs)
        ymin, ymax = np.min(ref_curve.ys), np.max(ref_curve.ys)
        zmin, zmax = np.min(ref_curve.zs), np.max(ref_curve.zs)

        # --- rebuild the figure without axes, and enlarged ---
        self.mw.curves.ff['flag_graphic'] = True
        self.mw.curves.ff['flag_add_text_plot'] = False
        fig_new, _, _ = cpr.plot(self.mw.curves)
        self.mw.curves.ff['flag_graphic'] = False
        self.mw.curves.ff['flag_add_text_plot'] = True

        # --- create external .eps figure (using Image Magic) ---
        name_plot_png = fname + GLO.ext_png
        fig_new.savefig(name_plot_png)

        name_plot_eps = fname + GLO.ext_eps
        command_line = 'convert -trim ' + name_plot_png + ' ' + name_plot_eps
        os.system(command_line)

        # --- read template to create a .tex file ---
        ff_template = open("ivis/template_plot_2d.txt", 'r')
        template_text = ff_template.read()
        ff_template.close()

        # --- save the corresponding .tex file ---
        file_name = fname + GLO.ext_latex
        ff_tex = open(file_name, 'w')
        result_text = template_text

        # set .eps file name to read
        result_text = mix.template_set(result_text, 1, name_plot_eps)

        # set min and max values of axes
        result_text = mix.template_set(result_text, 2, xmin)
        result_text = mix.template_set(result_text, 3, xmax)
        result_text = mix.template_set(result_text, 4, ymin)
        result_text = mix.template_set(result_text, 5, ymax)
        result_text = mix.template_set(result_text, 6, zmin)
        result_text = mix.template_set(result_text, 7, zmax)

        # set colormap
        resulting_colormap = ref_curve.ff['colormap'] \
            if ref_curve.ff['colormap'] is not None \
            else GLO.DEF_COLORMAP
        result_text = mix.template_set(result_text, 8, resulting_colormap)

        # set general figure properties:
        result_text = self.set_figure_general_properties(result_text)

        # save the resulting text into the .tex file
        ff_tex.write(result_text)
        ff_tex.close()

    def set_figure_general_properties(self, result_text):
        # set title (101)
        result_text = self.set_title(result_text)

        # set XY labels (102, 103)
        result_text = self.set_xy_labels(result_text)

        # set figure width (104)
        result_text = self.set_figure_width(result_text)

        # set additional XY ticks (105, 106):
        result_text = self.set_xy_extra_ticks(result_text)

        # set additional text (107):
        result_text = self.set_additional_text(result_text)

        return result_text

    def set_title(self, result_text):
        resulting_title = '' \
            if self.mw.curves.ff['title'] is None \
            else self.mw.curves.ff['title']
        return mix.template_set(result_text, 101, resulting_title)

    def set_xy_labels(self, result_text):
        result_text = mix.template_set(
            result_text, 102, self.mw.curves.ff['xlabel']
        )
        result_text = mix.template_set(
            result_text, 103, self.mw.curves.ff['ylabel']
        )
        return result_text

    def set_figure_width(self, result_text):
        resulting_w = GLO.PGFPLOT_WIDTH
        if self.mw.curves.ff['figure_width'] != GLO.FIG_SIZE_W:
            resulting_w = GLO.PGFPLOT_WIDTH * \
                          self.mw.curves.ff['figure_width'] / self.mw.curves.ff['figure_height']
        return mix.template_set(result_text, 104, resulting_w)

    def set_xy_extra_ticks(self, result_text):
        line_ticks = ''
        if not np.isnan(self.mw.curves.ff['xticks']).any():
            line_ticks = '{}'.format(line_ticks)
            line_ticks = line_ticks[1:-1]
        result_text = mix.template_set(result_text, 105, line_ticks)

        line_ticks = ''
        if not np.isnan(self.mw.curves.ff['yticks']).any():
            line_ticks = '{}'.format(line_ticks)
            line_ticks = line_ticks[1:-1]
        return mix.template_set(result_text, 106, line_ticks)

    def set_additional_text(self, result_text):
        resulting_command = 'node[]{};\n'

        ntexts = len(self.mw.curves.list_text)
        if ntexts != 0:
            resulting_command = ''

        for itext in range(ntexts):
            loc_text = self.mw.curves.list_text[itext]
            resulting_command += '\\node[{}] at (axis cs: {:0.3e}, {:0.3e})' \
                .format(loc_text.color, loc_text.x, loc_text.y)
            resulting_command += ' {' + '{}'.format(loc_text.line) + '};\n'
        resulting_command = resulting_command[:-1]
        return mix.template_set(result_text, 107, resulting_command)



