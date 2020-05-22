import Mix as mix
import ymath
import curve as crv
import Global_variables as GLO
import ControlPlot as cpr
import ivis.IVButtons as ivb
import ivis.IVFrames as ivf
import ivis.IVMenus as ivm

import tkinter as tk


def reload():
    mix.reload_module(mix)
    mix.reload_module(crv)
    mix.reload_module(ymath)
    mix.reload_module(GLO)
    mix.reload_module(cpr)
    mix.reload_module(ivb)
    mix.reload_module(ivf)
    mix.reload_module(ivm)


def plot_data(curves, fig=None, ax=None):
    # --- modify font size if ivis is used ---
    if curves.ff['flag_ivis']:
        curves.ff['fontS'] /= GLO.IVIS_COEF_FONT_SIZE
        curves.ff['flag_graphic'] = False

    # --- create a figure and plot data ---
    fig, ax, css = cpr.plot(
        curves, fig, ax,
        FIG_W=curves.ff['figure_width']/2,
        FIG_H=curves.ff['figure_height']/2,
    )

    # --- ivis: interactive visualisation ---
    if curves.ff['flag_ivis']:
        oi = Ivis(fig=fig, curves=curves)

    return fig, ax, css


class Ivis:
    root = None  # main window
    curves = None
    fig = None
    flag_2d = False
    curves_default = None

    # Menu bar
    menubar = None

    # styles:
    colorbg_root = mix.to_rgb((120, 120, 120))

    # frames:
    fLeft = None
    fFigure = None

    WINDOW_SIZE_W = 1900
    WINDOW_SIZE_H = 990
    WINDOW_POSITION_FROM_RIGHT = 2
    WINDOW_POSITION_FROM_DOWN = 2

    # --- Available shortcurts ---
    # Ctrl + S (Ctrl + Shift + s)  - save pgfplot

    def __init__(self, fig, curves, **kwargs):
        # Set curves and figure to plot:
        self.curves = curves
        self.fig = fig
        if curves.list_curves[0].zs is not None:
            self.flag_2d = True

        ax = self.get_ax()
        self.curves_default = crv.copy_curves(self.curves, ax)

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

        # to handle closure of the window
        self.root.protocol("WM_DELETE_WINDOW", self.delete_window)
        self.root.bind("<Destroy>", self.destroy)

        # Figure frame
        self.fFigure = ivf.FigureFrame(mw=self, master=self.root)

        # Left frame
        self.fLeft = ivf.LeftFrame(mw=self, master=self.root)

        # create Menu bar
        self.menubar = tk.Menu(self.root)

        # create a pulldown menu, and add it to the menu bar
        filemenu = ivm.FileMenu(
            mw=self,
            master=self.menubar,
            tearoff=0
        )
        self.menubar.add_cascade(label="File", menu=filemenu)
        self.root.config(menu=self.menubar)

        # --- create the window grid ---
        n_columns = 3
        self.fLeft.grid(
            row=0, column=0,
            sticky=tk.N + tk.S + tk.E + tk.W
        )
        self.fFigure.grid(
            row=0, column=1, columnspan=n_columns-1,
            sticky=tk.N + tk.S + tk.E + tk.W
        )
        self.root.rowconfigure(0, weight=1)
        for id_column in range(n_columns):
            self.root.columnconfigure(id_column, weight=1)

        # --- Launch tkinter loop ---
        tk.mainloop()

    def delete_window(self):
        print("delete_window")
        try:
            self.root.destroy()
        except:
            pass

    def destroy(self, event):
        pass
        # print("destroy application")

    def get_ax(self):
        return self.fig.axes[0]

    def draw(self):
        self.fig.canvas.draw()
