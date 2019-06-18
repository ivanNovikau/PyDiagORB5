import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import write_data as wr
import general
import gam_theory
import Geom as geom
import common
import numpy as np
from scipy import interpolate
import h5py as h5
from scipy.signal import correlate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(wr)
    mix.reload_module(general)
    mix.reload_module(gam_theory)
    mix.reload_module(geom)
    mix.reload_module(common)


# Calculate frequency (w) and damping/growth rate(g)
def calc_wg(oo_var, oo_wg):
    # --- MAIN FUNCTIONALITY ---
    # -> oo_wg.flag_E = True:
    #   calc. wg using peaks -> wgE
    # -> oo_wg.flag_A = True:
    #   calc. wg using nl-fitting -> wgA
    # -> oo_wg.flag_stat = True: (some printing flags will be switched off)
    #   vary time intervals in a give time range;
    #   compare wgE and wgA, takes only those pairs that are close enough
    #       to get rid of outliers, save those pairs;
    #   build histograms for w and g, and using Gaussian or Student's
    #       distribution functions find 95% confidence intervals as 1.96 sigma
    # -> oo_wg.flag_two_steps:
    #   first step  - find gamma;
    #   second step - find frequency of signal_init * exp(-g*t).
    # -----------------------------------------------------------------------------
    # --- REMARKS ---
    # -> does not work if there are less than three peaks;
    # -----------------------------------------------------------------------------
    # --- INPUT DICTIONARIES ---
    # oo_var - to choose variable (t)
    # oo_wg  - to adjust parameters of the wg calculation

    return