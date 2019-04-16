import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import zf_gam as zf
import ITG_gamma as itg
import write_data as wr
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(zf)
    mix.reload_module(itg)
    mix.reload_module(wr)


def save_aug_signals(dd, oo):
    # radial points, where the signals will be considered
    s_points = oo.get('s_points', [0.5])

    # zonal Er:
    zf.erbar(dd)

    # poloidal angles, where a full signal will be considered
    chi_s = [0.0,           np.pi/4,           3*np.pi/8,   np.pi/2,
                  2*np.pi - np.pi/4, 2*np.pi - 3*np.pi/8, 3*np.pi/2]

    # full Er:
    names_ersc = itg.ersc(dd, {'chi_s': chi_s})

    # open .h5 file
    ff = wr.create_open_file(dd)

    # save radial points to the result file:
    wr.save_data('s_points', s_points)
    
    # close .h5 file
    wr.close_file(ff)
