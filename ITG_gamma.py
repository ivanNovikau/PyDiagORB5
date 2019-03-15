from numpy.core.multiarray import ndarray

import Mix as mix
import read_data as rd
import ControlPlot as cpr
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)


def read_data(dd, path_to_folder):
    dd = rd.potsc(dd, path_to_folder)
    return dd


def plot_t(dd, s1, s2, chi1, oo={}):
    # coordinate axes
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    _, ids_s = mix.get_array(s, s1, s2)
    t, ids_t = mix.get_array_oo(oo, t, 't')
    id_chi, _ = mix.find(chi, chi1)

    # non-zonal (at phi = 0) Phi in chosen intervals
    pot_nz_chi = mix.get_slice(dd['potsc']['data'], ids_t, id_chi, ids_s)

    # averaging of the Phi on s
    pot_nz = np.mean(pot_nz_chi, axis=1)

    # plot data:
    cpr.plot_x1(t, pot_nz)


def plot_schi(dd, t1, oo={}):
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference
    r = dd['potsc']['r']  # create a new reference
    z = dd['potsc']['z']  # create a new reference

    id_t1 = np.argmax(t >= t1)
    pot_nz = np.array(dd['potsc']['data'][id_t1, :, :])  # actually copy data

    oop = {'xlabel': 's', 'ylabel': 'chi'}
    cpr.plot_x1x2(s, chi, pot_nz, oop)

    # oop = {'xlabel': 'r', 'ylabel': 'z'}
    # cpr.plot_x1x2(r, z, pot_nz, oop)






