import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import general as gn
import numpy as np
from scipy import interpolate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(gn)


# NEW:
def choose_one_var_tvpar(ovar, dd):
    opt_var = ovar[0]
    species_name = ovar[1]

    data, t, tit_var, vpar = [], [], [], []
    if opt_var == 'f_vel_1d':
        rd.distribution_1d(dd, species_name)
        data = dd[species_name].f_1d['f_vel_1d']
        tit_var = species_name + ':\ f(v_{\parallel})'
        t = dd[species_name].f_1d['t']
        vpar = dd[species_name].f_1d['vpar']
    if opt_var == 'df_vel_1d-dv':
        rd.distribution_1d(dd, species_name)
        data = dd[species_name].f_1d['f_vel_1d']
        tit_var = species_name + ':\ f(v_{\parallel})'
        t = dd[species_name].f_1d['t']
        vpar = dd[species_name].f_1d['vpar']

        data = np.gradient(data, vpar, axis=1)
    if opt_var == 'df_vel_1d':
        rd.distribution_1d(dd, species_name)
        data = dd[species_name].f_1d['df_vel_1d']
        tit_var = species_name + ':\ \delta f(v_{\parallel})'
        t = dd[species_name].f_1d['t']
        vpar = dd[species_name].f_1d['vpar']

    res = {
        'data': data,
        't':   t,
        'vpar': vpar,
        'tit': tit_var
    }

    return res
