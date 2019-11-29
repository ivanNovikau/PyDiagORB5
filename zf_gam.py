import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import gam_theory
import numpy as np
from scipy import interpolate
import h5py as h5


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(gam_theory)


def phibar(dd):
    if 'phibar' in dd:
        return
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')

    t = np.array(f['/data/var1d/generic/phibar/time'])
    s = np.array(f['/data/var1d/generic/phibar/coord1'])
    data = np.array(f['/data/var1d/generic/phibar/data'])
    dd['phibar'] = {
        't': t,
        's': s,
        'data': data
    }


def erbar(dd):
    if 'erbar' in dd:
        return
    phibar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']
    data = - np.gradient(dd['phibar']['data'], s, axis=1)
    dd['erbar'] = {
        't': t,
        's': s,
        'data': data
    }


def vorbar(dd):
    # shearing rate, vorticity
    if 'vorbar' in dd:
        return
    erbar(dd)
    t = dd['erbar']['t']
    s = dd['erbar']['s']

    data = dd['erbar']['data'] * s
    data = np.gradient(data, s, axis=1)
    if s[0] > 0:
        data = data / s
    else:
        data[:, 1:-1] = data[:, 1:-1] / s[1:-1]
        f = interpolate.interp2d(s[1:-1], t, data[:, 1:-1], kind='cubic')
        data = f(s, t)

    dd['vorbar'] = {
        't': t,
        's': s,
        'data': data
    }


def choose_one_var_ts(one_signal):
    dd = one_signal['dd']
    opt_var = one_signal['variable']

    vorbar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']
    if opt_var == 'phi':
        var_name = 'phibar'
        tit_var = '\overline{\Phi}'
    elif opt_var == 'er':
        var_name = 'erbar'
        tit_var = '\overline{E}_r'
    elif opt_var == 'vorticity':
        var_name = 'vorbar'
        tit_var = '\overline{\Omega}_r'
    else:
        mix.error_mes('Wrong name of a zonal signal.')

    vvar = dd[var_name]['data']
    res = {
        'data': vvar,
        's': s,
        't': t,
        'tit': tit_var
    }

    return res

















