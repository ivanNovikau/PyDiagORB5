import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import gam_theory
import numpy as np
from scipy import interpolate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(gam_theory)


def erbar(dd):
    if 'erbar' in dd:
        return
    rd.phibar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']
    data = - np.gradient(dd['phibar']['data'], s, axis=1)
    dd['erbar'] = {
        't': t,
        's': s,
        'data': data}


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
        'data': data}


def choose_var(dd, oo):
    vorbar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']

    opt_var = oo.get('opt_var', 'erbar')
    vvar = dd[opt_var]['data']

    tit_var = ''
    if opt_var == 'phibar':
        tit_var = '\overline{\Phi}'
    if opt_var == 'erbar':
        tit_var = '\overline{E}_r'
    if opt_var == 'vorbar':
        tit_var = '\overline{\Omega}_r'

    res = {
        'var': vvar,
        's': s,
        't': t,
        'tit': tit_var
    }

    return res


# NEW: take signal (t,s)
def choose_one_var_ts(ovar, dd):
    vorbar(dd)
    t = dd['phibar']['t']
    s = dd['phibar']['s']

    opt_var = ovar[0]

    vvar = dd[opt_var]['data']

    tit_var = ''
    if opt_var == 'phibar':
        tit_var = '\overline{\Phi}'
    if opt_var == 'erbar':
        tit_var = '\overline{E}_r'
    if opt_var == 'vorbar':
        tit_var = '\overline{\Omega}_r'

    res = {
        'data': vvar,
        's': s,
        't': t,
        'tit': tit_var
    }

    return res

















