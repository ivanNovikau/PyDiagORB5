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


def chi(species_name, dd):
    if 'chi' in dd:
        return

    # read heat flux
    rd.radial_heat_flux(dd)
    rd.nT_evol(dd, species_name)
    efluxw_rad = dd[species_name].efluxw_rad
    Lx = dd['Lx']
    rho_star_inv = Lx / 2.

    # signals:
    t        = efluxw_rad['t']  # the same for all 1d signals
    s_flux   = efluxw_rad['s']
    rad_flux = efluxw_rad['data']
    s = dd[species_name].nT_evol['s']
    n = dd[species_name].nT_evol['n']
    gradT = dd[species_name].nT_evol['gradT']  # grad is taken w.r.t s not rho

    T0    = dd[species_name].nT_equil['T']
    n0    = dd[species_name].nT_equil['n']
    s_equ = dd[species_name].nT_equil['s']

    func_interp = interpolate.interp1d(s_flux, rad_flux, axis=1)
    rad_flux    = func_interp(s)

    func_interp = interpolate.interp1d(s_equ, T0, fill_value='extrapolate')
    T0          = func_interp(s)

    func_interp = interpolate.interp1d(s_equ, n0, fill_value='extrapolate')
    n0          = func_interp(s)

    gradT0   = np.gradient(T0, s)

    chi_var  = - rad_flux / (n * gradT)
    chi_var0 = - rad_flux / (n0[None, :] * gradT0[None, :])

    chi_norm  = chi_var  * rho_star_inv**2
    chi_norm0 = chi_var0 * rho_star_inv ** 2

    dd[species_name].chi = {
        'data': chi_var,
        'data0': chi_var0,
        'data_norm': chi_norm,
        'data_norm0': chi_norm0,
        's': s,
        't': t,
    }


# NEW:
def choose_one_var_ts(ovar, dd):
    opt_var = ovar[0]
    species_name = ovar[1]

    data, t, s, tit_var = [], [], [], []
    if opt_var == 'chi_norm':
        chi(species_name, dd)
        data = dd[species_name].chi['data_norm']
        tit_var = species_name + ':\ \chi/\chi_B'
        t = dd[species_name].chi['t']
        s = dd[species_name].chi['s']
    if opt_var == 'chi_norm0':
        chi(species_name, dd)
        data = dd[species_name].chi['data_norm0']
        tit_var = species_name + ':\ \chi_0/\chi_B'
        t = dd[species_name].chi['t']
        s = dd[species_name].chi['s']
    if opt_var == 'efluxw_rad':
        rd.radial_heat_flux(dd)
        efluxw_rad = dd[species_name].efluxw_rad
        data = efluxw_rad['data']
        tit_var = species_name + ':\ efluxw\_rad'
        t = efluxw_rad['t']
        s = efluxw_rad['s']
    if opt_var == 'T':
        rd.nT_evol(dd, species_name)
        data = dd[species_name].nT_evol['T']
        s = dd[species_name].nT_evol['s']
        t = dd[species_name].nT_evol['t']
        tit_var = species_name + ':\ T'
    if opt_var == 'n':
        rd.nT_evol(dd, species_name)
        data = dd[species_name].nT_evol['n']
        s = dd[species_name].nT_evol['s']
        t = dd[species_name].nT_evol['t']
        tit_var = species_name + ':\ n'
    if opt_var == 'p':
        rd.nT_evol(dd, species_name)
        data = dd[species_name].nT_evol['p']
        s = dd[species_name].nT_evol['s']
        t = dd[species_name].nT_evol['t']
        tit_var = species_name + ':\ p'

    res = {
        'data': data,
        's': s,
        't': t,
        'tit': tit_var
    }

    return res


def choose_var(dd, oo):
    species_name = oo.get('species_name', 'deuterium')

    opt_var = oo.get('opt_var', 'chi_norm')
    vvar, t, s, tit_var = [], [], [], []
    if opt_var == 'chi_norm':
        chi(species_name, dd)
        vvar = dd[species_name].chi['data_norm']
        tit_var = species_name + ':\ \chi/\chi_B'
        t = dd[species_name].chi['t']
        s = dd[species_name].chi['s']
    if opt_var == 'chi_norm0':
        chi(species_name, dd)
        vvar = dd[species_name].chi['data_norm0']
        tit_var = species_name + ':\ \chi_0/\chi_B'
        t = dd[species_name].chi['t']
        s = dd[species_name].chi['s']
    if opt_var == 'efluxw_rad':
        rd.radial_heat_flux(dd)
        efluxw_rad = dd[species_name].efluxw_rad
        vvar = efluxw_rad['data']
        tit_var = species_name + ':\ efluxw\_rad'
        t = efluxw_rad['t']
        s = efluxw_rad['s']

    res = {
        'var': vvar,
        's': s,
        't': t,
        'tit': tit_var
    }

    return res


def plot_fft(dd, oo):
    out = choose_var(dd, oo)
    vvar, r, t, tit_var = out['var'], out['s'], out['t'], out['tit']

    # radial coordinate normalization
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'
    line_r = ''
    if sel_r == 's':
        r = r
        line_r = 's = \sqrt{\psi/\psi_{edge}}'
    if sel_r == 'psi':
        r = r ** 2
        line_r = '\psi/\psi_{edge}'

    # plotting
    oo_fft = dict(oo)
    oo_fft.update({
        'var': vvar, 't_wci': t, 'r': r,
        'tit': tit_var,
        'labr': line_r
    })
    gn.plot_fft(dd, oo_fft)


def plot_fft_1d(dd, oo):
    out = choose_var(dd, oo)
    vvar, s, t, tit_var = out['var'], out['s'], out['t'], out['tit']

    # plotting
    oo_fft = dict(oo)
    oo_fft.update({
        'vars': [vvar], 'ts_wci': [t], 'ss': [s],
        'tit': tit_var,
        'vars_names': ['']
    })
    gn.plot_fft_1d(dd, oo_fft)
