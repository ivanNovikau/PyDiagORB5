import Mix as mix
import ControlPlot as cpr
import ymath
import curve as crv
import Global_variables as GLO
import numpy as np
import h5py as h5


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(GLO)


def read_q(dd):
    if 'q' in dd:
        return
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')
    s = np.array(f['/equil/profiles/generic/sgrid_eq'])
    data = np.array(f['/equil/profiles/generic/q'])

    dd['q'] = {
        's': s,
        'data': data
    }


def choose_one_var_ts(one_signal):
    dd = one_signal['dd']
    opt_var = one_signal['variable']
    if opt_var == 'q':
        read_q(dd)
        var_name = 'q'
        tit_var = 'q'
        vvar = dd[var_name]['data']
        s = dd[var_name]['s']
        t = np.array([0])
    elif opt_var == 'T':
        sp_name = one_signal['species']
        tit_var = sp_name + ':\ T'
        vvar = dd[sp_name].nT_equil['T']
        s = dd[sp_name].nT_equil['s']
        t = np.array([0])
    elif opt_var == 'kT':
        sp_name = one_signal['species']
        tit_var = sp_name + ':\ T'
        vvar = dd[sp_name].nT_equil['T']
        s = dd[sp_name].nT_equil['s']
        vvar = - np.gradient(vvar, s)/vvar
        t = np.array([0])
    elif opt_var == 'T-keV':
        sp_name = one_signal['species']
        tit_var = sp_name + ':\ T(keV)'
        vvar = dd[sp_name].nT_equil['T_keV']
        s = dd[sp_name].nT_equil['s']
        t = np.array([0])
    elif opt_var == 'n':
        sp_name = one_signal['species']
        tit_var = sp_name + ':\ n'
        vvar = dd[sp_name].nT_equil['n']
        s = dd[sp_name].nT_equil['s']
        t = np.array([0])
    else:
        mix.error_mes('Wrong name of an equilibrium variable.')

    res = {
        'data': np.array([vvar]),
        's': s,
        't': t,
        'tit': tit_var
    }
    return res


def q_adhoc(dd, qcoefs):
    # qcoefs = np.array([1.8503, -0.074096, 0.95318, -5.703, 7.5779, -0.19182, -1.8253])

    ncoef = np.size(qcoefs)
    qa = 0
    rho = dd['deuterium'].nT_equil['rho']
    for icoef in range(ncoef):
        qa += qcoefs[icoef] * rho**icoef

    return qa


# Get electron density averaged in space:
def ne_avr_m3(dd):
    B0_gauss     = dd['B0'] * 1e4
    Te_speak_erg = dd['electrons'].T_speak(dd) * 1e7

    res = 1e6 * dd['beta'] * B0_gauss ** 2 / (4 * np.pi * Te_speak_erg)
    return res
