import Mix as mix
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np
import h5py as h5


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)


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


def q_prof(dd):
    read_q(dd)

    s = dd['q']['s']
    # rho = dd['deuterium'].nT_equil['rho']
    #
    # # analytical expression:
    # qcoefs = np.array([1.8503, -0.074096, 0.95318, -5.703, 7.5779, -0.19182, -1.8253])
    # ncoef = np.size(qcoefs)
    # qa = 0
    # for icoef in range(ncoef):
    #     qa += qcoefs[icoef] * rho**icoef

    curves = crv.Curves().xlab('s').ylab('q').tit('safety\ factor\ profile')
    curves.new().XS(s).YS(dd['q']['data']).leg('q')
    # curves.new().XS(s).YS(qa).leg('q-analytical')
    cpr.plot_curves(curves)


def B_equil(dd):
    # --- read name of the equilibrium .h5 file ---
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')

    equ_file = f['/parameters/equil/fname'].attrs
    ids_attr = list(equ_file)
    equ_file = [equ_file[name].decode("utf-8")
                           for name in ids_attr[0:len(ids_attr)]]
    equ_file = equ_file[0]

    f.close()

    # --- read parameters from the equilibrium file ---
    path_to_equ_file = dd['path'] + '/' + equ_file
    ff = h5.File(path_to_equ_file, 'r')

    chi = np.array(ff['/data/grid/CHI'])
    psi = np.array(ff['/data/grid/PSI'])
    B = np.array(ff['/data/var2d/B'])
    R = np.array(ff['/data/var2d/R'])
    Z = np.array(ff['/data/var2d/Z'])

    chi = np.append(chi, 2 * np.pi)

    nChi = np.shape(B)[0]
    nR = np.shape(B)[1]
    B_new = np.zeros([np.size(chi), np.size(psi)])
    Z_new = np.zeros([np.size(chi), np.size(psi)])
    R_new = np.zeros([np.size(chi), np.size(psi)])
    for i in range(nR):
        B_new[:, i] = np.append(B[:, i], B[0, i])
        Z_new[:, i] = np.append(Z[:, i], Z[0, i])
        R_new[:, i] = np.append(R[:, i], R[0, i])

    R_new *= dd['R0']
    B_new *= dd['B0']

    labx = 'R(m)'
    laby = 'Z(m)'
    tit = '|B|(T)'

    curves = crv.Curves().xlab(labx).ylab(laby).tit(tit)
    curves.new().XS(R_new).YS(Z_new).ZS(B_new.T).lev(160)
    cpr.plot_curves_3d(curves)

    f.close()


# Get electron density averaged in space:
def ne_avr_m3(dd):
    B0_gauss     = dd['B0'] * 1e4
    Te_speak_erg = dd['electrons'].T_speak(dd) * 1e7

    res = 1e6 * dd['beta'] * B0_gauss ** 2 / (4 * np.pi * Te_speak_erg)
    return res
