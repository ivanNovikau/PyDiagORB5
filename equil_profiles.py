import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
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


def q_prof(dd):
    rd.q(dd)

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


def nT_profs(dd):
    # output equilibrium temperature from ORB5 as it is
    curves_T = crv.Curves().xlab('s').ylab('T')\
        .tit('species\ temperature\ (ORB5\ output)').set_diff_styles()
    for sp_name in dd['species_names']:
        curves_T.new(sp_name)\
            .XS(dd[sp_name].nT_equil['s'])\
            .YS(dd[sp_name].nT_equil['T'])\
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # normalized temperature
    curves_T.new_tit('species\ norm.\ temperature').ylab('norm.\ T')
    curves_T.flag_norm = True
    cpr.plot_curves(curves_T)

    # temperature in keV
    curves_T = crv.Curves().xlab('s').ylab('T(keV)')\
        .tit('species\ temperature\ in\ keV').set_diff_styles()
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['T_keV']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # output equilibrium density from ORB5 as it is
    curves_n = crv.Curves().xlab('s').ylab('norm.\ n') \
        .tit('species\ density').set_diff_styles()
    for sp_name in dd['species_names']:
        curves_n.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['n']) \
            .leg(sp_name)
    cpr.plot_curves(curves_n)

    # normalized equilibrium density from ORB5
    curves_n = crv.Curves().xlab('s').ylab('norm.\ n') \
        .tit('species\ norm.\ density').set_diff_styles()
    curves_n.flag_norm = True
    for sp_name in dd['species_names']:
        curves_n.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['n']) \
            .leg(sp_name)
    cpr.plot_curves(curves_n)

    # normalized temperature gradient
    curves_T = crv.Curves().xlab('s').ylab('norm.\ \\nabla(T)') \
        .tit('species\ norm.\ temperature\ gradient')\
        .set_diff_styles()
    curves_T.flag_norm = True
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['gradT']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # normalized temperature gradient over temperature
    curves_T = crv.Curves().xlab('s').ylab('norm.\ \\nabla(T)/T') \
        .tit('species\ norm.\ temperature\ gradient\n over\ temperature') \
        .set_diff_styles()
    curves_T.flag_norm = True
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['gradToT']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # normalized density gradient
    curves_T = crv.Curves().xlab('s').ylab('norm.\ \\nabla(n)') \
        .tit('species\ norm.\ density\ gradient') \
        .set_diff_styles()
    curves_T.flag_norm = True
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['gradn']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # normalized log temperature gradient
    curves_T = crv.Curves().xlab('s').ylab('norm.\ \\nabla(log(T))') \
        .tit('species\ norm.\ log.\ temperature\ gradient') \
        .set_diff_styles()
    curves_T.flag_norm = True
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['grad_logT']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # normalized log density gradient
    curves_T = crv.Curves().xlab('s').ylab('norm.\ \\nabla(log(n))') \
        .tit('species\ norm.\ log.\ density\ gradient') \
        .set_diff_styles()
    curves_T.flag_norm = True
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['grad_logn']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)


def n_profs(dd):
    # electron density averaged in space:
    ne_avr_m3_loc = ne_avr_m3(dd)

    # output equilibrium density from ORB5 as it is
    curves_n = crv.Curves().xlab('s').ylab('norm.\ n') \
        .tit('species\ density').set_diff_styles()
    for sp_name in dd['species_names']:
        curves_n.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['n']) \
            .leg(sp_name)
    cpr.plot_curves(curves_n)

    # normalized equilibrium density from ORB5
    curves_n = crv.Curves().xlab('s').ylab('norm.\ n') \
        .tit('species\ norm.\ density').set_diff_styles()
    curves_n.flag_norm = True
    for sp_name in dd['species_names']:
        curves_n.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['n']) \
            .leg(sp_name)
    cpr.plot_curves(curves_n)

    # density in 1/m3:
    curves_n = crv.Curves().xlab('s').ylab('n(m^{-3})') \
        .tit('species\ density\ (m^{-3})').set_diff_styles()
    for sp_name in dd['species_names']:
        curves_n.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['n'] * ne_avr_m3_loc) \
            .leg(sp_name)
    cpr.plot_curves(curves_n)

    # normalized density gradient
    curves_T = crv.Curves().xlab('s').ylab('norm.\ \\nabla(n)') \
        .tit('species\ norm.\ density\ gradient') \
        .set_diff_styles()
    curves_T.flag_norm = True
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['gradn']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)


def vp_profs(dd):
    # output equilibrium temperature from ORB5 as it is
    curves_vp = crv.Curves().xlab('s').ylab('v_{\parallel}')\
        .tit('species\ parallel\ velocity\ (ORB5\ output)').set_diff_styles()
    for sp_name in dd['species_names']:
        curves_vp.new(sp_name)\
            .XS(dd[sp_name].nT_equil['s'])\
            .YS(dd[sp_name].nT_equil['vp'])\
            .leg(sp_name)
    cpr.plot_curves(curves_vp)


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


# NEW: electron density averaged in space:
def ne_avr_m3(dd):
    B0_gauss     = dd['B0'] * 1e4
    Te_speak_erg = dd['electrons'].T_speak(dd) * 1e7

    res = 1e6 * dd['beta'] * B0_gauss ** 2 / (4 * np.pi * Te_speak_erg)
    return res


# NEW: species density (in 1/m3) at a particular point:
def n_m3_s1(dd, name_sp, s1):
    s = dd[name_sp].nT_equil['s']
    n = dd[name_sp].nT_equil['n']

    id_s1, s1_res, _ = mix.get_ids(s, s1)
    n_res = n[id_s1] * ne_avr_m3(dd)
    return n_res, s1_res


# NEW: choose a variable (s):
def choose_one_var_ts(ovar, dd):
    opt_var = ovar[0]
    tit_var = ''
    vvar, s, t = [], None, None
    if opt_var == 'q':
        rd.q(dd)
        var_name = 'q'
        tit_var = 'q'
        vvar.append(dd[var_name]['data'])
        s = dd[var_name]['s']
        t = np.array([0])

    if opt_var == 'T-equ':
        sp_name = ovar[1]
        tit_var = sp_name + ':\ T'
        vvar.append(dd[sp_name].nT_equil['T'])
        s = dd[sp_name].nT_equil['s']
        t = np.array([0])

    if opt_var == 'T-equ-keV':
        sp_name = ovar[1]
        tit_var = sp_name + ':\ T(keV)'
        vvar.append(dd[sp_name].nT_equil['T_keV'])
        s = dd[sp_name].nT_equil['s']
        t = np.array([0])

    if opt_var == 'n-equ':
        sp_name = ovar[1]
        tit_var = sp_name + ':\ n'
        vvar.append(dd[sp_name].nT_equil['n'])
        s = dd[sp_name].nT_equil['s']
        t = np.array([0])

    if opt_var == 'n-equ-m3':
        sp_name = ovar[1]
        tit_var = sp_name + ':\ n'
        vvar.append(dd[sp_name].nT_equil['n'] * ne_avr_m3(dd))
        s = dd[sp_name].nT_equil['s']
        t = np.array([0])

    res = {
        'data': np.array(vvar),
        's': s,
        't': t,
        'tit': tit_var
    }
    return res


# NEW: get a safety factor value at a radial position s1:
def q_s1(dd, s1):
    rd.q(dd)

    s = dd['q']['s']
    q = dd['q']['data']

    id_s1, s1_new, line_s1 = mix.get_ids(s, s1, '{:0.2f}')
    value_q_s1 = q[id_s1]

    return value_q_s1, s1_new, line_s1
