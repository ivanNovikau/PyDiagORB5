import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np
from scipy import interpolate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)


def exp_AUG20787(dd, oo):
    # curves, where to add the data
    curves = oo.get('curves', crv.Curves())
    sel_norm = oo.get('sel_norm', 'wci')  # -> 'wci', 'khz', 'csa', 'csr'
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'
    col = oo.get('col', 'blue')

    # experimental GAM frequency:
    rho_fr = np.array([0.863, 0.879, 0.891, 0.902, 0.912, 0.912,
              0.922, 0.922, 0.932, 0.932, 0.941, 0.950,
              0.959, 0.967, 0.975, 0.984])
    rho_err = 0.005 * np.ones(np.size(rho_fr))
    fr_GAM_kHz = np.array([20.2, 20.1, 18.9, 18.9, 17.7, 15.5, 16.5,
                  14.0, 14.5, 13.4, 14.0, 13.4, 13.2, 12.2,
                  12.2, 12.5])  # kHz
    fr_err = 0.4 * np.ones(np.size(rho_fr))

    # experimental profiles:
    rho_T_AUG = np.array([0.863, 0.879, 0.891, 0.902, 0.912,
                 0.922, 0.932, 0.941, 0.950, 0.959,
                 0.967, 0.975, 0.984])
    Te_AUG = np.array([240,   210,   185,   165,   140,
              122,   105,    85,   70,     52,
               40,    25,    15])  # eV
    Ti_AUG = np.array([240,   210,   185,   165,   140,
              122,   105,    85,    70,    57,
               50,    38,    30])  # eV

    # choose normalization:
    coef_norm = None
    if sel_norm == 'khz':
        coef_norm = 1
    if sel_norm == 'wc':
        coef_norm = 1.e3 * 2 * np.pi / dd['wc']
    if sel_norm == 'csa':
        coef_norm = 1.e3 * 2 * np.pi / (dd['cs']/dd['a0'])
    if sel_norm == 'csr':
        coef_norm = 1.e3 * 2 * np.pi / (dd['cs']/dd['R0'])
    fr_res     = fr_GAM_kHz * coef_norm
    fr_err_res = fr_err     * coef_norm

    r, r_err = [], []
    if sel_r == 's':
        r = np.sqrt(rho_fr)
        r_err = np.zeros([2, np.size(rho_fr)])
        r_err[0] = np.sqrt(rho_fr) - np.sqrt(rho_fr - rho_err)
        r_err[1] = np.sqrt(rho_fr + rho_err) - np.sqrt(rho_fr)
    if sel_r == 'psi':
        r = rho_fr
        r_err = rho_err

    curves.new('aug20787').XS(r).YS(fr_res)\
        .XS_ERR(r_err).YS_ERR(fr_err_res).leg('EXP.\ AUG20787')\
        .sty('o').col(col).ms(1)

    return curves
