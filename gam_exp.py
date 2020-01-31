import Global_variables as GLO
import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(GLO)


def exp_AUG20787(dd, oo):
    # Some experimental data from AUG shot #20787,
    # described in [Conway 2008 Plasma Physics and Controlled Fusion]
    # https://iopscience.iop.org/article/10.1088/0741-3335/50/5/055009/pdf
    # ------------------------------------------------------------------------------------------
    # CONWAY, G.D. et al., “Frequency scaling and localization of geodesic acoustic modes in
    # ASDEX Upgrade”, Plasma Phys. Control. Fusion 50 (2008) 055009.
    # ------------------------------------------------------------------------------------------

    # curves, where to add the data
    curves = oo.get('curves', crv.Curves())
    sel_norm = oo.get('sel_norm', 'wc')  # -> 'wci', 'khz', 'csa', 'csr'
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'

    # experimental GAM frequency:
    rho_fr = np.array([0.863, 0.879, 0.891, 0.902, 0.912, 0.912,
              0.922, 0.922, 0.932, 0.932, 0.941, 0.950,
              0.959, 0.967, 0.975, 0.984])
    rho_err = 0.005 * np.ones(np.size(rho_fr))
    fr_GAM_kHz = np.array([20.2, 20.1, 18.9, 18.9, 17.7, 15.5, 16.5,
                  14.0, 14.5, 13.4, 14.0, 13.4, 13.2, 12.2,
                  12.2, 12.5])  # kHz
    fr_err_kHz = 0.4 * np.ones(np.size(rho_fr))

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

    # to wc normalization:
    norm_to_wc = 1.e3 * 2 * np.pi / dd['wc']
    fr_wc     = fr_GAM_kHz * norm_to_wc
    fr_err_wc = fr_err_kHz * norm_to_wc

    # choose a radial grid
    r, r_err, desc_r = [], [], None
    if sel_r == 's':
        r = np.sqrt(rho_fr)
        r_err = np.zeros([2, np.size(rho_fr)])
        r_err[0] = np.sqrt(rho_fr) - np.sqrt(rho_fr - rho_err)
        r_err[1] = np.sqrt(rho_fr + rho_err) - np.sqrt(rho_fr)
        desc_r = '(s)'
    elif sel_r == 'psi':
        r = rho_fr
        r_err = rho_err
        desc_r = '(psi)'
    else:
        mix.error_mes('Wrong selector of the radial grid normalization.')

    # normalization
    dict_norm = mix.normalization('frequency-' + sel_norm, dd)
    coef_norm = dict_norm['coef_norm']
    desc_norm = dict_norm['line_norm']

    fr_res     = fr_wc * coef_norm
    fr_err_res = fr_err_wc * coef_norm

    # description of the data
    print('Exp. AUG20787: w' + desc_r + desc_norm)

    # styling:
    ff = dict(GLO.DEF_CURVE_FORMAT)
    ff.update({
        'legend': 'Experiment:\ [Conway08\ PPCF]',
    })

    # curve format
    curves.new().XS(r).YS(fr_res)\
        .set_ff(ff)\
        .set_errorbar(True, ys=fr_err_res, xs=r_err)

    return curves


def exp_AUG31213():
    # Approximate experimental EGAM spectrogram:
    exp_w_kHZ = [48, 50.5, 52.25, 53.75, 55, 56, 56.8, 57.4]
    exp_t_s = [0.8410, 0.8415, 0.8420, 0.8425, 0.8430, 0.8435, 0.8440, 0.8445]
    exp_t_ms = (np.array(exp_t_s) - exp_t_s[0]) * 1e3
    exp_w_rel = np.array(exp_w_kHZ) / exp_w_kHZ[0]

    exp_t_ms_cont = np.linspace(exp_t_ms[0], exp_t_ms[-1], 101)

    exp_w_abs_kHz_cont = np.interp(exp_t_ms_cont, exp_t_ms, np.array(exp_w_kHZ))
    exp_w_rel_cont = np.interp(exp_t_ms_cont, exp_t_ms, exp_w_rel)

    res = {
        't_ms': exp_t_ms_cont,
        'w_rel': exp_w_rel_cont,
        'w_kHz': exp_w_abs_kHz_cont,
    }
    return res
