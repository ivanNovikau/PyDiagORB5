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


def omegagamma_GAM(q, tau_e, khat, elong):
    # ****************************************************
    # *** Author of the script: Alessandro Biancalani. ***
    # ****************************************************
    #
    # Normalization to sqrt(2)*vti / R0
    #
    # Here the GAM frequency and damping rate obtained from
    # Sugama-2008, Qiu-2009, Gao-2009 and Gao-2010 are calculated.
    # For this version, FOW and elongation effects only are kept,
    # but s_k, Delta', epsilon effects are neglected w.r.t. the papers of Gao.
    # q - safety factor;
    # tau_e = Te/Ti;
    # khat = k = kr * rho_i;
    # elong - elongation.

    T = 7/4 + tau_e

    # --- Sugama 2008 ---
    fhat = (1 + 2 * (23+16*tau_e+4*tau_e**2)
            / q ** 2 / (7 + 4 * tau_e) ** 2)
    # [np.sqrt(2)*vti/R]   where v_ti = np.sqrt(T_i/m_i)
    om_SW = np.sqrt(T) * np.sqrt(fhat)

    # It is omega_gam_hat defined like in Sugama - Watanabe
    om_SW_hat = om_SW * q
    num_main = np.exp(-om_SW_hat**2) * \
               (om_SW_hat**4 + (1+2*tau_e)*om_SW_hat**2)
    num_flr = 0.25 * khat ** 2 * q ** 2 * np.exp(-om_SW_hat ** 2 / 4)\
              * (om_SW_hat ** 6 / 128 +
              (1+tau_e) * om_SW_hat ** 4 / 16 + (3/8 + 7*tau_e/16 +
              5*tau_e**2/32) * om_SW_hat ** 2)
    ga_SW = np.sqrt(np.pi) / 2 * q * fhat ** (-1) * (num_main + num_flr)

    # --- Qiu 2009 ---
    f1_Qiu = 31/16 + 9*tau_e/4 + tau_e**2
    f2_Qiu = 747/32 + 481*tau_e/32 + 35*tau_e**2/8 + tau_e**3/2

    # there was a typo in the paper here
    f3_Qiu = 23/8 + 2*tau_e + tau_e**2 / 2

    om_Qiu = np.sqrt(T) * (1 - khat ** 2 / 4 * f1_Qiu / T
                           + khat ** 2 / 4 * f2_Qiu / T ** 2
                           + 1 / 2 / q ** 2 * f3_Qiu / T ** 2)
    g_Qiu  = 1 + 1 / 24 * om_Qiu / q ** 2 / khat ** 3 \
             * (-4 + om_Qiu / khat) + khat * tau_e / om_Qiu + \
              (1 + 5*tau_e/4 + tau_e**2)*khat**2/om_Qiu**2 - khat**2
    ga_Qiu = om_Qiu * np.sqrt(2)/khat**2 * np.exp(-om_Qiu/khat) \
             * (1 + khat**2 / 2 / om_Qiu**2 * f1_Qiu -
            khat**2 / om_Qiu**4 * f2_Qiu - 2 / q**2 / om_Qiu**4
                * f3_Qiu)*g_Qiu

    # --- Gao-2010 ---
    sfunct4 = 0  # to be modified
    epsfunct1 = 0  # to be modified
    deltafunct1 = 0  # to be modified
    om_gao_NEW_funct1 = 1  # to be modified

    # this is the Q function
    om_gao_NEW_funct2 = 1/T**2 * (elong**2+1)/2*(747/32
                + 481*tau_e/32 + 35*tau_e**2/8 + tau_e**3/2) -\
        1/T * ((39+13*elong**2)/16 + (25-elong**2)*tau_e/8
                + (5-elong**2) * tau_e ** 2 / 4) + \
         ((9+6*elong**2+9*elong**4)/16/elong**2/(elong**2+1)
           + (elong**2-1)**2 * tau_e / 4 / elong**2/(elong**2+1))

    om_gao_2010 = np.sqrt(T)*np.sqrt(2)/np.sqrt(elong**2+1)*om_gao_NEW_funct1 \
                  * (1 + (elong**2+1)/(4*q**2) *
                     (23/8 + 2*tau_e + tau_e**2/2)/T**2
                     + khat**2/(4*elong**2)*om_gao_NEW_funct2)
    ga_gao_NEW_funct1 = (1 + sfunct4 + epsfunct1 + deltafunct1)
    ga_gao_NEW_funct2 = 1  # to be modified
    ga_gao_2010 = 4*elong**2*np.sqrt(T)/khat**2/(elong**2+1)**(3/2) \
                  * ga_gao_NEW_funct1 * np.exp(- np.sqrt(T)/khat
                * np.sqrt(2*elong**2/(elong**2+1)) * ga_gao_NEW_funct2)

    # --- RESULTS ---
    out = {'SW-w': om_SW, 'SW-g': ga_SW,
           'Qiu-w': om_Qiu, 'Qiu-g': ga_Qiu,
           'Gao-w': om_gao_2010, 'Gao-g': ga_gao_2010}
    return out


def GK_linear_fitting(elong, q_safety_factor, khat):
    # FIT-w and FIT-g are normalized to sqrt(2) * vti / R_0
    def Fw(g, e, q, k):
        res = (g[0] + g[1]*k**2 + g[2]*k**4) / (1 + g[3]*k) *\
            np.exp(-g[4]*q**2) / (1 + g[5]*e)
        return res

    def Fg(h, e, q, k):
        res = (h[0] + h[1] * k**2) * np.exp(-h[2] * q**2) \
                            / (1 + h[3] * e**2) +\
              (h[4] + h[5] * k**2) * np.exp(-h[6] * q**2)\
                            / (1 + h[7] * e**4)
        return res

    Cw = np.array([3.7733, 6.3505, -1.9741e1, 1.3557e-1, 1.4620e-3, 1.1684])
    Cg = np.array([-1.2494e-2, -8.9688e-1, 4.5498e-2,
                   -1.9884e-1, -1.1248e-2, -2.5481, - 5.3340e-3, 7.7748e-1])
    out = {
        'FIT-w': Fw(Cw, elong, q_safety_factor, khat),
        'FIT-g': Fg(Cg, elong, q_safety_factor, khat)
    }
    return out


def get_gk_fit(dd, oo):
    # curves, where to add the data
    curves = oo.get('curves', crv.Curves())
    sel_norm = oo.get('sel_norm', 'wci')  # -> 'wci', 'khz', 'csa', 'csr'
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'
    sel_res = oo.get('sel_res', 'w')  # -> 'w' (frequency), 'g' (damping rate)
    col = oo.get('col', 'blue')

    # radial grid:
    r = oo.get('r', dd['pf'].nT_equil['s'])

    # assumption on the GAM wave-number
    kr = oo.get('kr', 2 * np.pi / (0.1 * dd['a0']))
    line_k = 'k = {:0.3f}'.format(kr * dd['rhoL_speak'])

    # equilibrium profiles
    rd.q(dd)

    if sel_r == 's':
        qr = dd['q']['s']
        Tr = dd['pf'].nT_equil['s']
    if sel_r == 'psi':
        qr = dd['q']['s']**2
        Tr = dd['pf'].nT_equil['psi']
    qs  = np.interp(r, qr, dd['q']['data'])
    T_J = np.interp(r, Tr, dd['pf'].nT_equil['T_J'])

    # possible plasma elongation:
    rd.elongation(dd)
    elong = oo.get('elong', dd['elong'])

    # first species parameters
    mf = dd['pf'].mass
    Zf = dd['pf'].Z
    tau_i = dd['pf'].tau

    # radial profile of the GAM frequency
    res = np.zeros(np.size(r))
    for id_s1 in range(np.size(r)):
        q1 = qs[id_s1]
        rhoL = ymath.find_rhoL(T_J[id_s1], dd['B0'], mf, Zf)
        khat = kr * rhoL

        # th = omegagamma_GAM(q1, 1/tau_i, khat, elong)
        fit = GK_linear_fitting(elong, q1, khat)

        vti = np.sqrt(T_J[id_s1]/mf)
        norm_vti = np.sqrt(2) * vti / dd['R0']

        res[id_s1] = fit['FIT-' + sel_res] * norm_vti

    # chose normalization:
    if sel_norm == 'khz':
        coef_norm = 1 / 1.e3
        if sel_res is 'w':
            coef_norm = coef_norm / (2*np.pi)
    if sel_norm == 'wci':
        coef_norm = 1. / dd['wc']
    if sel_norm == 'csa':
        coef_norm = 1. / (dd['cs']/dd['a0'])
    if sel_norm == 'csr':
        coef_norm = 1. / (dd['cs']/dd['R0'])
    res = res * coef_norm

    # result curve
    curves.new('aug20787').XS(r).YS(res) \
        .leg('LIN.\ GK\ fit.\ ' + line_k) \
        .sty(':').col(col)
    return curves


def get_gao(dd, oo):
    # curves, where to add the data
    curves = oo.get('curves', crv.Curves())
    sel_norm = oo.get('sel_norm', 'wci')  # -> 'wci', 'khz', 'csa', 'csr'
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'
    sel_res = oo.get('sel_res', 'w')  # -> 'w' (frequency), 'g' (damping rate)
    col = oo.get('col', 'blue')

    # radial grid:
    r = oo.get('r', dd['pf'].nT_equil['s'])

    # assumption on the GAM wave-number
    kr = oo.get('kr', 2 * np.pi / (0.1 * dd['a0']))
    line_k = 'k = {:0.3f}'.format(kr * dd['rhoL_speak'])

    # equilibrium profiles
    rd.q(dd)

    if sel_r == 's':
        qr = dd['q']['s']
        Tr = dd['pf'].nT_equil['s']
    if sel_r == 'psi':
        qr = dd['q']['s']**2
        Tr = dd['pf'].nT_equil['psi']
    qs  = np.interp(r, qr, dd['q']['data'])
    T_J = np.interp(r, Tr, dd['pf'].nT_equil['T_J'])

    # possible plasma elongation:
    rd.elongation(dd)
    elong = oo.get('elong', dd['elong'])

    # first species parameters
    mf = dd['pf'].mass
    Zf = dd['pf'].Z
    tau_i = dd['pf'].tau

    # radial profile of the GAM frequency
    res = np.zeros(np.size(r))
    for id_s1 in range(np.size(r)):
        q1 = qs[id_s1]
        rhoL = ymath.find_rhoL(T_J[id_s1], dd['B0'], mf, Zf)
        khat = kr * rhoL

        th = omegagamma_GAM(q1, 1/tau_i, khat, elong)

        vti = np.sqrt(T_J[id_s1]/mf)
        norm_vti = np.sqrt(2) * vti / dd['R0']

        res[id_s1] = th['Gao-' + sel_res] * norm_vti

    # chose normalization:
    if sel_norm == 'khz':
        coef_norm = 1 / 1.e3
        if sel_res is 'w':
            coef_norm = coef_norm / (2*np.pi)
    if sel_norm == 'wci':
        coef_norm = 1. / dd['wc']
    if sel_norm == 'csa':
        coef_norm = 1. / (dd['cs']/dd['a0'])
    if sel_norm == 'csr':
        coef_norm = 1. / (dd['cs']/dd['R0'])
    res = res * coef_norm

    # result curve
    curves.new('aug20787').XS(r).YS(res) \
        .leg('Gao\ 2010\ ' + line_k) \
        .sty('--').col(col)
    return curves
