import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import write_data as wr
import general as gn
import zf_gam as zf
import ITG_gamma as itg
import transport
import equil_profiles
import Distribution as distribution
import MPR
import gam_theory
import gam_exp
import Geom as geom
import numpy as np
import scipy.signal
from scipy import interpolate
import h5py as h5
from scipy.signal import correlate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(wr)
    mix.reload_module(gn)
    mix.reload_module(zf)
    mix.reload_module(itg)
    mix.reload_module(transport)
    mix.reload_module(equil_profiles)
    mix.reload_module(distribution)
    mix.reload_module(MPR)
    mix.reload_module(gam_theory)
    mix.reload_module(gam_exp)
    mix.reload_module(geom)


# find correlation between two signals:
def correlation_two_vars(dd, oo):
    def aver_var(vvar, s, s_interval, tit_var, flag_aver):
        ids_s, _, lines_s = \
            mix.get_interval(s, [s_interval], 's', '0.3f')
        tit_var += ':\ ' + lines_s[0]

        ids_s = ids_s[0]
        var_avr = vvar[:, ids_s[0]:ids_s[-1]+1]

        if flag_aver == 'mean':
            var_avr = np.mean(var_avr, axis=1)
        if flag_aver == 'rms':
            var_avr = np.sqrt(np.mean(var_avr ** 2, axis=1))
        return var_avr, tit_var

    # --- var1 ---
    oo_var = dict(oo)
    oo_var.update({
        'opt_var':      oo['opt_var1'],
        'opt_type':     oo['opt_type1'],
        'species_name': oo.get('species_name1', 'deuterium')
    })
    var1 = choose_signal_for_comparison(dd, oo_var)

    # --- var2 ---
    oo_var = dict(oo)
    oo_var.update({
        'opt_var':      oo['opt_var2'],
        'opt_type':     oo['opt_type2'],
        'species_name': oo.get('species_name2', 'deuterium')
    })
    var2 = choose_signal_for_comparison(dd, oo_var)

    # --- take the signals in some radial intervals:
    s_interval1 = oo.get('s_interval1', None)
    flag_aver1 = oo.get('flag_aver1', 'mean')
    s_interval2 = oo.get('s_interval2', None)
    flag_aver2 = oo.get('flag_aver2', 'mean')

    var1_avr, tit1_var = \
        aver_var(var1['var'], var1['s'], s_interval1, var1['tit'], flag_aver1)
    var2_avr, tit2_var = \
        aver_var(var2['var'], var2['s'], s_interval2, var2['tit'], flag_aver2)

    # --- calculate a time delay ---
    oo_dt = dict(oo)
    oo_dt.update({
        'var1':    var1_avr,  'var2':    var2_avr,
        'grid_t1': var1['t'], 'grid_t2': var2['t'],
        'vars_names': [tit1_var, tit2_var]
    })
    gn.find_time_delay(dd, oo_dt)
    n_vars = oo.get('nvars', 1)
    vvars, ts, tit_vars, lines_s = [], [], [], []
    for i_var in range(n_vars):
        line_id_var = '{:d}'.format(i_var + 1)

        # choose a signal
        opt_type = oo.get('opt_type' + line_id_var, None)
        opt_var = oo.get('opt_var' + line_id_var, None)
        species_name = oo.get('species_name' + line_id_var, None)
        oo_var = {
            'opt_type': opt_type,
            'opt_var': opt_var,
            'species_name': species_name,
        }
        res = choose_signal_for_comparison(dd, oo_var)
        data, s, t, tit = res['var'], res['s'], res['t'], res['tit']

        # average the signal or take it at particular radial points:
        opt_av = oo.get('opt_av' + line_id_var, 's1')
        line_s = ''
        if opt_av == 'rms' or opt_av == 'mean':
            oo_avs = dict(oo)
            s_interval = oo.get('s_interval' + line_id_var, None)
            oo_avs.update({
                'vars': [data], 'ts': [t], 'ss': [s],
                'opts_av': [opt_av],
                'tit': tit, 'vars_names': [tit],
                'ns_av': 1,
                's_av1_start': s_interval[0], 's_av1_end': s_interval[-1]
            })
            res = mix.find_avs(oo_avs)[0]
            data, line_s = res['data'][0], res['lines_avs'][0]
        if opt_av == 's1':
            s1 = oo.get('s_point', None)
            ids_s, s1 = mix.find(s, s1)
            data = data[:, ids_s[0]:ids_s[-1] + 1]
            line_s = 's = {:0.3f}'.format(s1)

        # save the signal
        vvars.append(data)
        ts.append(t)
        tit_vars.append(tit)
        lines_s.append(line_s)


# choose variables:
def choose_signals(dd, oo):
    n_vars = oo.get('nvars', 1)
    vvars, ts, ss, tit_vars, lines_avr = [], [], [], [], []
    for i_var in range(n_vars):
        line_id_var = '{:d}'.format(i_var + 1)

        # choose a signal
        opt_type     = oo.get('opt_type'     + line_id_var, None)
        opt_var      = oo.get('opt_var'      + line_id_var, None)
        species_name = oo.get('species_name' + line_id_var, None)
        oo_var = {
            'opt_type':     opt_type,
            'opt_var':      opt_var,
            'species_name': species_name,
        }
        res = choose_signal_for_comparison(dd, oo_var)
        data, s, t, tit = res['var'], res['s'], res['t'], res['tit']

        # average the signal or take it at particular radial points:
        opt_av = oo.get('opt_av' + line_id_var, '')
        line_avr = ''
        if opt_av is None:
            line_avr = ''
        if opt_av == 'rms' or opt_av == 'mean':
            oo_avs = dict(oo)
            s_interval = oo.get('s_interval' + line_id_var, None)
            oo_avs.update({
                'vars': [data], 'ts': [t], 'ss': [s],
                'opts_av': [opt_av],
                'tit': tit, 'vars_names': [tit],
                'ns_av': 1,
                's_av1_start': s_interval[0], 's_av1_end': s_interval[-1]
            })
            res = mix.find_avs(oo_avs)[0]
            data, line_avr = res['data'][0], res['lines_avs'][0]
        if opt_av == 's1':
            s1 = oo.get('s_point', None)
            ids_s, s1 = mix.find(s, s1)
            data = data[:, ids_s[0]:ids_s[-1] + 1]
            line_avr = 's = {:0.3f}'.format(s1)

        # save the signal
        vvars.append(data)
        ts.append(t)
        ss.append(s)
        tit_vars.append(tit)
        lines_avr.append(line_avr)

    # results:
    res = {
        'vars': vvars, 'ts': ts, 'ss': ss,
        'tit_vars': tit_vars, 'lines_avr': lines_avr
    }
    return res


# choose a signal for comparison:
def choose_signal_for_comparison(dd, oo):
    opt_type = oo.get('opt_type', 'zonal')
    oo_var = dict(oo)

    out = {}
    if opt_type == 'zonal':
        out = zf.choose_var(dd, oo_var)
    if opt_type == 'transport':
        out = transport.choose_var(dd, oo_var)
    if opt_type == 'nonzonal':
        out = itg.choose_var(dd, oo_var)
    if opt_type == 'nonzonal-s':
        out = itg.choose_var_s(dd, oo_var)
    if opt_type == 'equ-profile':
        out = equil_profiles.choose_var(dd, oo_var)

    return out


# NEW: 2d plot: (x,y)
def plot_vars_2d(oo):
    # oo.ovars = [[type, opt_var, opts], [], ...]
    # oo.dds = [dd1, dd2, dd3, ...]
    # oo.tit_plot
    # oo.labx, oo.laby
    oo_use = dict(oo)

    n_vars = len(oo['ovars'])

    # correct averaging parameter
    avrs_use = list(oo['avrs'])
    for ivar in range(n_vars):
        if len(avrs_use[ivar]) >= 2:
            avrs_use[ivar][1] = 'none-'
        else:
            avrs_use[ivar].append('none-')
    oo_use.update({
        'avrs': avrs_use,
    })

    vvars = choose_vars(oo_use)

    # additional data:
    labx = oo.get('labx', 't[\omega_{ci}^{-1}]')  # x-label
    laby = oo.get('labx', 's')  # y-label
    tit_plot   = oo.get('tit_plot', '')  # title

    flag_norm     = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)

    # plotting:
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        leg = vvar['legs'][0]
        if tit_plot is '':
            tit_plot_res = leg
        else:
            tit_plot_res = tit_plot

        labx = vvar.get('labx', labx)
        laby = vvar.get('laby', laby)
        curves = crv.Curves().xlab(labx).ylab(laby).tit(tit_plot_res)
        curves.flag_norm     = flag_norm
        curves.flag_semilogy = flag_semilogy
        x1, ids_x1 = mix.get_array_oo(oo, vvar[vvar['x1']], vvar['x1'])
        x2, ids_x2 = mix.get_array_oo(oo, vvar[vvar['x2']], vvar['x2'])
        data = mix.get_slice(vvars[ivar]['data'], ids_x1, ids_x2)
        curves.new().XS(x1).YS(x2).ZS(data).lev(60)
        cpr.plot_curves_3d(curves)


# NEW: 1d plot: plot vars along t or s:
def plot_vars_1d(oo):
    # oo.ovars = [[type, opt_var, opts], [], ...]
    # oo.dds = [dd1, dd2, dd3, ...]
    # oo.avrs = [[coords, type-coord_av, domains], [], ...]
    # oo.tit_plot
    # oo.labx
    # oo.laby
    # oo.func_mod

    vvars = choose_vars(oo)
    n_vars = len(vvars)

    # additional data:
    labx       = oo.get('labx', '')  # x-label
    laby       = oo.get('laby', '')  # y-label
    tit_plot   = oo.get('tit_plot', '')  # title

    flag_norm     = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)

    func_mod = oo.get('func_mod', None)

    # plotting:
    curves = crv.Curves().xlab(labx).ylab(laby).tit(tit_plot)
    curves.flag_norm     = flag_norm
    curves.flag_semilogy = flag_semilogy
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        datas = vvar['data']
        x_init, ids_x = mix.get_array_oo(oo, vvar['x'], 'x')
        legs  = vvar['legs']
        ns_av = len(datas)
        for is_av in range(ns_av):
            data = mix.get_slice(datas[is_av], ids_x)

            x = x_init
            if func_mod is not None:
                x, data = func_mod(x_init, data)

            curves.new().XS(x).YS(data).leg(legs[is_av])
    # curves.set_colors_styles()
    if len(curves.list_curves) is not 0:
        cpr.plot_curves(curves)


# NEW: FFT 1d:
def plot_fft_1d(oo):
    # oo.ovars = [[type, opt_var, opts], [], ...]
    # oo.dds = [dd1, dd2, dd3, ...]
    # oo.avrs = [[coords, type-coord_av, domains], [], ...]
    # oo.x_fft_domains = [domain1, domain2, ...]
    # oo.w_domain = [w_start, w_end],
    # oo.tit_plot
    # oo.labx
    # oo.laby
    # oo.sel_norm - normalization of frequency: 'khz', 'wc', 'csa', 'csr'
    # oo.flag_f2 - one or two-sided FFT

    # Domains (oo.x_fft_domains), where FFT will be found, are the same for all
    # variables.

    # choose variables and their averaging:
    vvars = get_fft_1d(oo)
    n_vars = len(vvars['datas'])

    dds = oo.get('dds', None)

    # frequency domain for plotting:
    w_domain = oo.get('w_domain', None)

    # - frequency normalization (notation) -
    sel_norm = oo.get('sel_norm', 'wc')

    line_w = ''
    if sel_norm == 'khz':
        line_w = '\omega,\ kHz'
    if sel_norm == 'wc':
        line_w = '\omega[\omega_c]'
    if sel_norm == 'csa':
        line_w = '\omega[c_s/a_0]'
    if sel_norm == 'csr':
        line_w = '\omega[c_s/R_0]'

    # additional data:
    labx = oo.get('labx', line_w)  # x-label
    laby = oo.get('laby', '')  # y-label
    tit_plot = oo.get('tit_plot', '')  # title

    flag_norm = oo.get('flag_norm', False)

    curves = crv.Curves().xlab(labx).ylab(laby).tit(tit_plot)
    curves.flag_norm = flag_norm
    for ivar in range(n_vars):
        dd = dds[ivar]
        data = vvars['datas'][ivar]
        w    = vvars['ws'][ivar]
        leg  = vvars['legs'][ivar]

        # - frequency normalization (coefficient) -
        coef_norm = None
        if sel_norm == 'khz':
            coef_norm = dd['wc'] / 1.e3
        if sel_norm == 'wc':
            coef_norm = 2 * np.pi
        if sel_norm == 'csa':
            coef_norm = 2 * np.pi * dd['wc'] / (dd['cs'] / dd['a0'])
        if sel_norm == 'csr':
            coef_norm = 2 * np.pi * dd['wc'] / (dd['cs'] / dd['R0'])
        w = w * coef_norm

        # - frequency interval for plotting -
        ids_w, w, _ = mix.get_ids(w, w_domain)
        data = mix.get_slice(data, ids_w)

        # - plotting -
        curves.new() \
            .XS(w) \
            .YS(data) \
            .leg(leg)

    curves.set_colors_styles()
    cpr.plot_curves(curves)


# NEW: FFT 2d:
def plot_fft_2d(oo):
    # oo.ovars = [[type, opt_var, opts], [], ...]
    # oo.dds = [dd1, dd2, dd3, ...]
    # oo.avrs = [[coords, type-coord_av, domains], [], ...]
    # oo.x_fft_domains = [domain1, domain2, ...]
    # oo.w_domain = [w_start, w_end],
    # oo.dep_domain = [dep_start, dep_end]
    # oo.tit_plot
    # oo.labx
    # oo.laby
    # oo.sel_cmp - colormap
    # oo.sel_norm - normalization of frequency: 'khz', 'wc', 'csa', 'csr'
    # oo.flag_f2 - one or two-sided FFT
    # oo.coord_fft - coordinate axis, along which FFT will be taken

    # Domains (oo.x_fft_domains), where FFT will be found, are the same for all
    # variables.

    # choose variables and their averaging:
    vvars = get_fft_2d(oo)
    n_vars = len(vvars['datas'])

    dds = oo.get('dds', None)

    # theoretical and experimental
    flag_gao = oo.get('flag_gao', True)
    flag_gk_fit = oo.get('flag_gk_fit', False)
    flag_aug20787 = oo.get('flag_aug20787', False)

    # - frequency normalization (notation) -
    sel_norm = oo.get('sel_norm', 'wc')

    line_w = ''
    if sel_norm == 'khz':
        line_w = '\omega,\ kHz'
    if sel_norm == 'wc':
        line_w = '\omega[\omega_c]'
    if sel_norm == 'csa':
        line_w = '\omega[c_s/a_0]'
    if sel_norm == 'csr':
        line_w = '\omega[c_s/R_0]'

    # additional data:
    laby = oo.get('laby', line_w)  # y-label
    sel_cmp = oo.get('sel_cmp', 'jet')

    for ivar in range(n_vars):
        dd = dds[ivar]
        data      = vvars['datas'][ivar]  # [coord_dep, w]
        w         = vvars['ws'][ivar]
        coord_dep = vvars['ys'][ivar]
        name_dep  = vvars['names_y'][ivar]
        tit       = vvars['legs'][ivar]

        labx = oo.get('labx', name_dep)  # x-label
        tit = oo.get('tit_plot', tit)  # title

        # - create a figure -
        curves = crv.Curves().xlab(labx).ylab(laby).tit(tit)

        # intervals for plotting:
        w_domain = oo.get('w_domain', [w[0], w[-1]])
        dep_domain = oo.get('dep_domain', [coord_dep[0], coord_dep[-1]])

        # - frequency normalization (coefficient) -
        coef_norm = None
        if sel_norm == 'khz':
            coef_norm = dd['wc'] / 1.e3
        if sel_norm == 'wc':
            coef_norm = 2 * np.pi
        if sel_norm == 'csa':
            coef_norm = 2 * np.pi * dd['wc'] / (dd['cs'] / dd['a0'])
        if sel_norm == 'csr':
            coef_norm = 2 * np.pi * dd['wc'] / (dd['cs'] / dd['R0'])
        w = w * coef_norm

        # - frequency and dependence-coordinate interval for plotting -
        ids_w,           w, _ = mix.get_ids(w,         w_domain)
        ids_dep, coord_dep, _ = mix.get_ids(coord_dep, dep_domain)
        data = mix.get_slice(data, ids_dep, ids_w)

        # - limits, legend position -
        curves\
            .xlim([coord_dep[0], coord_dep[-1]])\
            .ylim([w[0], w[-1]])\
            .leg_pos('upper left')

        # - numerical results -
        curves.new().XS(coord_dep).YS(w).ZS(data)\
            .lev(60).cmp(sel_cmp)

        # theoretical and experimental curves:
        oo_th = {'curves': curves, 'sel_norm': sel_norm,
                 'r': coord_dep, 'col': 'white'}

        if flag_aug20787:
            curves = gam_exp.exp_AUG20787(dd, oo_th)
        if flag_gao:
            curves = gam_theory.get_gao(dd, oo_th)
        if flag_gk_fit:
            curves = gam_theory.get_gk_fit(dd, oo_th)

        # - build 2d figure -
        cpr.plot_curves_3d(curves)


# NEW: find dynamic rate for several signals and in several time intervals
def find_gamma(oo):
    # the same time domains for different signals:
    # oo.t_work_domains - where to find gamma
    # oo.ovars, oo.avrs

    vvars = choose_vars(oo)
    nvars = len(vvars)

    tit_plot = oo.get('tit_plot', '')  # title
    labx = oo.get('labx', 't[wci^{-1}]')
    laby = oo.get('laby', '')  # y-label
    flag_semilogy = oo.get('flag_semilogy', True)
    line_g = '\gamma[\omega_{ci}]'

    t_work_domains = oo.get('t_work_domains', None)
    nt_domains = len(t_work_domains)

    # --- find dynamic rates ---
    vvars_work_vars, ts_work_vars, lines_t_work_vars, gs_est_vars = [], [], [], []
    for id_var in range(nvars):
        # current data
        data = vvars[id_var]['data'][0]
        t    = vvars[id_var]['x']

        vvars_work, ts_work, lines_t_work, gs_est = [], [], [], []
        for id_t_interval in range(nt_domains):
            t_work_domain = t_work_domains[id_t_interval]

            # time domain where the gamma will be computed
            ids_t_work, t_work, line_t_work = mix.get_ids(t, t_work_domain)
            data_work = np.squeeze(data[ids_t_work])

            # find gamma in this time domain
            g_est = ymath.estimate_g(t_work, data_work)

            # results for a given time interval
            vvars_work.append(data_work)
            ts_work.append(t_work)
            lines_t_work.append(line_t_work)
            gs_est.append(g_est)

        # results for a given variable
        vvars_work_vars.append(vvars_work)
        ts_work_vars.append(ts_work)
        lines_t_work_vars.append(lines_t_work)
        gs_est_vars.append(gs_est)

    # --- PLOT INITIAL SIGNALS and WORKING DOMAINS (curves_work) ---
    styles_loc = ['-.'] * 20

    curves_work = crv.Curves().xlab(labx).ylab(laby).tit(tit_plot)
    curves_work.flag_semilogy = flag_semilogy

    curves_fit = crv.Curves().xlab(labx).tit(tit_plot)
    curves_fit.flag_semilogy = flag_semilogy

    for id_var in range(nvars):
        data = vvars[id_var]['data'][0]
        t    = vvars[id_var]['x']
        leg = vvars[id_var]['legs'][0]

        vvars_work   = vvars_work_vars[id_var]
        ts_work      = ts_work_vars[id_var]
        lines_t_work = lines_t_work_vars[id_var]
        gs_est       = gs_est_vars[id_var]

        # intervals
        t, ids_t = mix.get_array_oo(oo, t, 't')
        data = mix.get_slice(data, ids_t)

        # plot initial signals
        curves_work.new().XS(t).YS(data).leg(leg)
        curves_fit.new().XS(t).YS(data).leg(leg)

        # plot data for every time interval
        for id_t_interval in range(len(ts_work)):
            data_work   = vvars_work[id_t_interval]
            t_work      = ts_work[id_t_interval]
            line_t_work = lines_t_work[id_t_interval]
            g_est       = gs_est[id_t_interval]

            line_legend_fit = leg + ':\ ' + line_t_work

            # work domain
            curves_work.new() \
                .XS(t_work) \
                .YS(data_work) \
                .sty(styles_loc[id_t_interval])

            # fitting
            # curves_fit.new() \
            #     .XS(g_est['x_peaks']) \
            #     .YS(g_est['y_peaks']) \
            #     .leg(line_legend_fit).sty('o').col('green')
            curves_fit.new() \
                .XS(g_est['x_fit']) \
                .YS(g_est['y_fit']) \
                .leg(line_legend_fit).sty(styles_loc[id_t_interval])

            print(leg + ': ' + line_t_work + ': ' + line_g + ' = {:0.3e}'.format(g_est['g']))

    cpr.plot_curves(curves_work)
    cpr.plot_curves(curves_fit)

    return


# NEW: choose variables:
def choose_vars(oo):
    # oo.ovars = [[type1 var1 species1 ...], ...]
    # oo.dds = [dd1, dd2, ...]
    # oo.var_legs1, oo.var_legs2, ...
    # oo.sel_legs1, oo.sel_legs2, ...   # 'full', 'woPr', 'woVN', 'woPrVN', 'woAVR', 'Pr'
    # ------------------
    # return vvars
    # vvars[i].[data x lines_avr legs opt_av tit labx laby]

    ovars = oo.get('ovars', [])
    dds   = oo.get('dds', [])
    avrs  = oo.get('avrs', [])

    n_vars = len(ovars)

    vvars, xs, tit_vars, lines_avr = [], [], [], []
    for i_var in range(n_vars):
        dd = dds[i_var]
        line_id_var = '{:d}'.format(i_var + 1)
        vvar = average_one_var(avrs[i_var], ovars[i_var], dds[i_var])
        var_legs = oo.get('var_legs' + line_id_var, None)
        sel_legs = oo.get('sel_legs' + line_id_var, 'full')
        pr_name = dd['project_name']
        if pr_name is not '':
            pr_name += ':\ '

        nlegs = len(vvar['lines_avr'])
        legs = []
        if var_legs is not None:
            legs = var_legs
        elif sel_legs == 'full':
            for id_leg in range(nlegs):
                legs.append(
                    pr_name +
                    vvar['lines_avr'][id_leg] + ':\ ' +
                    vvar['tit']
                )
        elif sel_legs == 'woPr':
            for id_leg in range(nlegs):
                legs.append(
                    vvar['lines_avr'][id_leg] + ':\ ' +\
                    vvar['tit']
                )
        elif sel_legs == 'woVN':
            for id_leg in range(nlegs):
                legs.append(
                    pr_name + vvar['lines_avr'][id_leg]
                )
        elif sel_legs == 'woPrVN':
            for id_leg in range(nlegs):
                legs.append(
                    vvar['lines_avr'][id_leg]
                )
        elif sel_legs == 'woAVR':
            for id_leg in range(nlegs):
                legs.append(
                    pr_name + vvar['tit']
                )
        elif sel_legs == 'Pr':
            for id_leg in range(nlegs):
                legs.append(
                    pr_name
                )
        vvar['legs'] = legs

        # save the signal
        vvars.append(vvar)
    return vvars


# NEW: averaging of a variable:
def average_one_var(avr, ovar, dd):
    # oo.avrs = [[sel_coords type opt1 opt2 ...], [], ...]

    sel_coords  = avr[0]  # chosen coordinate system
    oavr = avr[1:len(avr)]

    vvar_res = {}

    # choose coordinate system, where the system will be considered:
    if sel_coords == 'ts':
        vvar_res = choose_one_var_ts(ovar, dd)
    if sel_coords == 'tchi':
        vvar_res = choose_one_var_tchi(ovar, dd)
    if sel_coords == 'tvpar':
        vvar_res = choose_one_var_tvpar(ovar, dd)
    if sel_coords == 'vparmu':
        vvar_res = choose_one_var_vparmu(ovar, dd)
    if sel_coords == 't':
        vvar_res = choose_one_var_t(ovar, dd)
        oavr = ['none-']
        vvar_res['laby'] = ''

    # averaging of the chosen system
    vvar_avr = ymath.avr_x1x2(vvar_res, oavr)

    vvar_avr['tit']  = vvar_res['tit']
    vvar_avr['labx'] = vvar_res['labx']
    vvar_avr['laby'] = vvar_res['laby']

    return vvar_avr


# NEW: take variable(t,s)
def choose_one_var_ts(ovar, dd):
    opt_type = ovar[0]
    oovar = ovar[1:len(ovar)]

    vvar = {}
    if opt_type == 'zonal':
        vvar = zf.choose_one_var_ts(oovar, dd)
    if opt_type == 'transport':
        vvar = transport.choose_one_var_ts(oovar, dd)
    if opt_type == 'nonzonal':
        vvar = itg.choose_one_var_ts(oovar, dd)
    if opt_type == 'equ-profile':
        vvar = equil_profiles.choose_one_var_ts(oovar, dd)

    # save information about coordinate system:
    vvar['x1'], vvar['fx1'], vvar['labx'] = 't', '{:0.3e}', 't'
    vvar['x2'], vvar['fx2'], vvar['laby'] = 's', '{:0.3f}', 's'

    return vvar


# NEW: take variable(t,chi)
def choose_one_var_tchi(ovar, dd):
    opt_type = ovar[0]
    oovar = ovar[1:len(ovar)]

    vvar = {}
    if opt_type == 'nonzonal':
        vvar = itg.choose_one_var_tchi(oovar, dd)

    # save information about coordinate system:
    vvar['x1'], vvar['fx1'], vvar['labx'] = 't', '{:0.3e}', 't'
    vvar['x2'], vvar['fx2'], vvar['laby'] = 'chi', '{:0.3f}', '\chi'

    return vvar


# NEW: take variable(t,vpar)
def choose_one_var_tvpar(ovar, dd):
    opt_type = ovar[0]
    oovar = ovar[1:len(ovar)]

    vvar = {}
    if opt_type == 'distribution':
        vvar = distribution.choose_one_var_tvpar(oovar, dd)
    # if opt_type == 'mpr':
    #     vvar = mpr.choose_one_var_tvpar(oovar, dd)

    # save information about coordinate system:
    vvar['x1'], vvar['fx1'], vvar['labx'] = 't',    '{:0.3e}', 't'
    vvar['x2'], vvar['fx2'], vvar['laby'] = 'vpar', '{:0.3f}', 'v_{\parallel}'

    return vvar


# NEW: take variable(t,vpar)
def choose_one_var_vparmu(ovar, dd):
    opt_type = ovar[0]
    oovar = ovar[1:len(ovar)]

    vvar = {}
    if opt_type == 'mpr':
        vvar = MPR.choose_one_var_vparmu(oovar, dd)

    # save information about coordinate system:
    vvar['x1'], vvar['fx1'], vvar['labx'] = 'mu',   '{:0.3f}', '\mu'
    vvar['x2'], vvar['fx2'], vvar['laby'] = 'vpar', '{:0.3f}', 'v_{\parallel}'

    return vvar


# NEW: take variable(t,s)
def choose_one_var_t(ovar, dd):
    opt_type = ovar[0]
    oovar = ovar[1:len(ovar)]

    vvar = {}
    if opt_type == 'mpr':
        vvar = MPR.choose_one_var_t(oovar, dd)

    # save information about coordinate system:
    vvar['fx'], vvar['labx'] = '{:0.3e}', 't'

    return vvar


# NEW: get FFT 1d
def get_fft_1d(oo):
    # oo.x_fft_domains = [domain1, domain2, ...]
    # oo.w_domain = [w_start, w_end],
    # oo.flag_f2 - one or two sided FFT

    # Domains (oo.x_fft_domains), where FFT will be found, are the same for all
    # variables.

    # choose variables and their averaging:
    vvars = choose_vars(oo)
    n_vars = len(vvars)

    # domains, where the FFT will be calculated:
    x_fft_domains = oo.get('x_fft_domains', None)
    ndomains_fft = len(x_fft_domains)

    # one or two-sided FFT
    flag_f2 = oo.get('flag_f2', False)

    vars_fft, ws, legs_fft = [], [], []
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        datas = vvar['data']
        legs = vvar['legs']
        x        = vvar['x']
        name_x   = vvar['name_x']
        format_x = vvar['format_x']

        # calculate domains, where FFT will be calculated:
        xs_fft, ids_xs_fft, lines_x_fft = [], [], []
        for idomain in range(ndomains_fft):
            ids_x_fft, x_fft, line_x_fft = mix.get_ids(
                x, x_fft_domains[idomain], format_x
            )
            xs_fft.append(x_fft)
            ids_xs_fft.append(ids_x_fft)
            lines_x_fft.append(line_x_fft)

        ns_av = len(datas)
        for is_av in range(ns_av):
            data_av = datas[is_av]
            leg = legs[is_av]
            for idomain in range(ndomains_fft):
                x_fft = xs_fft[idomain]
                ids_x_fft = ids_xs_fft[idomain]
                line_x_fft = lines_x_fft[idomain]

                data = mix.get_slice(data_av, ids_x_fft)

                # - frequency grid -
                ffres = ymath.fft_y(x_fft)
                if flag_f2:
                    w = ffres['w2']
                else:
                    w = ffres['w']

                # - FFT -
                var_fft = ymath.fft_y(x_fft, data, {'flag_f2_arranged': flag_f2})

                # - result FFT -
                if flag_f2:
                    var_fft = var_fft['f2_arranged']
                else:
                    var_fft = var_fft['f']

                leg_line = leg + ':\ ' + name_x + ' = ' + line_x_fft

                # save
                vars_fft.append(var_fft)
                ws.append(w)
                legs_fft.append(leg_line)

    res = {
        'datas': vars_fft,
        'ws': ws,
        'legs': legs_fft
    }
    return res


# NEW: get FFT 2d
def get_fft_2d(oo):
    # oo.x_fft_domains = [domain1, domain2, ...]
    # oo.w_domain = [w_start, w_end],
    # oo.flag_f2 - one or two sided FFT
    # oo.coord_fft - coordinate axis, along which FFT will be taken

    # Domains (oo.x_fft_domains), where FFT will be found, are the same for all
    # variables.

    # choose variables and their averaging:
    vvars = choose_vars(oo)
    n_vars = len(vvars)

    # coordinate axis, along which FFT will be taken
    coord_fft = oo.get('coord_fft', None)

    # domains, where the FFT will be calculated:
    x_fft_domains = oo.get('x_fft_domains', None)
    ndomains_fft = len(x_fft_domains)

    # one or two-sided FFT
    flag_f2 = oo.get('flag_f2', False)

    vars_fft, ws, ys, legs_fft, names_y = [], [], [], [], []
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        data = vvar['data']
        leg = vvar['legs'][0]
        x    = vvar[coord_fft]
        if coord_fft == vvar['x1']:
            name_x = vvar['x1']
            format_x = vvar['fx1']
            y = vvar[vvar['x2']]
            name_y = vvar['x2']
        else:
            name_x = vvar['x2']
            format_x = vvar['fx2']
            y = vvar[vvar['x1']]
            name_y = vvar['x1']

        # calculate domains, where FFT will be calculated:
        xs_fft, ids_xs_fft, lines_x_fft = [], [], []
        for idomain in range(ndomains_fft):
            ids_x_fft, x_fft, line_x_fft = mix.get_ids(
                x, x_fft_domains[idomain], format_x
            )
            xs_fft.append(x_fft)
            ids_xs_fft.append(ids_x_fft)
            lines_x_fft.append(line_x_fft)

        for idomain in range(ndomains_fft):
            x_fft = xs_fft[idomain]
            ids_x_fft  = ids_xs_fft[idomain]
            line_x_fft = lines_x_fft[idomain]

            if coord_fft == vvar['x1']:
                data = np.squeeze(data[ids_x_fft, :])
                data = data.T
            else:
                data = np.squeeze(data[:, ids_x_fft])

            # - frequency grid -
            ffres = ymath.fft_y(x_fft)
            if flag_f2:
                w = ffres['w2']
            else:
                w = ffres['w']

            # - FFT -
            var_fft = ymath.fft_y(
                x_fft, data, {'flag_f2_arranged': flag_f2, 'axis': 1}
            )

            # - result FFT -
            if flag_f2:
                var_fft = var_fft['f2_arranged']
            else:
                var_fft = var_fft['f']

            leg_line = leg + ':\ ' + name_x + ' = ' + line_x_fft

            # save
            vars_fft.append(var_fft)
            ws.append(w)
            ys.append(y)
            legs_fft.append(leg_line)
            names_y.append(name_y)

    res = {
        'datas': vars_fft,  # datas_i[y, w]
        'ws': ws,
        'ys': ys,
        'names_y': names_y,
        'legs': legs_fft
    }
    return res


# NEW: contribution of different Fourier components:
def contribution_fft_1d(oo):
    # oo.w_area1 = [w_left, w_right]
    # oo.w_area2 = [w_left, w_right]
    # find a ratio fft[oo.w_area2]/fft[oo.w_area1]

    vvars = get_fft_1d(oo)
    dd = oo.get('dds', None)[0]

    # frequency domain for plotting:
    w_domain = oo.get('w_domain', None)

    # - frequency normalization -
    sel_norm = oo.get('sel_norm', 'wc')

    # additional data:
    laby = oo.get('laby', '')  # y-label
    tit_plot = oo.get('tit_plot', '')  # title

    flag_norm = oo.get('flag_norm', False)  # plotting of normalized or non-normalized data

    # - data -
    data = vvars['datas'][0]
    w = vvars['ws'][0]
    leg = vvars['legs'][0]

    # - frequency normalization -
    coef_norm, line_w = None, ''
    if sel_norm == 'khz':
        coef_norm = dd['wc'] / 1.e3
        line_w = '\omega,\ kHz'
    if sel_norm == 'wc':
        coef_norm = 2 * np.pi
        line_w = '\omega[\omega_c]'
    if sel_norm == 'csa':
        coef_norm = 2 * np.pi * dd['wc'] / (dd['cs'] / dd['a0'])
        line_w = '\omega[c_s/a_0]'
    if sel_norm == 'csr':
        coef_norm = 2 * np.pi * dd['wc'] / (dd['cs'] / dd['R0'])
        line_w = '\omega[c_s/R_0]'
    labx = oo.get('labx', line_w)  # x-label
    w = w * coef_norm

    # - two areas to compare -
    w_area1 = np.array(oo.get('w_area1', [w[0], w[-1]]))
    w_area2 = np.array(oo.get('w_area2', [w[0], w[-1]]))

    # - frequency interval for plotting -
    ids_w, w_plot, _ = mix.get_ids(w, w_domain)
    data_plot = mix.get_slice(data, ids_w)

    # - create shaded areas -
    area1 = geom.Fill()
    area1.xs = [w_area1[0], w_area1[-1], w_area1[-1], w_area1[0]]
    area1.ys = ['limb', 'limb', 'limu', 'limu']
    area1.color = 'grey'

    area2 = geom.Fill()
    area2.xs = [w_area2[0], w_area2[-1], w_area2[-1], w_area2[0]]
    area2.ys = ['limb', 'limb', 'limu', 'limu']
    area2.color = 'green'

    # - plotting -
    curves = crv.Curves().xlab(labx).ylab(laby).tit(tit_plot)
    curves.flag_norm = flag_norm
    curves.new() \
        .XS(w_plot) \
        .YS(data_plot) \
        .leg(leg)
    curves.newg(area1)
    curves.newg(area2)
    cpr.plot_curves(curves)

    # find ratio:
    id_wa1, _, line_wa1 = mix.get_ids(w, w_area1)
    id_wa2, _, line_wa2 = mix.get_ids(w, w_area2)

    fft_a1 = np.squeeze(np.mean(data[id_wa1]))
    fft_a2 = np.squeeze(np.mean(data[id_wa2]))

    # print result:
    line_res  = 'area1: ' + line_w + ' = ' + line_wa1 + '\n'
    line_res += 'area2: ' + line_w + ' = ' + line_wa2 + '\n'
    line_res += 'fft(area2)/fft(area1) = {:0.3f}'.format(fft_a2/fft_a1)
    print(line_res)


# NEW: contribution of different Fourier components:
def contribution_fft_2d(oo):
    # oo.coord_fft - coordinate axis, along which FFT will be taken
    # oo.w_area1 = [w_left, w_right]
    # oo.w_area2 = [w_left, w_right]
    # find a ratio fft[oo.w_area2]/fft[oo.w_area1](oo.coord_fft)

    # choose variables and their averaging:
    vvars = get_fft_2d(oo)
    dd = oo.get('dds', None)[0]

    # theoretical and experimental
    flag_gao      = oo.get('flag_gao', True)
    flag_gk_fit   = oo.get('flag_gk_fit', False)
    flag_aug20787 = oo.get('flag_aug20787', False)

    # - frequency normalization (notation) -
    sel_norm = oo.get('sel_norm', 'wc')

    # colormap
    sel_cmp = oo.get('sel_cmp', 'jet')

    data      = vvars['datas'][0]  # [coord_dep, w]
    w         = vvars['ws'][0]
    coord_dep = vvars['ys'][0]
    name_dep  = vvars['names_y'][0]
    tit       = vvars['legs'][0]

    labx = oo.get('labx', name_dep)  # x-label
    tit = oo.get('tit_plot', tit)  # title

    # intervals for plotting:
    w_domain = oo.get('w_domain', [w[0], w[-1]])
    dep_domain = oo.get('dep_domain', [coord_dep[0], coord_dep[-1]])

    # - frequency normalization (coefficient) -
    coef_norm, line_w = None, ''
    if sel_norm == 'khz':
        coef_norm = dd['wc'] / 1.e3
        line_w = '\omega,\ kHz'
    if sel_norm == 'wc':
        coef_norm = 2 * np.pi
        line_w = '\omega[\omega_c]'
    if sel_norm == 'csa':
        coef_norm = 2 * np.pi * dd['wc'] / (dd['cs'] / dd['a0'])
        line_w = '\omega[c_s/a_0]'
    if sel_norm == 'csr':
        coef_norm = 2 * np.pi * dd['wc'] / (dd['cs'] / dd['R0'])
        line_w = '\omega[c_s/R_0]'
    w = w * coef_norm
    laby = oo.get('laby', line_w)  # y-label

    # - two areas to compare -
    w_area1 = np.array(oo.get('w_area1', [w[0], w[-1]]))
    w_area2 = np.array(oo.get('w_area2', [w[0], w[-1]]))

    # - frequency and dependence-coordinate interval for plotting -
    ids_w, w_plot, _ = mix.get_ids(w, w_domain)
    ids_dep, coord_dep_plot, _ = mix.get_ids(coord_dep, dep_domain)
    data_plot = mix.get_slice(data, ids_dep, ids_w)

    # - create shaded areas -
    area1 = geom.Fill()
    area1.xs = ['liml', 'limr', 'limr', 'liml']
    area1.ys = [w_area1[0], w_area1[0], w_area1[-1], w_area1[-1]]
    area1.color = 'grey'
    area1.alpha = 0.3

    area2 = geom.Fill()
    area2.xs = ['liml', 'limr', 'limr', 'liml']
    area2.ys = [w_area2[0], w_area2[0], w_area2[-1], w_area2[-1]]
    area2.color = 'green'
    area2.alpha = 0.3

    # - create a figure -
    curves = crv.Curves().xlab(labx).ylab(laby).tit(tit)
    curves.xlim([coord_dep_plot[0], coord_dep_plot[-1]]) \
        .ylim([w_plot[0], w_plot[-1]]) \
        .leg_pos('upper left')

    # - numerical results -
    curves.new().XS(coord_dep_plot).YS(w_plot).ZS(data_plot) \
        .lev(60).cmp(sel_cmp)
    curves.newg(area1)
    curves.newg(area2)

    # - theoretical and experimental curves -
    oo_th = {'curves': curves, 'sel_norm': sel_norm,
             'r': coord_dep, 'col': 'red'}

    if flag_aug20787:
        curves = gam_exp.exp_AUG20787(dd, oo_th)
    if flag_gao:
        curves = gam_theory.get_gao(dd, oo_th)
    if flag_gk_fit:
        curves = gam_theory.get_gk_fit(dd, oo_th)

    # - build 2d figure -
    cpr.plot_curves_3d(curves)

    # - find ratio -
    id_wa1, _, line_wa1 = mix.get_ids(w, w_area1)
    id_wa2, _, line_wa2 = mix.get_ids(w, w_area2)

    fft_a1 = np.squeeze( np.mean(data[:, id_wa1], axis=1) )
    fft_a2 = np.squeeze( np.mean(data[:, id_wa2], axis=1) )
    ratio21 = fft_a2 / fft_a1
    ratio12 = fft_a1 / fft_a2

    # - plot result -
    line_res1 = 'area1: ' + line_w + ' = ' + line_wa1
    line_res2 = 'area2: ' + line_w + ' = ' + line_wa2

    curves = crv.Curves().xlab(labx).ylab('fft(area2)/fft(area1)')\
        .tit(line_res1)\
        .titn(line_res2)
    curves.new().XS(coord_dep).YS(ratio21)
    cpr.plot_curves(curves)

    curves = crv.Curves().xlab(labx).ylab('fft(area1)/fft(area2)') \
        .tit(line_res1) \
        .titn(line_res2)
    curves.new().XS(coord_dep).YS(ratio12)
    cpr.plot_curves(curves)


# NEW: check max on s:
def check_absmax_s(oo, t_point):
    # variable along s:
    oo_var = dict(oo)
    oo_var.update({
        'avrs': [
            ['ts', 'point-t', [t_point]]
        ]
    })
    vvar = choose_vars(oo_var)[0]

    oo_var.update({
        'avrs': [
            ['ts', 'max-s']
        ]
    })
    vvar_max = choose_vars(oo_var)[0]

    # additional data:
    labx = oo.get('labx', 's')  # x-label
    tit_plot = oo.get('tit_plot', '')  # title

    flag_norm     = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)

    # data
    s, ids_s = mix.get_array_oo(oo, vvar['x'], 'x')
    data = vvar['data'][0][ids_s[0]:ids_s[-1]+1]

    # max data
    id_t, _, _ = mix.get_ids(vvar_max['x'], t_point)
    data_max   = vvar_max['data'][0][id_t]
    s_max      = vvar_max['opt_av'][0][id_t]

    # plotting:
    curves = crv.Curves().xlab(labx).tit(tit_plot)
    curves.flag_norm = flag_norm
    curves.flag_semilogy = flag_semilogy
    curves.new()\
        .XS(s)\
        .YS(data)\
        .leg('data')
    curves.new() \
        .XS(s_max) \
        .YS(data_max) \
        .leg('max').sty('o')
    cpr.plot_curves(curves)


# Plot main 1d figures for the MPR diagnostic:
def MPR_plot_1d(dd, oo):
    # Local data
    DL = {}

    # --- Names ---
    sp_names = dd['kin_species_names']
    total_name = 'total'
    part_names = sp_names + [total_name]

    je_name = 'jdote_es'
    ef_name = 'efield'
    var_names = [je_name, ef_name]

    # --- GET JE and Efield (t) ---
    oo_vvar = {
        'avrs': [['t']],
        'dds': [dd],
    }
    for var_name in var_names:
        DL[var_name] = {}
        for part_name in part_names:
            oo_vvar.update({
                'ovars': [ ['mpr', var_name, part_name] ],
            })
            DL[var_name][part_name] = choose_vars(oo_vvar)[0]

    del oo_vvar

    # # Plot species signals:
    # vvar = DL[je_name]['deuterium']
    # curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit(vvar['tit'])
    # curves.new().XS(vvar['x']).YS(vvar['data']).leg(vvar['line_sum'])
    # cpr.plot_curves(curves)
    #
    # vvar = DL[ef_name]['deuterium']
    # curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit(vvar['tit'])
    # curves.new().XS(vvar['x']).YS(vvar['data'])
    # cpr.plot_curves(curves)

    # --- ELIMINATE CONTRIBUTION of Zero-Frequency Zonal Flows (ZFZF) ---

    # - Get GAM and ZFZF components of the zonal electric field -
    def obtain_e0e1(g, e, t, t_point_inf):
        id_t_inf, t_point_inf, _ = mix.get_ids(t, t_point_inf)
        id_t_inf = id_t_inf[0]

        A_zfzf = np.mean(e[id_t_inf:])
        e_gam = e - A_zfzf
        e_gam_cos = e_gam * np.exp(-g * t)

        ids_peaks, e_gam_peaks = scipy.signal.find_peaks(e_gam_cos, height=0)
        e_gam_peaks = e_gam_peaks['peak_heights']
        t_peaks = t[ids_peaks]

        A_gam = np.mean(e_gam_peaks)

        res = {
            'A_zfzf': A_zfzf, 'A_gam': A_gam,
            'e_gam': e_gam,
            'e_gam_cos': e_gam_cos,
            't_peaks': t_peaks, 'e_gam_peaks': e_gam_peaks,
        }

        return res

    # - Eliminate ZFZF from field energy and energy transfer signals -
    def obtain_gam(w, g, E, P, t, id_peak, eta):
        # E - field energy
        # P - J*E

        ids_peaks, _ = scipy.signal.find_peaks(E, height=0)

        E_peak1 = E[ids_peaks[id_peak]]
        E_peak2 = E[ids_peaks[id_peak+1]]
        if eta > 0:
            idd = np.array([E_peak1, E_peak2]).argmax()
        else:
            # in this case you should take at the beginning where
            # there are still alternating peaks
            idd = np.array([E_peak1, E_peak2]).argmin()
        if idd == 1:
            id_peak += 1

        t_peak = t[ids_peaks[id_peak]]
        Eref   = E[ids_peaks[id_peak]]

        e12 = Eref / (0.5 * eta**2 + eta * np.exp(g*t_peak) + 0.5 * np.exp(2*g*t_peak))
        e1 = np.sqrt(e12)
        e0 = e1 * eta

        Egam = E - 0.5 * e0 ** 2 - e0 * e1 * np.cos(w * t) * np.exp(g * t)
        Pgam = P + e0*e1 * np.exp(g*t) * (g * np.cos(w*t) - w * np.sin(w*t))

        res = {
            'ids_peaks': ids_peaks,
            'id_peak': id_peak,
            'Egam': Egam,
            'Pgam': Pgam
        }
        return res

    # - FREQUENCY and DAMPING RATE -
    w = 3.9e-3
    g = -1.1e-4

    t_point_inf = 1.0e4

    # - Field energy and time -
    Efield = DL[ef_name]['deuterium']['data']
    t = DL[ef_name]['deuterium']['x']

    # - Energy transfer signal -
    P = DL[je_name]['deuterium']['data']
    P = np.interp(t, DL[je_name]['deuterium']['x'], P)

    # - Zonal electric field at s1 -
    s1 = 0.76
    oo_erbar = {
        'ovars': [ ['zonal', 'erbar'] ],
        'avrs': [ ['ts', 'point-s', [s1]] ],
        'dds': [dd]
    }
    erbar_dict = choose_vars(oo_erbar)[0]
    erbar = erbar_dict['data'][0]
    erbar = np.interp(t, erbar_dict['x'], erbar)

    # find ZFZF and GAM amplitudes in zonal electri field:
    res_e = obtain_e0e1(g, erbar, t, t_point_inf)

    # - Exclude ZFZF from electric field, field energy and energy transfer signals -
    id_peak = 4
    eta = res_e['A_zfzf'] / res_e['A_gam']
    res_f = obtain_gam(w, g, Efield, P, t, id_peak, eta)
    id_peak = res_f['id_peak']

    # - PLOTTING: electric field -
    curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit('Electric\ field')
    curves.flag_norm = False
    curves.new() \
        .XS(t).YS(erbar) \
        .leg('ZF + GAM')
    curves.new() \
        .XS(t).YS(res_e['e_gam_cos']) \
        .leg('GAM:\ cos')
    curves.new() \
        .XS(res_e['t_peaks']).YS(res_e['e_gam_peaks']) \
        .leg('GAM:\ cos:\ peaks').sty('o')
    curves.new() \
        .XS(t).YS(res_e['e_gam']) \
        .leg('GAM')
    cpr.plot_curves(curves)

    # - PLOTTING: field energy -
    curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit('Field\ energy')
    curves.flag_norm = False
    curves.new().XS(t).YS(Efield).leg('ZF+GAM')
    curves.new() \
        .XS(t[res_f['ids_peaks']]) \
        .YS(Efield[res_f['ids_peaks']]) \
        .leg('peaks').sty('o')
    curves.new() \
        .XS( t[res_f['ids_peaks'][id_peak]] ) \
        .YS( Efield[res_f['ids_peaks'][id_peak]] ) \
        .leg('chosen\ peak').sty('s')
    curves.new().XS(t).YS(res_f['Egam']).leg('GAM')
    cpr.plot_curves(curves)

    # - PLOTTING: energy transfer signal -
    curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit('Field\ energy')
    curves.flag_norm = False
    curves.new().XS(t).YS(P).leg('ZF+GAM')
    curves.new().XS(t).YS(res_f['Pgam']).leg('GAM')
    cpr.plot_curves(curves)

    return


# Plot main 2d figures for the MPR diagnostic:
def MPR_plot_2d(dd, oo):
    return


# Find GAM/EGAM damping/growth rates:
def MPR_gamma(dd, oo_vars, oo_wg, oo_plot):
    # --- DESCRIPTION ---
    # Find GAM/EGAM damping/growth rates:
    # total, for different species, in different velocity domains
    # ----------------------------------------------------------------------------
    # --- INPUT DICTIONARIES ---
    # -> oo_vars - dictionary:
    #   'sel_species': 'total', 'deuterium', 'electrons' etc.;
    #   'mu_domain', 'vpar_domain' - 1d arrays with two elements:
    #           velocity domain, where JE integrated;
    # -> oo_plot - dict.:
    #   't_plot' - 1d array [t_start, t_end]:
    #       time domain, where data will be plotted
    # -> oo_wg - dict.: parameters for MPR.GAM_calc_wg() diagnostic:
    #   't_work' - 1d array [t_start, t_end]:
    #       time interval, where gamma will be calculated
    # ---
    #   'flag_naive_t_peaks' = True: only for naive calculation
    #       take a time interval between peaks to calculate the damping/growth rate,
    #       otherwise, take an arbitrary time interval
    #   'naive_n_periods' - integer:
    #       number of GAM periods that will define a length of the time interval
    #
    # ---
    #   'flag_ZFZF' = True:
    #        eliminate contribution of the zero-frequency zonal flows;
    #   'gam-s1' - float:
    #       radial point at which Erbar should be taken
    #   'gam-w', 'gam-g' - floats:
    #       estimation of the GAM frequency and damping/growth rate
    #       normalization: wci
    #       'gam-w' have to include 2*pi
    #   'zfzf-t1' - float:
    #       time point to estimate ZFZF amplitude in Erbar
    #   'id-efield-peak' - integer:
    #       peak in Efield to exclude ZFZF from Efield and JE
    # ---
    #   'flag_t_right' = True:
    #        vary the right time point to estimate better GAM period;
    #   'right-n' - integer:
    #       number of points during variation
    # ---
    #   'flag_stat' = True: (some printing flags will be switched off)
    #        vary time intervals in a give time range;
    #        build histograms for g, and using Gaussian or Student's
    #            distribution functions find 95% confidence intervals as 1.96 sigma;
    #   'n_samples' - integer:
    #       total number of the variations of calculation time intervals
    #   'min_je_n_periods' - integer:
    #       minimum number of GAM periods to estimate in one time interval

    # --- Parameters ---
    sel_species = oo_vars.get('sel_species', 'total')
    mu_domain   = oo_vars.get('mu_domain', [])
    vpar_domain = oo_vars.get('vpar_domain', [])

    flag_ZFZF = oo_wg.get('flag_ZFZF', False)

    # --- Names ---
    je_name = 'jdote_es'
    ef_name = 'efield'

    # --- GET JE and Efield (t) ---
    oo_vvar = {
        'avrs': [['t']],
        'dds': [dd],
    }

    ovar_array = [['mpr', je_name, sel_species, mu_domain, vpar_domain]]
    oo_vvar.update({'ovars': ovar_array})
    je_dict = choose_vars(oo_vvar)[0]

    ovar_array[0][1] = ef_name
    oo_vvar.update({'ovars': ovar_array})
    ef_dict = choose_vars(oo_vvar)[0]

    del oo_vvar, ovar_array

    # --- INTERPOLATE Efield to the t-axis of JE ---
    ef_dict['data'] = np.interp(je_dict['x'], ef_dict['x'], ef_dict['data'])
    ef_dict['x'] = je_dict['x']

    oo_wg.update({
        'je_dict': je_dict,
        'ef_dict': ef_dict
    })

    # --- EXCLUDE ZFZF from signals ---
    if flag_ZFZF:
        t  = je_dict['x']
        je = je_dict['data']
        ef = ef_dict['data']

        s1 = oo_wg.get('gam-s1', None)
        t_point_inf = oo_wg.get('zfzf-t1', None)
        gam_w = oo_wg.get('gam-w', None)
        gam_g = oo_wg.get('gam-g', None)
        id_peak = oo_wg.get('id-efield-peak', None)

        # - find Erbar -
        oo_erbar = {
            'ovars': [['zonal', 'erbar']],
            'avrs': [['ts', 'point-s', [s1]]],
            'dds': [dd]
        }
        erbar_dict = choose_vars(oo_erbar)[0]
        erbar = np.interp(t, erbar_dict['x'], erbar_dict['data'][0])

        # - estimate GAM, ZFZF amplitdes in Erbar -
        res_e = estimate_GAM_ZFZF_amplitudes(gam_g, erbar, t, t_point_inf)

        # - Exclude ZFZF from Ef and JE -
        eta = res_e['A_zfzf'] / res_e['A_gam']
        res_f = MPR_exclude_ZFZF(gam_w, gam_g, ef, je, t, id_peak, eta)
        id_peak = res_f['id_peak']

        # - PLOTTING: electric field -
        curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit(erbar_dict['tit'])
        curves.flag_norm = False
        curves.new() \
            .XS(t).YS(erbar) \
            .leg('ZF + GAM')
        curves.new() \
            .XS(t).YS(res_e['e_gam_cos']) \
            .leg('GAM:\ cos')
        curves.new() \
            .XS(res_e['t_peaks']).YS(res_e['e_gam_peaks']) \
            .leg('GAM:\ cos:\ peaks').sty('o')
        curves.new() \
            .XS(t).YS(res_e['e_gam']) \
            .leg('GAM')
        cpr.plot_curves(curves)

        # - PLOTTING: field energy -
        curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit(ef_dict['tit'])
        curves.flag_norm = False
        curves.new().XS(t).YS(ef).leg('ZF+GAM')
        curves.new() \
            .XS(t[res_f['ids_peaks']]) \
            .YS(ef[res_f['ids_peaks']]) \
            .leg('peaks').sty('o')
        curves.new() \
            .XS(t[res_f['ids_peaks'][id_peak]]) \
            .YS(ef[res_f['ids_peaks'][id_peak]]) \
            .leg('chosen\ peak').sty('s')
        curves.new().XS(t).YS(res_f['Egam']).leg('GAM')
        cpr.plot_curves(curves)

        # - PLOTTING: energy transfer signal -
        curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit(je_dict['tit'])
        curves.flag_norm = False
        curves.new().XS(t).YS(je).leg('ZF+GAM')
        curves.new().XS(t).YS(res_f['Pgam']).leg('GAM')
        cpr.plot_curves(curves)

        # - UPDATE je and efield signals -
        je_dict_gam = dict(je_dict)
        je_dict_gam.update({
            'data': res_f['Pgam']
        })

        ef_dict_gam = dict(ef_dict)
        ef_dict_gam.update({
            'data': res_f['Egam']
        })

        oo_wg.update({
            'je_dict': je_dict_gam,
            'ef_dict': ef_dict_gam
        })

    # --- CALCULATION ---
    oo_desc = {
        'species': sel_species
    }
    MPR.calc_gamma(oo_wg, oo_plot, oo_desc)

    return


# NEW: test GAM from ZF:
def MPR_test_gam_zf():
    # --- TEST CONTRIBUTION OF ZFs ---
    def obtain_gam(w, g, e, E, t, t_point_inf):
        id_t_inf, t_point_inf, _ = mix.get_ids(t, t_point_inf)
        id_t_inf = id_t_inf[0]

        # electric field
        e0 = np.mean(e[id_t_inf:])
        e_gam = e - e0
        e_gam_cos = e_gam * np.exp(-g * t)

        ids_peaks, e_gam_peaks = scipy.signal.find_peaks(e_gam_cos, height=0)
        e_gam_peaks = e_gam_peaks['peak_heights']
        t_peaks = t[ids_peaks]

        e1 = np.mean(e_gam_peaks)

        E_gam = E - 0.5 * e0 ** 2 - e0 * e1 * np.cos(w * t) * np.exp(g * t)

        res = {
            'E_gam': E_gam, 'e_gam': e_gam,
            'e_gam_cos': e_gam_cos,
            't_peaks': t_peaks, 'e_gam_peaks': e_gam_peaks,
        }

        return res

    # - FREQUENCY and DAMPING RATE -
    w = 3.9e-3
    g = -1.1e-4
    t = np.linspace(0, 2e4, 20000)

    # --- TEST SIGNALS: use electric field itself to eliminate ZFZF ---
    e1 = 1
    e0 = 0.2
    e_zf = e0 + e1 * np.cos(w * t) * np.exp(g * t)
    E_gam = 0.25 * e1 ** 2 * (1 + np.cos(2 * w * t)) * np.exp(2 * g * t)
    E_zf = 0.5 * e0 ** 2 + e0 * e1 * np.cos(w * t) * np.exp(g * t) \
           + 0.25 * e1 ** 2 * (1 + np.cos(2 * w * t)) * np.exp(2 * g * t)

    t_point_inf = 1.5e4
    res = obtain_gam(w, g, e_zf, E_zf, t, t_point_inf)

    curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit('Eelctric field')
    curves.flag_norm = False
    curves.new() \
        .XS(t).YS(e_zf) \
        .leg('ZF + GAM')
    curves.new() \
        .XS(t).YS(res['e_gam_cos']) \
        .leg('GAM:\ cos')
    curves.new() \
        .XS(res['t_peaks']).YS(res['e_gam_peaks']) \
        .leg('GAM:\ cos:\ peaks').sty('o')
    curves.new() \
        .XS(t).YS(res['e_gam']) \
        .leg('GAM')
    cpr.plot_curves(curves)

    curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit('Field\ energy')
    curves.flag_norm = False
    curves.new().XS(t).YS(E_gam).leg('GAM')
    curves.new().XS(t).YS(E_zf).leg('ZF+GAM')
    curves.new().XS(t).YS(res['E_gam']).leg('GAM\ from\ ZF').sty(':')
    cpr.plot_curves(curves)

    # --- TEST SIGNALS: several components on frequency, E0 and E1 ---

    return


# NEW: test GAM from ZF:
def MPR_test_gam_zf_efield():
    # --- TEST CONTRIBUTION OF ZFs ---
    def obtain_gam(w, g, E, t, id_peak, eta):
        ids_peaks, _ = scipy.signal.find_peaks(E, height=0)

        E_peak1 = E[ids_peaks[id_peak]]
        E_peak2 = E[ids_peaks[id_peak+1]]

        if eta > 0:
            idd = np.array([E_peak1, E_peak2]).argmax()
        else:
            # in this case you should take at the beginning where
            # there are still alternating peaks
            idd = np.array([E_peak1, E_peak2]).argmin()
        if idd == 1:
            id_peak += 1

        t_peak = t[ids_peaks[id_peak]]
        Eref   = E[ids_peaks[id_peak]]

        e12 = Eref / (0.5 * eta**2 + eta * np.exp(g*t_peak) + 0.5 * np.exp(2*g*t_peak))
        e1 = np.sqrt(e12)
        e0 = e1 * eta

        Egam = E - 0.5 * e0 ** 2 - e0 * e1 * np.cos(w * t) * np.exp(g * t)

        res = {
            'ids_peaks': ids_peaks,
            'id_peak': id_peak,
            'Egam': Egam
        }
        return res

    # - FREQUENCY and DAMPING RATE -
    w = 3.9e-3
    g = -1.1e-4
    t = np.linspace(0, 2e4, 20000)

    # --- TEST SIGNALS: use electric field itself to eliminate ZFZF ---
    e1 = 1
    e0 = -0.2
    E_gam = 0.25 * e1 ** 2 * (1 + np.cos(2 * w * t)) * np.exp(2 * g * t)
    E_zf = 0.5 * e0 ** 2 + e0 * e1 * np.cos(w * t) * np.exp(g * t) \
           + 0.25 * e1 ** 2 * (1 + np.cos(2 * w * t)) * np.exp(2 * g * t)

    eta = e0 / e1
    id_peak = 10

    res = obtain_gam(w, g, E_zf, t, id_peak, eta)

    curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit('Field\ energy')
    curves.flag_norm = False
    curves.flag_semilogy = False
    curves.new().XS(t).YS(E_gam).leg('GAM')
    curves.new().XS(t).YS(E_zf).leg('ZF+GAM')
    curves.new()\
        .XS(t[res['ids_peaks']])\
        .YS(E_zf[res['ids_peaks']])\
        .leg('peaks').sty('o')
    curves.new() \
        .XS( t[ res['ids_peaks'][res['id_peak']] ]  ) \
        .YS( E_zf[ res['ids_peaks'][res['id_peak']] ] ) \
        .leg('chosen\ peak').sty('s')
    curves.new().XS(t).YS(res['Egam']).leg('GAM\ from\ ZF').sty(':')
    cpr.plot_curves(curves)

    # --- TEST SIGNALS: several components on frequency, E0 and E1 ---

    return


# NEW: estimate GAM and ZFZF amplitudes in Erbar:
def estimate_GAM_ZFZF_amplitudes(g, e, t, t_point_inf):
    # g - estimated damping rate of the GAM
    # e - erbar
    # t - time array
    # t_point_inf - time point somewhere in the middle to estimate amplitude of ZFZF
    id_t_inf, t_point_inf, _ = mix.get_ids(t, t_point_inf)

    A_zfzf = np.mean(e[id_t_inf:])
    e_gam = e - A_zfzf
    e_gam_cos = e_gam * np.exp(-g * t)

    ids_peaks, e_gam_peaks = scipy.signal.find_peaks(e_gam_cos, height=0)
    e_gam_peaks = e_gam_peaks['peak_heights']
    t_peaks = t[ids_peaks]

    A_gam = np.mean(e_gam_peaks)

    res = {
        'A_zfzf': A_zfzf, 'A_gam': A_gam,
        'e_gam': e_gam,
        'e_gam_cos': e_gam_cos,
        't_peaks': t_peaks, 'e_gam_peaks': e_gam_peaks,
    }

    return res


# NEW: exclude ZFZF from Efield and JE:
def MPR_exclude_ZFZF(w, g, E, P, t, id_peak, eta):
    # w, g - estimated GAM frequency and damping rate
    # E - field energy
    # P - J*E
    # t - time array
    # id_peak - peak to estimate properly space averaged Erbar_gam and Erbar_zfzf
    # eta - Erbar_zfzf(s1) / Erbar_gam(s1)

    ids_peaks, _ = scipy.signal.find_peaks(E, height=0)

    E_peak1 = E[ids_peaks[id_peak]]
    E_peak2 = E[ids_peaks[id_peak+1]]
    if eta > 0:
        idd = np.array([E_peak1, E_peak2]).argmax()
    else:
        # in this case you should take the peak at the beginning where
        # there are still alternating peaks
        idd = np.array([E_peak1, E_peak2]).argmin()
    if idd == 1:
        id_peak += 1

    t_peak = t[ids_peaks[id_peak]]
    Eref   = E[ids_peaks[id_peak]]

    e12 = Eref / (0.5 * eta**2 + eta * np.exp(g*t_peak) + 0.5 * np.exp(2*g*t_peak))
    e1 = np.sqrt(e12)
    e0 = e1 * eta

    Egam = E - 0.5 * e0 ** 2 - e0 * e1 * np.cos(w * t) * np.exp(g * t)
    Pgam = P + e0*e1 * np.exp(g*t) * (g * np.cos(w*t) - w * np.sin(w*t))

    res = {
        'ids_peaks': ids_peaks,
        'id_peak': id_peak,
        'Egam': Egam,
        'Pgam': Pgam
    }
    return res


