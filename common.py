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
import MPR as mpr
import numpy as np
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
    mix.reload_module(mpr)


def plot_chit(dd, oo):
    oo_res = dict(oo)
    oo_res.update({
        's_points': [oo['s_point']]
    })

    out = choose_signal_for_comparison(dd, oo)

    oo_chit = dict(oo)
    oo_chit.update(out)
    gn.plot_chit(dd, oo_chit)


def plot_t1(dd, oo):
    # plot radial structures of different signals at the same time moments

    oo_choice = dict(oo)
    oo_choice.update({
        'opt_av': None
    })
    res = choose_signals(dd, oo)

    oo_gamma = dict(oo)
    oo_gamma.update({
        'vars': res['vars'],
        'ts':   res['ts'],
        'ss':   res['ss'],
        'tit_vars': res['tit_vars'],
        'lines_s':  res['lines_avr']
    })
    gn.plot_t1_vars(dd, oo_gamma)


def find_gamma_max_s1(dd, oo):
    # take max along chi of non-zonal phi at s1
    # take zonal er at s1

    # time work domains:
    t_work_domains_z  = oo.get('t_work_domains_z',  [])
    t_work_domains_nz = oo.get('t_work_domains_nz', [])

    # radial point
    s_point = oo.get('s_point', None)

    # zonal signal:
    oo_erbar = dict(oo)
    oo_erbar = {
        'opt_var': 'erbar'
    }
    var_erbar = zf.choose_var(dd, oo_erbar)

    id_s, s_point = mix.find(var_erbar['s'], s_point)
    line_s_zf = 's = {:0.3f}'.format(s_point)
    data_erbar_s1 = var_erbar['var'][:, id_s]

    # non-zonal signal:
    name_var_max = itg.phinz_abs_max_chi(dd, oo)
    var_phinz = dd[name_var_max]

    id_s, s_point = mix.find(var_phinz['s'], s_point)
    line_s_nzf = 's = {:0.3f}'.format(s_point)
    data_phinz_s1 = var_phinz['data'][:, id_s]

    # find dynamic rate
    oo_gamma = dict(oo)
    oo_gamma.update({
        'vars':     [data_erbar_s1,    data_phinz_s1],
        'ts':       [var_erbar['t'],   var_phinz['t']],
        'tit_vars': [var_erbar['tit'], 'max_{\chi}:\ \Phi - overline{\Phi}'],
        'lines_s':  [line_s_zf,        line_s_nzf],
        't_work_domains1': t_work_domains_z,
        't_work_domains2': t_work_domains_nz,
        'var_names': ['\overline{E}_r', '\Phi - \overline{\Phi}_r']
    })
    gn.find_gamma_adv(dd, oo_gamma)


# find dynamic rate for several signals and in several time intervals
def find_gamma_adv(dd, oo):
    res = choose_signals(dd, oo)

    oo_gamma = dict(oo)
    oo_gamma.update({
        'vars':     res['vars'],
        'ts':       res['ts'],
        'tit_vars': res['tit_vars'],
        'lines_s':  res['lines_avr']
    })
    gn.find_gamma_adv(dd, oo_gamma)


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
        vvar = mpr.choose_one_var_vparmu(oovar, dd)

    # save information about coordinate system:
    vvar['x1'], vvar['fx1'], vvar['labx'] = 'mu',   '{:0.3f}', '\mu'
    vvar['x2'], vvar['fx2'], vvar['laby'] = 'vpar', '{:0.3f}', 'v_{\parallel}'

    return vvar


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
