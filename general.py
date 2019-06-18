import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import gam_theory
import gam_exp
import write_data as wr
import numpy as np
from scipy import interpolate
from scipy.signal import correlate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(gam_theory)
    mix.reload_module(gam_exp)
    mix.reload_module(wr)


# plot time evolutino of a signal in several radial points:
def plot_s1(dd, oo):
    # additional data:
    flag_save = oo.get('flag_save', False)
    save_name = oo.get('save_name', None)

    project_name = dd.get('project_name', '')
    var_name = oo.get('var_name', '')

    vvar, s, t, tit = oo['var'], oo['s'], oo['t'], oo['tit']

    labx = oo.get('labx', 't[\omega_{ci}^{-1}]')  # x-label
    flag_norm = oo.get('flag_norm', False)

    # basic legend:
    line_leg_0 = ''
    if len(project_name) > 0:
        line_leg_0 += project_name + ':\ '
    if var_name is not '':
        line_leg_0 += var_name + ':\ '

    # additional data to compare with:
    list_curves_load = oo.get('list_curves_load_s1', {None})

    # radial points:
    s_points = np.array(oo.get('s_points', [0.5]))
    ns_points = np.size(s_points)
    ids_s = np.zeros(ns_points, dtype=np.int64)
    lines_s = []
    for id_s in range(ns_points):
        s1 = s_points[id_s]
        ids_s[id_s], s_points[id_s] = mix.find(s, s1)
        line_s1 = 's = {:0.3f}'.format(s_points[id_s])
        lines_s.append(line_s1)

    # time interval:
    t, ids_t = mix.get_array_oo(oo, t, 't')

    # plotting
    curves = crv.Curves().xlab(labx).tit(tit)
    curves.flag_norm = flag_norm
    for id_s in range(ns_points):
        vvar_s1 = mix.get_slice(vvar, ids_t, ids_s[id_s])
        line_s1 = lines_s[id_s]
        curves.new() \
            .XS(t) \
            .YS(vvar_s1) \
            .leg(line_leg_0 + line_s1)
    for curves_load in list_curves_load:
        curves.load(curves_load)
    curves.set_colors_styles()
    cpr.plot_curves(curves)

    # save data
    if flag_save:
        if 'saved_data' not in dd:
            dd['saved_data'] = {}
        dd['saved_data'][save_name + '-s1'] = curves


# plot radial structure of signals at several time moments:
def plot_t1_vars(dd, oo):
    # additional data:
    project_name = dd.get('project_name', '')

    vvars, ss, ts, tit_vars = oo['vars'], oo['ss'], oo['ts'], oo['tit_vars']
    nvars = len(vvars)

    labx = oo.get('labx', 's')
    flag_norm = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)

    # time points, common for all signals:
    t_points = np.array(oo.get('t_points', [0.0]))
    nt_points = np.size(t_points)

    # basic legend:
    line_leg_0 = ''
    # if len(project_name) > 0:
    #     line_leg_0 += project_name + ':\ '

    # additional data to compare with:
    list_curves_load = oo.get('list_curves_load_t1', {None})

    # find radial structures and plot them:
    curves = crv.Curves().xlab(labx).tit(project_name)
    curves.flag_norm     = flag_norm
    curves.flag_semilogy = flag_semilogy
    for id_var in range(nvars):
        vvar, s, t, tit = vvars[id_var], ss[id_var], ts[id_var], tit_vars[id_var]
        line_leg = line_leg_0 + tit + ':\ '

        if len(np.shape(vvar)) == 1:
            curves.new() \
                .XS(s) \
                .YS(vvar) \
                .leg(line_leg)
            continue

        # time points:
        ids_t = np.zeros(nt_points, dtype=np.int64)
        lines_t = []
        for id_t in range(nt_points):
            t1 = t_points[id_t]
            ids_t[id_t], t_points[id_t] = mix.find(t, t1)
            line_t1 = 't[\omega_{ci}^{-1}]' + ' = {:0.3e}'.format(t_points[id_t])
            lines_t.append(line_t1)

        # radial interval:
        s, ids_s = mix.get_array_oo(oo, s, 's')

        for id_t in range(nt_points):
            vvar_t1 = mix.get_slice(vvar, ids_t[id_t], ids_s)
            line_t1 = lines_t[id_t]
            curves.new() \
                .XS(s) \
                .YS(vvar_t1) \
                .leg(line_leg + line_t1)

    for curves_load in list_curves_load:
        curves.load(curves_load)
    curves.set_colors_styles()
    cpr.plot_curves(curves)


# plot radial structure of a signal at several time moments:
def plot_t1(dd, oo):
    # additional data:
    flag_save = oo.get('flag_save', False)
    save_name = oo.get('save_name', None)

    project_name = dd.get('project_name', '')
    var_name = oo.get('var_name', '')

    vvar, s, t, tit = oo['var'], oo['s'], oo['t'], oo['tit']

    labx = oo.get('labx', 's')  # x-label
    flag_norm = oo.get('flag_norm', False)

    # basic legend:
    line_leg_0 = ''
    if len(project_name) > 0:
        line_leg_0 += project_name + ':\ '
    if var_name is not '':
        line_leg_0 += var_name + ':\ '

    # additional data to compare with:
    list_curves_load = oo.get('list_curves_load_t1', {None})

    # time points:
    t_points = np.array(oo.get('t_points', [0.0]))
    nt_points = np.size(t_points)
    ids_t = np.zeros(nt_points, dtype=np.int64)
    lines_t = []
    for id_t in range(nt_points):
        t1 = t_points[id_t]
        ids_t[id_t], t_points[id_t] = mix.find(t, t1)
        line_t1 = 't[\omega_{ci}^{-1}]' + ' = {:0.3e}'.format(t_points[id_t])
        lines_t.append(line_t1)

    # radial interval:
    s, ids_s = mix.get_array_oo(oo, s, 's')

    # plotting
    curves = crv.Curves().xlab(labx).tit(tit)
    curves.flag_norm = flag_norm
    for id_t in range(nt_points):
        vvar_t1 = mix.get_slice(vvar, ids_t[id_t], ids_s)
        line_t1 = lines_t[id_t]
        curves.new() \
            .XS(s) \
            .YS(vvar_t1) \
            .leg(line_leg_0 + line_t1)
    for curves_load in list_curves_load:
        curves.load(curves_load)
    curves.set_colors_styles()
    cpr.plot_curves(curves)

    # save data
    if flag_save:
        if 'saved_data' not in dd:
            dd['saved_data'] = {}
        dd['saved_data'][save_name + '-t1'] = curves


# plot signals, averaged in space
def plot_avs(dd, oo):
    # additional data:
    flag_save = oo.get('flag_save', False)
    save_name = oo.get('save_name', None)

    project_name = dd.get('project_name', '')
    if len(project_name) > 0:
        project_name = project_name + ':\ '

    labx = oo.get('labx', 't[wci^{-1}]')  # x-label
    tit = oo.get('tit', '')  # title
    vars_names = oo.get('vars_names', [])  # names of the signals

    flag_norm = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)

    # additional data to compare with:
    list_curves_load_avs = oo.get('list_curves_load_avs', {None})

    # find averaging in space:
    avs = mix.find_avs(oo)
    n_avs = len(avs)  # number of signals

    # plotting:
    curves = crv.Curves().xlab(labx).tit(tit)
    curves.flag_norm = flag_norm
    curves.flag_semilogy = flag_semilogy
    for ivar in range(n_avs):
        t = avs[ivar]['t']

        if len(vars_names) == 0:
            var_name = ''
        else:
            var_name = vars_names[ivar]
            if var_name is not '':
                var_name = var_name + ':\ '

        vvar_s, lines_avs = avs[ivar]['data'], avs[ivar]['lines_avs']
        ns_av = len(lines_avs)
        for is_av in range(ns_av):
            vvar, line_av = vvar_s[is_av], lines_avs[is_av]
            curves.new() \
                .XS(t) \
                .YS(vvar)

            if len(project_name + var_name) == 0:
                curves.list_curves[-1] \
                    .leg(line_av)
            else:
                curves.list_curves[-1] \
                    .leg(project_name + var_name) \
                    .legn(line_av)
    for curves_load_avs in list_curves_load_avs:
        curves.load(curves_load_avs)
    curves.set_colors_styles()

    if len(curves.list_curves) is not 0:
        cpr.plot_curves(curves)

    # save data
    if flag_save:
        if 'saved_data' not in dd:
            dd['saved_data'] = {}
        dd['saved_data'][save_name + '-avs'] = curves


# plot signals, averaged in time
def plot_avt(dd, oo):
    # additional data:
    flag_save = oo.get('flag_save', False)
    save_name = oo.get('save_name', None)

    project_name = dd.get('project_name', '')
    if project_name is not '':
        project_name = project_name + ':\ '

    labx = oo.get('labx', 's')  # x-label
    tit = oo.get('tit', '')  # title
    vars_names = oo.get('vars_names', [])  # names of the signals

    flag_norm = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)

    # additional data to compare with:
    list_curves_load_avt = oo.get('list_curves_load_avt', {None})

    # find averaging in time:
    avt = mix.find_avt(oo)
    n_avt = len(avt)  # number of signals

    # plotting:
    curves = crv.Curves().xlab(labx).tit(tit)
    curves.flag_norm = flag_norm
    curves.flag_semilogy = flag_semilogy
    for ivar in range(n_avt):
        s = avt[ivar]['s']

        if len(vars_names) == 0:
            var_name = ''
        else:
            var_name = vars_names[ivar]
            if var_name is not '':
                var_name = var_name + ':\ '

        vvar_s, lines_avt = avt[ivar]['data'], avt[ivar]['lines_avt']
        ns_av = len(lines_avt)
        for it_av in range(ns_av):
            vvar, line_av = vvar_s[it_av], lines_avt[it_av]
            curves.new() \
                .XS(s) \
                .YS(vvar)

            if len(project_name + var_name) == 0:
                curves.list_curves[-1]\
                    .leg(line_av)
            else:
                curves.list_curves[-1] \
                    .leg(project_name + var_name) \
                    .legn(line_av)
    for curves_load_avt in list_curves_load_avt:
        curves.load(curves_load_avt)
    curves.set_colors_styles()

    if len(curves.list_curves) is not 0:
        cpr.plot_curves(curves)

    # save data
    if flag_save:
        if 'saved_data' not in dd:
            dd['saved_data'] = {}
        dd['saved_data'][save_name + '-avt'] = curves


# plot (t,s) structure:
def plot_st(dd, oo):
    vvar = oo.get('var', [])  # signal (t,s)
    t = oo.get('t', [])
    s = oo.get('s', [])
    labx = oo.get('labx', 't[wci^{-1}]')
    laby = oo.get('laby', 's')
    tit = oo.get('tit', '')

    project_name = dd.get('project_name', '')
    if project_name is not '':
        project_name = project_name + ':\ '

    # intervals
    t, ids_t = mix.get_array_oo(oo, t, 't')
    s, ids_s = mix.get_array_oo(oo, s, 's')
    vvar = mix.get_slice(vvar, ids_t, ids_s)

    # plotting
    curves = crv.Curves().xlab(labx).ylab(laby).tit(project_name + tit)
    curves.new().XS(t).YS(s).ZS(vvar).lev(60)
    cpr.plot_curves_3d(curves)


# plot (t,chi) structure:
def plot_chit(dd, oo):
    vvar = oo.get('var', [])  # signal (t,s)
    t = oo.get('t', [])
    chi = oo.get('chi', [])
    labx = oo.get('labx', 't[wci^{-1}]')
    laby = oo.get('laby', '\chi')
    tit = oo.get('tit', '')

    project_name = dd.get('project_name', '')
    if project_name is not '':
        project_name = project_name + ':\ '

    # intervals
    t, ids_t = mix.get_array_oo(oo, t, 't')
    chi, ids_chi = mix.get_array_oo(oo, chi, 'chi')
    vvar = mix.get_slice(vvar, ids_t, ids_chi)

    # plotting
    curves = crv.Curves().xlab(labx).ylab(laby).tit(project_name + tit)
    curves.new().XS(t).YS(chi).ZS(vvar).lev(60)
    cpr.plot_curves_3d(curves)


# plot fft:
def plot_fft(dd, oo):
    sel_norm = oo.get('sel_norm', 'wci')  # -> 'wci', 'kHz', 'csa', 'csr'
    sel_cmp = oo.get('sel_cmp', 'pink_r')  # -> 'jet', 'hot' etc.
    sel_r = oo.get('sel_r', 's')  # -> 's', 'psi'
    curves_to_load = oo.get('curves_to_load', None)

    flag_gao = oo.get('flag_gao', True)
    flag_gk_fit = oo.get('flag_gk_fit', False)
    flag_aug20787 = oo.get('flag_aug20787', False)

    labr = oo.get('labr', 's')
    tit = oo.get('tit', '')

    vvar = oo.get('var', [])
    t = oo.get('t_wci', [])  # should be normalized to wci^{-1}
    r = oo.get('r', [])

    project_name = dd.get('project_name', '')
    if project_name is not '':
        project_name = project_name + ':\ '

    # time interval
    t, ids_t = mix.get_array_oo(oo, t, 't')

    # radial coordinate normalization
    r, ids_r = mix.get_array_oo(oo, r, sel_r)

    # zonal radial electric field
    vvar = mix.get_slice(vvar, ids_t, ids_r)

    # information about the time-interval where FFT is found
    line_t_wci = mix.test_array(t, 't[wci^{-1}]', ':0.2e')
    line_t_seconds = mix.test_array(t / dd['wc'], 't(seconds)', ':0.3e')
    line_t = line_t_wci

    # change normalization of the time grid -> seconds
    t = t / dd['wc']

    # --- FIND FFT ---
    ffvar = ymath.fft_y(t, vvar, oo={'axis': 0})

    # frequency normalization
    coef_norm, line_w = None, ''
    if sel_norm == 'khz':
        coef_norm = 1. / 1.e3
        line_w = '\omega,\ kHz'
    if sel_norm == 'wc':
        coef_norm = 2 * np.pi / dd['wc']
        line_w = '\omega[\omega_c]'
    if sel_norm == 'csa':
        coef_norm = 2 * np.pi / (dd['cs'] / dd['a0'])
        line_w = '\omega[c_s/a_0]'
    if sel_norm == 'csr':
        coef_norm = 2 * np.pi / (dd['cs'] / dd['R0'])
        line_w = '\omega[c_s/R_0]'
    w = ffvar['w'] * coef_norm

    # intervals for the frequency:
    w, ids_w = mix.get_array_oo(oo, w, 'w')
    f_var = mix.get_slice(ffvar['f'], ids_w)

    # --- PLOTTING ---
    curves = crv.Curves().xlab(labr).ylab(line_w) \
        .tit(project_name + tit) \
        .titn(line_t) \
        .xlim([r[0], r[-1]])\
        .ylim([w[0], w[-1]]) \
        .leg_pos('lower left')
    curves.new()\
        .XS(r)\
        .YS(w)\
        .ZS(f_var.T)\
        .lev(20).cmp(sel_cmp)

    # parameters for the analytical and experimental data:
    oo_th = {'curves': curves, 'sel_norm': sel_norm,
             'sel_r': sel_r, 'r': r, 'col': 'white'}

    if flag_aug20787:
        curves = gam_exp.exp_AUG20787(dd, oo_th)
    if flag_gao:
        curves = gam_theory.get_gao(dd, oo_th)
    if flag_gk_fit:
        curves = gam_theory.get_gk_fit(dd, oo_th)

    # load saved data:
    curves.load(curves_to_load)

    # plot curves
    cpr.plot_curves_3d(curves)


# plot fft 1d:
def plot_fft_1d(dd, oo):
    # --- load necessary parameters ---
    vvars = oo.get('vars', [])
    nvars = len(vvars)

    ts = oo.get('ts_wci', [])
    ss = oo.get('ss', [])

    s_points, s_intervals = np.nan, np.nan
    ns = 0

    flag_points = oo.get('flag_points', True)
    if flag_points:
        s_points = oo.get('s_points', [])
        ns = ns + np.size(s_points)

    flag_intervals = oo.get('flag_intervals', False)
    flag_aver = None
    if flag_intervals:
        flag_aver = oo.get('flag_aver', 'mean')
        s_intervals = oo.get('s_intervals', [[0.0, 1.0]])
        ns = ns + np.shape(s_intervals)[0]

    sel_norm = oo.get('sel_norm', 'wc')

    tit = oo.get('tit', '')
    vars_names = oo.get('vars_names', [])
    flag_save = oo.get('flag_save', False)
    save_name = oo.get('save_name', None)

    flag_norm = oo.get('flag_norm', False)

    project_name = dd.get('project_name', '')
    if project_name is not '':
        project_name = project_name + ':\ '

    list_curves_load = oo.get('list_curves_load', {None})

    # --- frequency normalization ---
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

    # --- build 1d fft ---
    ws, ffvvars, lines_ss, lines_t = [], [], [], []
    for ivar in range(nvars):
        vvar = vvars[ivar]
        s, t = ss[ivar], ts[ivar]
        t, ids_t = mix.get_array_oo(oo, t, 't')
        lines_t.append('t[\omega_{ci}^{-1}]' + ' = [{:0.2e}, {:0.2e}]'.format(t[0], t[-1]))

        vvar_s, lines_s = [], []
        if flag_points:
            for s1 in s_points:
                id_s1, s1_res = mix.find(s, s1)
                lines_s.append('s = {:0.3f}'.format(s1_res))
                vvar_s.append(mix.get_slice(vvar, ids_t, id_s1))

        if flag_intervals:
            ids_s_intervals, _, lines_s_int = mix.get_interval(s, s_intervals, 's', '0.3f')
            count_int = -1
            for ids_s in ids_s_intervals:
                count_int += 1

                vvar_one = mix.get_slice(vvar, ids_t, ids_s)
                if flag_aver == 'mean':
                    vvar_one = np.mean(vvar_one, axis=1)
                if flag_aver == 'rms':
                    vvar_one = np.sqrt(np.mean(vvar_one**2, axis=1))

                vvar_s.append(vvar_one)
                lines_s.append(lines_s_int[count_int])

        # - frequency grid -
        ffres = ymath.fft_y(t)
        w = ffres['w'] * coef_norm

        # - FFT -
        ffvvar = []
        for id_s1 in range(ns):
            ffvvar.append(ymath.fft_y(t, vvar_s[id_s1]))

        # intervals for the frequency:
        w, ids_w = mix.get_array_oo(oo, w, 'w')
        for id_s1 in range(ns):
            ffvvar[id_s1]['f'] = mix.get_slice(ffvvar[id_s1]['f'], ids_w)

        # - save calculated data for every signal -
        ws.append(w)
        ffvvars.append(ffvvar)
        lines_ss.append(lines_s)

    # --- plot fourier spectrum ---
    curves = crv.Curves().xlab(line_w)\
        .tit('FFT:\ ' + tit)
    curves.flag_norm = flag_norm
    for ivar in range(nvars):
        w, ffvvar, lines_s = ws[ivar], ffvvars[ivar], lines_ss[ivar]
        var_name = vars_names[ivar]
        for id_s1 in range(ns):
            curves.new() \
                .XS(w) \
                .YS(ffvvar[id_s1]['f']) \

            if len(project_name + var_name) == 0:
                curves.list_curves[-1] \
                    .leg(lines_s[id_s1])
            else:
                curves.list_curves[-1] \
                    .leg(project_name + var_name) \
                    .legn(lines_s[id_s1])
            curves.list_curves[-1].legn(lines_t[id_s1])
    for curves_load in list_curves_load:
        curves.load(curves_load)
    curves.set_colors_styles()
    cpr.plot_curves(curves)

    # --- save data ---
    if flag_save:
        if 'saved_data' not in dd:
            dd['saved_data'] = {}
        dd['saved_data'][save_name + '-fft1d'] = curves


# find time delay between two signals:
def find_time_delay(dd, oo):
    # --- load necessary parameters ---
    y1 = oo.get('var1', [])
    y2 = oo.get('var2', [])
    grid_t1 = oo.get('grid_t1', [])  # must be in wc^{-1} units
    grid_t2 = oo.get('grid_t2', [])  # must be in wc^{-1} units

    flag_norm = oo.get('flag_norm', False)
    sel_norm = oo.get('sel_norm', 'wc')
    vars_names = oo.get('vars_names', [])

    # normalize if necessary the signals:
    if flag_norm:
        y1 = ymath.find_norm(y1)
        y2 = ymath.find_norm(y2)

    # interpolate the second signal to the first time grid:
    t = grid_t1
    f_interp = interpolate.interp1d(grid_t2, y2)
    y2_interp = f_interp(t)

    # time interval
    t, ids_t = mix.get_array_oo(oo, t, 't')
    y1        = mix.get_slice(y1, ids_t)
    y2_interp = mix.get_slice(y2_interp, ids_t)

    grid_t2, ids_t = mix.get_array_oo(oo, grid_t2, 't')
    y2 = mix.get_slice(y2, ids_t)
    del ids_t

    # --- time normalization ---
    coef_norm, line_t = None, ''
    if sel_norm == 'ms':
        coef_norm = 1. / (dd['wc'] * 1e3)
        line_t = 't,\ ms'
    if sel_norm == 'wc':
        coef_norm = 1
        line_t = 't[\omega_c^{-1}]'
    if sel_norm == 'csa':
        coef_norm = (dd['cs'] / dd['a0']) / dd['wc']
        line_t = 't[(c_s/a_0)^{-1}]'
    if sel_norm == 'csr':
        coef_norm = (dd['cs'] / dd['R0']) / dd['wc']
        line_t = 't[(c_s/R_0)^{-1}]'
    t = t * coef_norm
    grid_t2 = grid_t2 * coef_norm

    # plot signals:
    curves = crv.Curves().xlab(line_t).tit('signals')
    curves.flag_norm = True
    curves.new().XS(t).YS(y1).leg(vars_names[0])
    curves.new().XS(t).YS(y2_interp).leg(vars_names[1])
    # curves.new().XS(grid_t2).YS(y2).leg('orig.\ ' + vars_names[1]).sty('--')
    cpr.plot_curves(curves)

    # find correlation:
    # *** positive shift_calculated means that y2 is ahead of y1 ***
    # *** in general, y1(t + shift_calculated) = y2(t) ***
    nt = np.size(t)
    ids_for_correlation = np.arange(1 - nt, nt)
    coef_norm = 1.0 * (t[-1] - t[0]) / nt

    cross_correlation = correlate(y1, y2_interp)
    shift_calculated = coef_norm * ids_for_correlation[cross_correlation.argmax()]
    print(line_t + ' = ' + '{:0.3e}'.format(shift_calculated))

    # plot correlation:
    curves = crv.Curves().xlab(line_t).ylab('a.u')
    curves.new() \
        .XS(coef_norm * ids_for_correlation) \
        .YS(cross_correlation)
    cpr.plot_curves(curves)


# find growth rate:
def find_gamma(dd, oo):
    # additional data:
    project_name = dd.get('project_name', '')
    if len(project_name) > 0:
        project_name = project_name + ':\ '

    labx = oo.get('labx', 't[wci^{-1}]')  # x-label
    tit = oo.get('tit', '')  # title
    var_name = oo.get('var_name', '')  # names of the signals
    flag_norm = oo.get('flag_norm', False)

    vvar = oo.get('var', [])
    t = oo.get('t', [])
    s = oo.get('s', [])
    opt_av = oo.get('opt_av', 'rms')
    t_work_domain = oo.get('t_work_domain', [t[0], t[-1]])

    # find averaging in space:
    data, line_s = [], ''
    if opt_av == 'rms' or opt_av == 'mean':
        oo_avs = dict(oo)
        oo_avs.update({
            'vars': [vvar], 'ts': [t], 'ss': [s],
            'opts_av': [opt_av],
            'tit': tit, 'vars_names': [var_name]
        })
        res = mix.find_avs(oo_avs)[0]
        data, line_s = res['data'][0], res['lines_avs'][0]

    # choose time domain:
    ids_t_work, t_work, line_t_work = \
        mix.get_interval(t, [t_work_domain], 't[\omega_{ci}^{-1}]', '0.3e')
    ids_t_work, t_work, line_t_work = ids_t_work[0], t_work[0], line_t_work[0]
    data_work = mix.get_slice(data, ids_t_work)

    # plotting: chosen time interval:
    curves = crv.Curves().xlab(labx).ylab(var_name).tit(project_name + tit)
    curves.flag_norm = flag_norm
    curves.new() \
        .XS(t) \
        .YS(data)\
        .leg(line_s)
    curves.new() \
        .XS(t_work) \
        .YS(data_work) \
        .leg('work\ domain')

    curves.set_colors_styles()
    cpr.plot_curves(curves)

    # find gamma in this time domain:
    g_est = ymath.estimate_g(t_work, data_work)
    line_g = '\gamma[\omega_{ci}]'
    t_fit = g_est['x_fit']

    # plot fitted signals
    curves = crv.Curves().xlab(labx).ylab(var_name).tit(project_name + tit)
    curves.flag_semilogy = True
    curves.new() \
        .XS(t_work) \
        .YS(data_work) \
        .leg('signal').col('blue')
    curves.new() \
        .XS(g_est['x_peaks']) \
        .YS(g_est['y_peaks']) \
        .leg('peaks').sty('o').col('green')
    curves.new() \
        .XS(t_fit) \
        .YS(g_est['y_fit']) \
        .leg('fitting').col('red').sty('--')
    cpr.plot_curves(curves)

    # Print growth rate:
    print('--- Estimation ---')
    print('*** Growth rate: ' + line_s + ':\ ' + line_t_work + ' ***')
    print('E -> ' + line_g + ' = {:0.3e}'.format(g_est['g']))
    return


# find dynamic rate for several signals and in several time intervals
def find_gamma_adv(dd, oo):
    # several signals -> one can build dependencies between theses signals
    # several time intervals -> so far, for more informative plots and more compact results

    project_name = dd.get('project_name', '')
    if len(project_name) > 0:
        project_name = project_name + ':\ '

    labx = oo.get('labx', 't[wci^{-1}]')
    flag_semilogy = oo.get('flag_semilogy', True)

    vvars = oo.get('vars', [])
    nvars = len(vvars)
    ts = oo.get('ts', [])
    tit_vars = oo.get('tit_vars', [])
    var_names = oo.get('var_names', [])
    lines_s = oo.get('lines_s', [])

    cond_vars = oo.get('cond_vars', {})

    line_g = '\gamma[\omega_{ci}]'

    # read working time domains for every signal:
    t_work_domains_vars = []
    for id_var in range(nvars):
        line_id_var = '{:d}'.format(id_var + 1)
        t_work_domains = oo.get('t_work_domains' + line_id_var, None)
        t_work_domains_vars.append(t_work_domains)

    # --- find dynamic rates ---
    vvars_work_vars, ts_work_vars, lines_t_work_vars, gs_est_vars = [], [], [], []
    for id_var in range(nvars):
        # current data
        data = vvars[id_var]
        t = ts[id_var]
        t_work_domains = t_work_domains_vars[id_var]

        vvars_work, ts_work, lines_t_work, gs_est = [], [], [], []
        for id_t_interval in range(len(t_work_domains)):
            t_work_domain = t_work_domains[id_t_interval]

            # time domain where the gamma will be computed
            ids_t_work, t_work, line_t_work = \
                mix.get_interval(t, [t_work_domain], 't[\omega_{ci}^{-1}]', '0.3e')
            ids_t_work, t_work, line_t_work = ids_t_work[0], t_work[0], line_t_work[0]
            data_work = mix.get_slice(data, ids_t_work)

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
    styles_loc = ['--', '-.', ':']

    line_title = project_name
    curves_work = crv.Curves().xlab(labx).tit(line_title)
    curves_work.flag_semilogy = flag_semilogy
    curves_fit = crv.Curves().xlab(labx).tit(line_title)
    curves_fit.flag_semilogy = flag_semilogy
    for id_var in range(nvars):
        line_id_var = '{:d}'.format(id_var + 1)
        if len(var_names) is 0:
            var_names.append(line_id_var)

        data = vvars[id_var]
        t = ts[id_var]
        tit_var = tit_vars[id_var]
        line_s = lines_s[id_var]

        vvars_work   = vvars_work_vars[id_var]
        ts_work      = ts_work_vars[id_var]
        lines_t_work = lines_t_work_vars[id_var]
        gs_est       = gs_est_vars[id_var]

        # intervals
        t, ids_t = mix.get_array_oo(oo, t, 't')
        data = mix.get_slice(data, ids_t)

        # print information about a current signal:
        print(line_id_var + ': ' + tit_var + ': ' + line_s)

        # plot initial signals
        line_legend_init = line_id_var
        curves_work.new() \
            .XS(t) \
            .YS(data) \
            .leg(line_legend_init)
        curves_fit.new() \
            .XS(t) \
            .YS(data)

        # plot data for every time interval
        for id_t_interval in range(len(ts_work)):
            data_work = vvars_work[id_t_interval]
            t_work = ts_work[id_t_interval]
            line_t_work = lines_t_work[id_t_interval]
            g_est = gs_est[id_t_interval]

            # line_legend_fit   = line_id_var + ':\ ' + line_t_work
            line_legend_fit = var_names[id_var] + ':\ ' + line_t_work

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

            print(line_id_var + ': ' + line_t_work + ': ' + line_g + ' = {:0.3e}'.format(g_est['g']))

    cpr.plot_curves(curves_work)
    cpr.plot_curves(curves_fit)

    return


# save data
def save_data_to_external_file(dd, oo):
    var_name = oo.get('var_name', '')
    path_to_save = oo.get('path_to_save', '')
    var_desc = oo.get('var_desc', '')
    flag_new_file = oo.get('flag_new_file', False)

    # get a necessary signal:
    if var_name in dd:
        vvar = dd[var_name]
    else:
        print('There is not such a signal')
        return

    # save results:
    wr.open_file(path_to_save, flag_new_file)
    wr.save_data_adv(path_to_save, vvar, {'name': var_name, 'desc': var_desc})




























