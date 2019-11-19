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
import fields3d
import numpy as np
import scipy.signal
from scipy.stats import norm as stat_norm
import sys
import pywt
from scipy import constants


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
    mix.reload_module(fields3d)


GLOBAL_NONE_FILTER = {'sel_filt': None}


# OLD: find correlation between two signals:
def correlation_two_vars(dd, oo):
    # def aver_var(vvar, s, s_interval, tit_var, flag_aver):
    #     ids_s, _, lines_s = \
    #         mix.get_interval(s, [s_interval], 's', '0.3f')
    #     tit_var += ':\ ' + lines_s[0]
    #
    #     ids_s = ids_s[0]
    #     var_avr = vvar[:, ids_s[0]:ids_s[-1]+1]
    #
    #     if flag_aver == 'mean':
    #         var_avr = np.mean(var_avr, axis=1)
    #     if flag_aver == 'rms':
    #         var_avr = np.sqrt(np.mean(var_avr ** 2, axis=1))
    #     return var_avr, tit_var
    #
    # # --- var1 ---
    # oo_var = dict(oo)
    # oo_var.update({
    #     'opt_var':      oo['opt_var1'],
    #     'opt_type':     oo['opt_type1'],
    #     'species_name': oo.get('species_name1', 'deuterium')
    # })
    # # var1 = choose_signal_for_comparison(dd, oo_var)
    #
    # # --- var2 ---
    # oo_var = dict(oo)
    # oo_var.update({
    #     'opt_var':      oo['opt_var2'],
    #     'opt_type':     oo['opt_type2'],
    #     'species_name': oo.get('species_name2', 'deuterium')
    # })
    # # var2 = choose_signal_for_comparison(dd, oo_var)
    #
    # # --- take the signals in some radial intervals:
    # s_interval1 = oo.get('s_interval1', None)
    # flag_aver1 = oo.get('flag_aver1', 'mean')
    # s_interval2 = oo.get('s_interval2', None)
    # flag_aver2 = oo.get('flag_aver2', 'mean')
    #
    # var1_avr, tit1_var = \
    #     aver_var(var1['var'], var1['s'], s_interval1, var1['tit'], flag_aver1)
    # var2_avr, tit2_var = \
    #     aver_var(var2['var'], var2['s'], s_interval2, var2['tit'], flag_aver2)
    #
    # # --- calculate a time delay ---
    # oo_dt = dict(oo)
    # oo_dt.update({
    #     'var1':    var1_avr,  'var2':    var2_avr,
    #     'grid_t1': var1['t'], 'grid_t2': var2['t'],
    #     'vars_names': [tit1_var, tit2_var]
    # })
    # gn.find_time_delay(dd, oo_dt)
    # n_vars = oo.get('nvars', 1)
    # vvars, ts, tit_vars, lines_s = [], [], [], []
    # for i_var in range(n_vars):
    #     line_id_var = '{:d}'.format(i_var + 1)
    #
    #     # choose a signal
    #     opt_type = oo.get('opt_type' + line_id_var, None)
    #     opt_var = oo.get('opt_var' + line_id_var, None)
    #     species_name = oo.get('species_name' + line_id_var, None)
    #     oo_var = {
    #         'opt_type': opt_type,
    #         'opt_var': opt_var,
    #         'species_name': species_name,
    #     }
    #     res = choose_signal_for_comparison(dd, oo_var)
    #     data, s, t, tit = res['var'], res['s'], res['t'], res['tit']
    #
    #     # average the signal or take it at particular radial points:
    #     opt_av = oo.get('opt_av' + line_id_var, 's1')
    #     line_s = ''
    #     if opt_av == 'rms' or opt_av == 'mean':
    #         oo_avs = dict(oo)
    #         s_interval = oo.get('s_interval' + line_id_var, None)
    #         oo_avs.update({
    #             'vars': [data], 'ts': [t], 'ss': [s],
    #             'opts_av': [opt_av],
    #             'tit': tit, 'vars_names': [tit],
    #             'ns_av': 1,
    #             's_av1_start': s_interval[0], 's_av1_end': s_interval[-1]
    #         })
    #         res = mix.find_avs(oo_avs)[0]
    #         data, line_s = res['data'][0], res['lines_avs'][0]
    #     if opt_av == 's1':
    #         s1 = oo.get('s_point', None)
    #         ids_s, s1 = mix.find(s, s1)
    #         data = data[:, ids_s[0]:ids_s[-1] + 1]
    #         line_s = 's = {:0.3f}'.format(s1)
    #
    #     # save the signal
    #     vvars.append(data)
    #     ts.append(t)
    #     tit_vars.append(tit)
    #     lines_s.append(line_s)

    return


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
    labx = oo.get('labx', None)  # x-label
    laby = oo.get('laby', None)  # y-label
    tit_plot   = oo.get('tit_plot', None)  # title
    legs = oo.get('legs', None)

    dds = oo.get('dds', None)

    flag_norm     = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)
    flag_colorbar = oo.get('flag_colorbar', True)
    sel_norm_x1 = oo.get('sel_norm_x1', 'orig')
    oo_text_vars     = oo.get('text_vars', None)

    # normalization
    coef_x1_norm = 1
    line_x1_norm = ''
    if sel_norm_x1 == 'orig':
        coef_x1_norm = 1
        line_x1_norm = ''
    if sel_norm_x1 == 't-mili-seconds':
        line_x1_norm = '(ms)'

    # plotting:
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        leg = vvar['legs'][0]
        texts_var = oo_text_vars[ivar] if oo_text_vars is not None else []

        # title:
        if legs is not None:
            if len(legs) > ivar:
                leg = legs[ivar]
        if tit_plot is '':
            tit_plot_res = leg
        else:
            tit_plot_res = tit_plot

        # X,Y labels:
        if labx is None:
            labx = vvar.get('labx', None)
        if laby is None:
            laby = vvar.get('laby', None)

        # create curves:
        curves = crv.Curves().xlab(labx + line_x1_norm).ylab(laby).tit(tit_plot_res)
        curves.flag_norm     = flag_norm
        curves.flag_semilogy = flag_semilogy
        curves.flag_colorbar = flag_colorbar

        # original or cartesian 2d grid:
        if vvar['x1'] is 'r' and vvar['x2'] is 'z':
            _, ids_s   = mix.get_array_oo(oo, vvar['s'],   's')
            _, ids_chi = mix.get_array_oo(oo, vvar['chi'], 'chi')
            x1 = mix.get_slice(vvar['r'], ids_chi, ids_s)
            x2 = mix.get_slice(vvar['z'], ids_chi, ids_s)
            data = mix.get_slice(vvars[ivar]['data'], ids_s, ids_chi)
        else:
            x1, ids_x1 = mix.get_array_oo(oo, vvar[vvar['x1']], vvar['x1'])
            x2, ids_x2 = mix.get_array_oo(oo, vvar[vvar['x2']], vvar['x2'])
            data = mix.get_slice(vvars[ivar]['data'], ids_x1, ids_x2)

        # normalization
        dd_one = dds[ivar]
        if sel_norm_x1 == 't-mili-seconds':
            coef_x1_norm = 1. / dd_one['wc'] * 1e3

        # additional text:
        for oo_text in texts_var:
            oText = crv.PlText(oo_text)
            curves.newt(oText)

        curves.new().XS(x1 * coef_x1_norm).YS(x2).ZS(data).lev(60)
        cpr.plot_curves_3d(curves)


# NEW: 1d plot: plot vars along t or s:
def plot_vars_1d(oo):
    # - choose a variable to work with -
    # oo.ovars = [[type, opt_var, opts], [], ...]
    # oo.avrs = [[coords, type-coord_av, domains], [], ...]
    vvars = choose_vars(oo)
    n_vars = len(vvars)

    # - input data -
    labx       = oo.get('labx', None)      # x-label
    laby       = oo.get('laby', None)      # y-label
    leg_pos    = oo.get('leg_pos', None)   # position of a legend
    tit_plot   = oo.get('tit_plot', None)  # title
    stys       = oo.get('stys', None)
    cols       = oo.get('cols', None)
    ylim       = oo.get('ylim', None)
    oo_texts   = oo.get('text', [])
    xticks     = oo.get('xticks', None)

    flag_norm     = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)
    opt_legs      = oo.get('legs', [])

    sel_norm_x = oo.get('sel_norm', 'orig')
    sel_norm_ys = oo.get('sel_norm_ys', None)
    oo_filt = oo.get('oo_filt', {})
    dds = oo.get('dds', None)  # array with projects [dd1, dd2, dd3, ...]

    oo_postprocessing = oo.get('oo_postprocessing', None)  # postprocessing (one operation for every var)

    # normalization (first stage):
    coef_x_norm, line_x_norm = 1, ''
    line_y_norm = ''
    if sel_norm_x == 't-mili-seconds':
        line_x_norm = '(ms)'

    if sel_norm_ys is not None:
        if len(sel_norm_ys) == 1:
            if sel_norm_ys[0] == 'energy-transfer':
                line_y_norm = '\ [kW]'
            if sel_norm_ys[0] == 'energy-J':
                line_y_norm = '\ [J]'
            if sel_norm_ys[0] == 'n':
                line_y_norm = '\ [m^{-3}]'

    # -- PLOTTING --
    if labx is not None:
        labx += line_x_norm
    if laby is not None:
        laby += line_y_norm
    curves = crv.Curves().xlab(labx).ylab(laby).tit(tit_plot)\
        .leg_pos(leg_pos).ylim(ylim)
    curves.flag_norm     = flag_norm
    curves.flag_semilogy = flag_semilogy
    curves.xt(xticks)
    count_line = -1
    count_leg = -1

    # additional text:
    for oo_text in oo_texts:
        oText = crv.PlText(oo_text)
        curves.newt(oText)

    # - different variables -
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        datas = vvar['data']
        x_init = vvar['x']
        legs  = vvar['legs']
        ns_av = len(datas)
        dd_one = dds[ivar]

        oo_var_operations = oo_postprocessing[ivar] \
            if oo_postprocessing is not None else None

        # normalization (second stage):
        if sel_norm_x == 't-mili-seconds':
            coef_x_norm = 1. / dd_one['wc'] * 1e3

        line_leg_norm, coef_y_norm = '', 1
        if sel_norm_ys is not None:
            if sel_norm_ys[ivar] == 'energy-transfer':
                line_leg_norm = '\ [W]'
                # T_speak_eV = dd_one['T_speak'] / constants.elementary_charge
                # coef_y_norm = dd_one['Lx'] / 2. * dd_one['cs'] * T_speak_eV / (dd_one['a0'])
                coef_y_norm = dd_one['T_speak'] * dd_one['wc']
                # coef_y_norm *= 1.e-3
            if sel_norm_ys[ivar] == 'energy-J':
                line_leg_norm = '\ [J]'
                # T_speak_eV = dd_one['T_speak'] / constants.elementary_charge
                # coef_y_norm = dd_one['cs'] * T_speak_eV / (dd_one['a0'] * dd_one['wc'])
                coef_y_norm = dd_one['T_speak']
            if sel_norm_ys[ivar] == 'n':
                line_leg_norm = '\ [m^3]'
                coef_y_norm = equil_profiles.ne_avr_m3(dd_one)

        # - different averaging of the same variable -
        for is_av in range(ns_av):
            count_line += 1
            count_leg += 1

            # get data
            x = np.array(x_init)
            data = np.array(datas[is_av])

            # filtering:
            data, x_filt = global_filter_signal(x, data, oo_filt)
            data = np.interp(x, x_filt, data)

            # - post-processing -
            data, x = post_processing(data, x, oo_var_operations)

            # domain of plotting
            x, ids_x = mix.get_array_oo(oo, x, 'x')
            data = mix.get_slice(data, ids_x)

            # x normalization
            x = x * coef_x_norm
            data = data * coef_y_norm

            # styles
            sty_current = '-'
            if stys is not None:
                if count_line < len(stys):
                    sty_current = stys[count_line]

            # color
            color_current = None
            if cols is not None:
                if count_line < len(cols):
                    color_current = cols[count_line]

            # legends
            one_leg = legs[is_av] + [line_leg_norm] if isinstance(legs[is_av], list) \
                else legs[is_av] + line_leg_norm
            if len(opt_legs) > count_leg:
                one_leg = opt_legs[count_leg]

            # - add a new curve -
            curves.new().XS(x).YS(data).leg(one_leg).sty(sty_current).col(color_current)

    # - plot the curves -
    if len(curves.list_curves) is not 0:
        cpr.plot_curves(curves)


def post_processing(data, x, oo_operations):
    # - additional functions -
    def get_x_data_interval(x_domain_loc, x_loc, data_loc):
        ids_x, _, _ = mix.get_ids(x_loc, x_domain_loc)
        ids_x = range(ids_x[0], ids_x[-1] + 1)
        return  data_loc[ids_x], x_loc[ids_x]

    def incr(ii):
        ii[0] += 1
        return ii[0]

    # - check operations -
    if oo_operations is None:
        return data, x

    # adjust the structure of the oo_operations variable
    if not isinstance(oo_operations[0], list):
        oo_operations = [oo_operations]

    # -- PERFORM CONSEQUENTLY ALL OPERATIONS --
    data_work, x_work = np.array(data), np.array(x)
    for oo_operation in oo_operations:
        # - index of an option in the oo_operations -
        idoo = [-1]

        # - operation type and working domain -
        sel_operation = oo_operation[incr(idoo)]
        x_domain      = oo_operation[incr(idoo)] if len(oo_operation) > 1 else [x_work[0], x_work[-1]]
        if x_domain is None:
            x_domain = [x_work[0], x_work[-1]]
        data_work, x_work = get_x_data_interval(x_domain, x_work, data_work)
        data_temp, x_temp = [], []

        # * integration from beginning till the point for every selected x-steps *
        if sel_operation is None:
            data_temp, x_temp = data_work, x_work
        elif sel_operation == 'integration-accumulation':
            width_domain = oo_operation[incr(idoo)]

            x_start = x_work[0]
            x_end = x_start + width_domain
            while x_end <= x_work[-1]:
                data_curr, x_curr = get_x_data_interval([x_start, x_end], x_work, data_work)
                value_curr = np.trapz(data_curr, x=x_curr)

                data_temp.append(value_curr)
                x_temp.append(x_end)

                x_end += width_domain

        # * integration within every x-step *
        elif sel_operation == 'integration-step':
            width_domain = oo_operation[incr(idoo)]
            step_domain  = oo_operation[incr(idoo)]

            x_start = x_work[0]
            x_end   = x_start + width_domain
            while x_end <= x_work[-1]:
                data_curr, x_curr = get_x_data_interval([x_start, x_end], x_work, data_work)
                value_curr = np.trapz(data_curr, x=x_curr)

                x_point = (x_end + x_start) / 2.
                data_temp.append(value_curr)
                x_temp.append(x_point)

                x_start += step_domain
                x_end    = x_start + width_domain

        # * get peaks of absolute signal *
        elif sel_operation == 'abs_peaks':
            ids_abs_peaks, _ = scipy.signal.find_peaks(np.abs(data_work))
            x_temp    = x_work[ids_abs_peaks]
            data_temp = data_work[ids_abs_peaks]

        # * take absolute signal *
        elif sel_operation == 'abs':
            data_temp, x_temp = np.abs(data_work), x_work

        # * multiply by a value *
        elif sel_operation == 'mult':
            coef = oo_operation[incr(idoo)]
            data_temp, x_temp = coef * data_work, x_work
        else:
            print('Wrong operation name')
            sys.exit(-1)

        # - save modification -
        data_work, x_work = np.array(data_temp), np.array(x_temp)

    return data_work, x_work


# NEW: plot Continuous Wavelet transform
def plot_cwt(oo):
    non_filt = {'sel_filt': None}

    # - FUNCTION: filtering at one stage -
    def several_filters(x, y, oo_filt_loc):
        oo_filts_res = []
        if oo_filt_loc is None or len(oo_filt_loc) is 0:
            oo_filts_res.append(non_filt)
        elif isinstance(oo_filt_loc, dict):
            oo_filts_res.append(oo_filt_loc)
        else:
            oo_filts_res = oo_filt_loc  # array of filters

        # apply one by one the filters in the array :
        dict_loc = {'x': np.array(x), 'filt': np.array(y)}
        count_filt = -1
        for one_filt in oo_filts_res:
            count_filt += 1
            dict_loc = ymath.filtering(
                dict_loc['x'], dict_loc['filt'], one_filt
            )
        return dict_loc

    # -FUNCTION: Filter a signal -
    def filter_signal(t, y, oo_filt):
        # initial fft
        dict_fft_ini = ymath.filtering(t, y, non_filt)
        init_dict = {'x': t, 'data': y}

        # filtering
        res_dict = dict(init_dict)
        filt = several_filters(init_dict['x'], init_dict['data'], oo_filt)
        res_dict['data'] = np.array(filt['filt'])
        res_dict['x']    = np.array(filt['x'])

        return res_dict['data'], res_dict['x']

    vvars = choose_vars(oo)  # should be 1d variable
    n_vars = len(vvars)

    # additional data:
    labx       = oo.get('labx', None)  # x-label
    laby       = oo.get('laby', None)  # y-label
    leg_pos    = oo.get('leg_pos', None)
    tit_plot   = oo.get('tit_plot', None)  # title
    ylim       = oo.get('ylim', None)
    width_t    = oo.get('width_t', None)
    oo_filt    = oo.get('oo_filt', None)
    w_norm_domain = oo.get('w_norm_domain', None)

    domain_fft = oo.get('domain_fft', None)

    dds = oo.get('dds', None)

    # - frequency normalization (notation) -
    sel_norm_w = oo.get('sel_norm_w', 'wc')

    line_norm_w = None
    if sel_norm_w == 'khz':
        line_norm_w = '\omega,\ kHz'
    if sel_norm_w == 'wc':
        line_norm_w = '\omega[\omega_c]'
    if sel_norm_w == 'csa':
        line_norm_w = '\omega[c_s/a_0]'
    if sel_norm_w == 'csr':
        line_norm_w = '\omega[c_s/R_0]'

    # time normalization (notation):
    sel_norm_t = oo.get('sel_norm_t', 'wc')

    coef_norm_t, line_norm_t = None, None
    if sel_norm_t == 'orig':
        coef_norm_t = 1
        line_norm_t = ''
    if sel_norm_t == 't-mili-seconds':
        line_norm_t = '(ms)'

    for ivar in range(n_vars):
        curves_t = crv.Curves()\
            .xlab(labx + line_norm_t)\
            .ylab(laby)\
            .tit(tit_plot) \
            .leg_pos(leg_pos).ylim(ylim)
        curves_w = crv.Curves() \
            .xlab(labx + line_norm_t) \
            .ylab(line_norm_w) \
            .tit('CWT:\ ' + tit_plot)
        curves_fft = crv.Curves() \
            .xlab('\omega(kHz)') \
            .tit('FFT')

        dd_one = dds[ivar]
        vvar = vvars[ivar]
        x, ids_x = mix.get_array_oo(oo, vvar['x'], 'x')
        data = mix.get_slice(vvar['data'][0], ids_x)
        leg  = vvar['legs'][0]

        # filtering
        data_filt, x_filt = filter_signal(x, data, oo_filt)
        data_filt = np.interp(x, x_filt, data_filt)

        # result signal to analyse:
        data_res = data - data_filt

        # find CWT
        ids_x_width, _, _ = mix.get_ids(x, [x[0], x[0] + width_t])
        n_width = np.size(ids_x_width)
        cwr_res, w = ymath.cwt_y(x, data_res, {'width': np.arange(1, n_width)})

        # time normalization (second stage):
        if sel_norm_t == 't-mili-seconds':
            coef_norm_t = 1. / dd_one['wc'] * 1e3
        # - frequency normalization (coefficient) -
        coef_norm_w = None
        if sel_norm_w == 'khz':
            coef_norm_w = dd_one['wc'] / 1.e3
        if sel_norm_w == 'wc':
            coef_norm_w = 2 * np.pi
        if sel_norm_w == 'csa':
            coef_norm_w = 2 * np.pi * dd_one['wc'] / (dd_one['cs'] / dd_one['a0'])
        if sel_norm_w == 'csr':
            coef_norm_w = 2 * np.pi * dd_one['wc'] / (dd_one['cs'] / dd_one['R0'])

        # Fourier transform:
        ids_x_fft, x_fft, _ = mix.get_ids(x, domain_fft)
        ids_x_fft = np.arange(ids_x_fft[0], ids_x_fft[-1] + 1)
        res_fft = ymath.fft_y(x_fft, data_res[ids_x_fft])

        ids_x_fft, x_fft, _ = mix.get_ids(x, np.array(domain_fft) + 3 * domain_fft[0])
        ids_x_fft = np.arange(ids_x_fft[0], ids_x_fft[-1] + 1)
        res_fft2 = ymath.fft_y(x_fft, data_res[ids_x_fft])

        ids_x_fft, x_fft, _ = mix.get_ids(x, np.array(domain_fft) + 6 * domain_fft[0])
        ids_x_fft = np.arange(ids_x_fft[0], ids_x_fft[-1] + 1)
        res_fft3 = ymath.fft_y(x_fft, data_res[ids_x_fft])

        curves_t.new()\
            .XS(x * coef_norm_t)\
            .YS(data)\
            .leg(leg + ':\ original').sty('-')
        curves_t.new() \
            .XS(x * coef_norm_t) \
            .YS(data_filt) \
            .leg(leg + ':\ filtered').sty(':')
        curves_t.new() \
            .XS(x * coef_norm_t) \
            .YS(data_res) \
            .leg(leg + ':\ original - filtered').sty('-')

        w_norm_flipped = np.flip(w * coef_norm_w)
        ids_w_plot, w_norm_plot_flipped, _ = mix.get_ids(w_norm_flipped, w_norm_domain)
        w_norm_plot = np.flip(w_norm_plot_flipped)
        ids_w_plot = np.arange(ids_w_plot[0], ids_w_plot[-1] + 1)

        cwr_res_flipped = np.flip(cwr_res, axis=0)
        cwr_res_plot = np.flip(cwr_res_flipped[ids_w_plot, :])

        curves_w.new()\
            .XS(x * coef_norm_t)\
            .YS(w_norm_plot)\
            .ZS(cwr_res_plot.T)

        curves_fft.new()\
            .XS(res_fft['w'] * dd_one['wc']/1.e3)\
            .YS(res_fft['f'])
        curves_fft.new() \
            .XS(res_fft2['w'] * dd_one['wc'] / 1.e3) \
            .YS(res_fft2['f'])
        curves_fft.new() \
            .XS(res_fft3['w'] * dd_one['wc'] / 1.e3) \
            .YS(res_fft3['f'])

        if len(curves_t.list_curves) is not 0:
            cpr.plot_curves(curves_t)
        cpr.plot_curves_3d(curves_w)
        cpr.plot_curves(curves_fft)


# NEW: plot FFT in time
def plot_fft_in_time(oo):
    vvars = choose_vars(oo)  # should be 1d variable
    n_vars = len(vvars)

    # additional data:
    labx       = oo.get('labx', None)  # x-label
    laby       = oo.get('laby', None)  # y-label
    leg_pos    = oo.get('leg_pos', None)
    tit_plot   = oo.get('tit_plot', None)  # title
    ylim       = oo.get('ylim', None)
    width_x    = oo.get('width_x', None)
    x_fft_domain = oo.get('x_fft_domain', None)  # x-domain, where FFT will be performed, does not influence
                                                 # the domain where filtering is perfomed
    oo_filt    = oo.get('oo_filt', None)  # filter for the signal; to be performed in a whole x-domain
    w_norm_domain = oo.get('w_norm_domain', None)
    dds = oo.get('dds', None)

    flag_fv = oo.get('flag_fv', False)  # to calculate frequency from clumps and holes
    species_name = oo.get('species_name', None)  # species name for distribution function
    s1_q = oo.get('s1_q', None)  # radial point where a value of q to be taken

    # - frequency normalization (notation) -
    sel_norm_w = oo.get('sel_norm_w', 'wc')
    line_norm_w = None
    if sel_norm_w == 'khz':
        line_norm_w = '\omega,\ kHz'
    if sel_norm_w == 'wc':
        line_norm_w = '\omega[\omega_c]'
    if sel_norm_w == 'csa':
        line_norm_w = '\omega[c_s/a_0]'
    if sel_norm_w == 'csr':
        line_norm_w = '\omega[c_s/R_0]'

    # time normalization (notation):
    sel_norm_t = oo.get('sel_norm_t', 'wc')
    coef_norm_t, line_norm_t = None, None
    if sel_norm_t == 'orig':
        coef_norm_t = 1
        line_norm_t = ''
    if sel_norm_t == 't-mili-seconds':
        line_norm_t = '(ms)'

    # consider every variable
    for ivar in range(n_vars):
        curves_t = crv.Curves()\
            .xlab(labx + line_norm_t)\
            .ylab(laby)\
            .tit(tit_plot) \
            .leg_pos(leg_pos).ylim(ylim)
        curves_w = crv.Curves() \
            .xlab(labx + line_norm_t) \
            .ylab(line_norm_w) \
            .tit(tit_plot)

        dd_one = dds[ivar]
        vvar = vvars[ivar]
        x    = vvar['x']
        data = vvar['data'][0]
        leg  = vvar['legs'][0]

        # filtering (performed in the whole domain)
        data_filt, x_filt = global_filter_signal(x, data, oo_filt)
        data_filt = np.interp(x, x_filt, data_filt)

        # shifted signal:
        # data_shifted = data - data_filt
        data_shifted = data

        # x-domain, where the FFT will be performed:
        if x_fft_domain is None:
            x_work = np.array(x)
            data_work = np.array(data_shifted)
        else:
            ids_x_work, x_work, _ = mix.get_ids(x, x_fft_domain)
            ids_x_work = np.arange(ids_x_work[0], ids_x_work[-1] + 1)
            data_work = data_shifted[ids_x_work]

        # get frequency grid
        _, x_interval, _ =  mix.get_ids(x_work, [0, width_x])
        w = ymath.fft_y(x_interval)['w']

        # Calculate time intervals:
        n_intervals = 0
        for x1 in x_work:
            n_intervals += 1
            x_right_bound = x1 + width_x
            if x_right_bound > x_work[-1]:
                continue

        # Time evolution of the Fourier transform:
        res_fft = np.zeros([n_intervals, np.size(w)])
        for id_x in range(np.size(x_work)):
            x_right_bound = x_work[id_x] + width_x
            if x_right_bound > x_work[-1]:
                continue

            ids_x_interval, x_interval, _ = \
                mix.get_ids(x_work, [x_work[id_x], x_right_bound])
            ids_x_interval = np.arange(ids_x_interval[0], ids_x_interval[-1] + 1)
            res_fft[id_x, :] = ymath.fft_y(x_interval, data_work[ids_x_interval])['f']

        # time normalization (second stage):
        if sel_norm_t == 't-mili-seconds':
            coef_norm_t = 1. / dd_one['wc'] * 1e3

        # Clumps and Holes and zeros: parallel velocities, result frequencies ---
        if flag_fv:
            # delta f
            oo_get_f = {'avrs': [['tvpar', 'none-']], 'dds': [dd_one],
                        'ovars': [['distribution', 'df_vel_1d', species_name]]}
            vvar_delta_f = choose_vars(oo_get_f)[0]
            deltaf = vvar_delta_f['data']
            t_f = vvar_delta_f['t']
            vpar_f = vvar_delta_f['vpar']

            # d \delta f / dv
            oo_get_f['ovars'] = [['distribution', 'ddeltaf_vel_1d-dv', species_name]]
            vvar_d_delta_f = choose_vars(oo_get_f)[0]
            ddeltaf = vvar_d_delta_f['data']

            # clumps, holes, zeros:
            ids_pos, _, _ = mix.get_ids(vpar_f, [0.0, np.max(vpar_f[-1])])
            ids_pos = range(ids_pos[0], ids_pos[-1] + 1)
            vpar_f_pos = vpar_f[ids_pos]

            id_f_max = np.argmax(deltaf[:, ids_pos], axis=1)
            id_f_min = np.argmin(deltaf[:, ids_pos], axis=1)
            id_f_zero = np.argmin(ddeltaf[:, ids_pos], axis=1)

            vpar_max, vpar_min, vpar_zero = \
                np.zeros(np.size(id_f_max)), np.zeros(np.size(id_f_min)), np.zeros(np.size(id_f_zero))
            for id_t in range(np.size(id_f_max)):
                vpar_max[id_t] = vpar_f_pos[id_f_max[id_t]]
                vpar_min[id_t] = vpar_f_pos[id_f_min[id_t]]
                vpar_zero[id_t] = vpar_f_pos[id_f_zero[id_t]]
            w_from_max = get_w_from_v_res(dd_one, s1_q, vpar_max, species_name)
            w_from_min = get_w_from_v_res(dd_one, s1_q, vpar_min, species_name)
            w_from_zero = get_w_from_v_res(dd_one, s1_q, vpar_zero, species_name)

            curves_draw_vres_max = geom.Curve()
            curves_draw_vres_min = geom.Curve()
            curves_draw_vres_zero = geom.Curve()
            curves_draw_vres_max.add_curve(t_f * coef_norm_t, vpar_max)
            curves_draw_vres_min.add_curve(t_f * coef_norm_t, vpar_min)
            curves_draw_vres_zero.add_curve(t_f * coef_norm_t, vpar_zero)
            curves_draw_vres_max.color, curves_draw_vres_min.color, curves_draw_vres_zero.color = \
                'black', 'black', 'black'
            curves_draw_vres_max.style, curves_draw_vres_min.style, curves_draw_vres_zero.style = \
                ':', ':', ':'
            curves_draw_vres_max.width, curves_draw_vres_min.width, curves_draw_vres_zero.width = \
                4, 4, 4

            label_deltaf = species_name + ':\ \delta f(t, v_{res})'
            label_ddeltaf = species_name + ':\ \partial \delta f(t, v_{res})/\partial v'

        # - frequency normalization (coefficient) -
        coef_norm_w = None
        if sel_norm_w == 'khz':
            coef_norm_w = dd_one['wc'] / 1.e3
        if sel_norm_w == 'wc':
            coef_norm_w = 2 * np.pi
        if sel_norm_w == 'csa':
            coef_norm_w = 2 * np.pi * dd_one['wc'] / (dd_one['cs'] / dd_one['a0'])
        if sel_norm_w == 'csr':
            coef_norm_w = 2 * np.pi * dd_one['wc'] / (dd_one['cs'] / dd_one['R0'])

        # working frequency domain
        ids_w_work, w_norm_work, _ = mix.get_ids(w * coef_norm_w, w_norm_domain)
        ids_w_work = np.arange(ids_w_work[0], ids_w_work[-1] + 1)
        fft_work = res_fft[:, ids_w_work]

        # x domain of plotting
        x_plot, ids_x_plot = mix.get_array_oo(oo, x_work, 'x')
        ids_x_plot = np.arange(ids_x_plot[0], ids_x_plot[-1] + 1)
        fft_plot = fft_work[ids_x_plot, :]

        # plot x-evolution of the signal
        curves_t.new()\
            .XS(x * coef_norm_t)\
            .YS(data)\
            .leg(leg + ':\ original').sty('-')
        curves_t.new() \
            .XS(x * coef_norm_t) \
            .YS(data_filt) \
            .leg(leg + ':\ filtered').sty(':')
        curves_t.new() \
            .XS(x * coef_norm_t) \
            .YS(data_work) \
            .leg(leg + ':\ original - filtered').sty('-')

        # plot x-evolution of the signal FFT
        curves_w.new()\
            .XS(x_plot * coef_norm_t)\
            .YS(w_norm_work)\
            .ZS(fft_plot).lev(60)

        if flag_fv:
            curves_w.new()\
                .XS(t_f * coef_norm_t)\
                .YS(w_from_max * coef_norm_w / (2 * np.pi))
            curves_w.new() \
                .XS(t_f * coef_norm_t) \
                .YS(w_from_min * coef_norm_w / (2 * np.pi))

            curves_f = crv.Curves() \
                .xlab(labx + line_norm_t) \
                .ylab('v_{\parallel}') \
                .tit(label_deltaf)
            curves_f.new() \
                .XS(t_f * coef_norm_t) \
                .YS(vpar_f)\
                .ZS(deltaf)
            curves_f.newg(curves_draw_vres_max)

        # if len(curves_t.list_curves) is not 0:
        #     cpr.plot_curves(curves_t)
        cpr.plot_curves_3d(curves_w)
        if flag_fv:
            cpr.plot_curves_3d(curves_f)


# NEW: 1d plot: localisation:
def plot_vars_1d_localisation(oo):
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

    vvars  = choose_vars(oo_use)
    n_vars = len(vvars)

    # additional data:
    name_axis_loc = oo.get('name_axis_loc', None)

    labx       = oo.get('labx', None)  # x-label
    laby       = oo.get('laby', None)  # y-label
    leg_pos    = oo.get('leg_pos', None)
    tit_plot   = oo.get('tit_plot', None)  # title
    stys       = oo.get('stys', None)
    ylim = oo.get('ylim', None)

    flag_norm     = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)
    sel_norm      = oo.get('sel_norm', 'orig')

    dds = oo.get('dds', None)

    # normalization (first stage):
    coef_x_norm, line_x_norm = None, None
    if sel_norm == 'orig':
        coef_x_norm = 1
        line_x_norm = ''
    if sel_norm == 't-mili-seconds':
        line_x_norm = '(ms)'

    # plotting:
    curves = crv.Curves().xlab(labx + line_x_norm).ylab(laby).tit(tit_plot)\
        .leg_pos(leg_pos).ylim(ylim)
    curves.flag_norm     = flag_norm
    curves.flag_semilogy = flag_semilogy
    count_line = -1
    for ivar in range(n_vars):
        count_line += 1

        vvar = vvars[ivar]
        data_2d = vvar['data']
        leg_one = vvar['legs']

        x1, ids_x1 = mix.get_array_oo(oo, vvar[vvar['x1']], vvar['x1'])
        x2, ids_x2 = mix.get_array_oo(oo, vvar[vvar['x2']], vvar['x2'])
        data_2d = mix.get_slice(data_2d, ids_x1, ids_x2)

        # find localisation
        x_loc, dir_loc, x_evo, dir_evo = None, None, None, None
        # func_data = None
        if name_axis_loc == vvar['x1']:
            x_loc, dir_loc = x1, 0
            x_evo, dir_evo = x2, 1
            func_data = lambda f, id_ev, id_loc: f[id_loc, id_ev]
        if name_axis_loc == vvar['x2']:
            x_loc, dir_loc = x2, 1
            x_evo, dir_evo = x1, 0
            func_data = lambda f, id_ev, id_loc: f[id_ev, id_loc]

        # find absolute maximums at every evolution-moment
        ids_max = np.argmax(np.abs(data_2d), axis=dir_loc)

        # find locations of these maximums at every evolution-moment
        locs = np.zeros(len(x_evo))
        for id_evol in range(len(x_evo)):
            locs[id_evol] = x_loc[ids_max[id_evol]]

        # maximum values at every evolution-moment
        data_max = np.zeros(len(x_evo))
        data_2d_abs = np.abs(data_2d)
        for id_evol in range(len(x_evo)):
            data_max[id_evol] = func_data(data_2d_abs, id_evol, ids_max[id_evol])

        # find peaks among maximums:
        ids_max_peaks, _ = scipy.signal.find_peaks(data_max)
        x_evo_peaks = x_evo[ids_max_peaks]
        data_max_peaks = data_max[ids_max_peaks]

        # find space locations of these peaks:
        locs_peaks = locs[ids_max_peaks]

        # normalization (second stage):
        dd_one = dds[ivar]
        if sel_norm == 't-mili-seconds':
            coef_x_norm = 1. / dd_one['wc'] * 1e3

        # plotting
        sty_current = 'o'
        if stys is not None:
            if count_line < len(stys):
                sty_current = stys[count_line]

        # curves.new().XS(x_evo * coef_x_norm).YS(locs).leg(leg_one).sty(sty_current)
        curves.new().XS(x_evo_peaks * coef_x_norm).YS(locs_peaks) \
            .leg(leg_one).sty(sty_current)

        # curves.new().XS(x_evo * coef_x_norm).YS(data_max) \
        #     .leg(leg_one).sty(sty_current)
        # curves.new().XS(x_evo_peaks * coef_x_norm).YS(data_max_peaks) \
        #     .leg(leg_one).sty('o')

    if len(curves.list_curves) is not 0:
        cpr.plot_curves(curves)


# Get some properties of a signal:
def get_value_signal(oo):
    vvars = choose_vars(oo)
    n_vars = len(vvars)
    oo_operations = oo.get('oo_operations', None)  # for every var (not averaging)
    sel_norm = oo.get('sel_norm', 'orig')
    dds = oo.get('dds', None)

    def perform_operation(data, x, oo_one_operation):
        sel_operation = oo_one_operation[0]
        x_domain      = oo_one_operation[1]

        ids_x_op, x_op, line_x_op = mix.get_ids(x, x_domain)
        ids_x_op = range(ids_x_op[0], ids_x_op[-1] + 1)
        x_work = x[ids_x_op]

        value_res, var_res = np.nan, np.nan
        if sel_operation == 'integration':
            value_res = np.trapz(data[ids_x_op], x=x_work)

            curves = crv.Curves()
            curves.new().XS(x).YS(data)
            curves.new().XS(x_work).YS(data[ids_x_op]).sty(':')
            # curves.flag_semilogy = True
            cpr.plot_curves(curves)
        if sel_operation == 'mean':
            value_res = np.mean(data[ids_x_op])
        if sel_operation == 'max':
            value_res = np.max(data[ids_x_op])
        if sel_operation == 'sum':
            value_res = np.sum(data[ids_x_op])
        if sel_operation == 'abs-max':
            data_work = data[ids_x_op]
            id_max = np.argmax(np.abs(data_work))
            value_res = np.abs(data_work[id_max])

            curves = crv.Curves()
            curves.new().XS(x_work).YS(data_work)
            curves.new().XS(x_work[id_max]).YS(value_res).sty('o')
            curves.flag_semilogy = True
            cpr.plot_curves(curves)
        if sel_operation == 'abs-mean':
            value_res = np.mean(np.abs(data_work))
        if sel_operation == 'abs-peaks-mean':
            data_work = data[ids_x_op]
            ids_peaks, _ = scipy.signal.find_peaks(np.abs(data_work))
            data_abs_peaks = np.abs(data_work[ids_peaks])
            value_res = np.mean(data_abs_peaks)
            var_res   = np.var(data_abs_peaks)

            curves = crv.Curves()
            curves.new().XS(x_work).YS(data_work)
            curves.new().XS(x_work[ids_peaks]).YS(data_work[ids_peaks]).sty('o')
            curves.new().XS(x_work[ids_peaks]).YS(data_abs_peaks).sty('s')
            curves.flag_semilogy = True
            cpr.plot_curves(curves)
        if sel_operation == 'abs-peaks-peaks-mean':
            data_work = data[ids_x_op]
            ids_peaks, _ = scipy.signal.find_peaks(np.abs(data_work))
            data_abs_peaks = np.abs(data_work[ids_peaks])

            ids_peaks_peaks, _ = scipy.signal.find_peaks(data_abs_peaks)
            data_abs_peaks_peaks = data_abs_peaks[ids_peaks_peaks]
            value_res = np.mean(data_abs_peaks_peaks)
            var_res = np.var(data_abs_peaks_peaks)

            curves = crv.Curves()
            curves.new().XS(x_work).YS(data_work)
            curves.new()\
                .XS(x_work[ids_peaks])\
                .YS(data_abs_peaks)\
                .sty('o').leg('peaks')
            curves.new()\
                .XS(x_work[ids_peaks][ids_peaks_peaks])\
                .YS(data_abs_peaks_peaks)\
                .sty('s').leg('peaks-peaks')
            curves.flag_semilogy = True
            cpr.plot_curves(curves)
        if sel_operation == 'filter-mean':
            oo_filt = oo.get('oo_filt', None)
            data_work = data[ids_x_op]
            data_filt, x_filt = global_filter_signal(x_work, data_work, oo_filt)
            data_filt = np.interp(x_work, x_filt, data_filt)
            value_res = np.mean(data_filt)
            var_res   = np.var(data_filt)

            curves = crv.Curves()
            curves.new().XS(x_work).YS(data_work).leg('init')
            curves.new().XS(x_work).YS(data_filt).leg('filtered').sty(':')
            curves.flag_semilogy = True
            cpr.plot_curves(curves)

        return value_res, var_res, line_x_op

    count_line = -1
    props, desc, vars_res = [], [], []
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        datas = vvar['data']
        x     = vvar['x']
        legs  = vvar['legs']
        ns_av = len(datas)
        oo_one_operation = oo_operations[ivar]
        dd = dds[ivar]

        norm_coef = 1
        norm_line = ''
        if sel_norm == 'energy-transfer':
            norm_line = 'kW/m3'
            T_speak_J  = dd['T_speak']
            T_speak_eV = T_speak_J / constants.elementary_charge
            norm_coef = dd['cs'] * T_speak_eV / ( dd['a0'])
            norm_coef *= 1.e-3
        if sel_norm == 'energy-Jpm':
            norm_line = 'J/m3'
            T_speak_J  = dd['T_speak']
            T_speak_eV = T_speak_J / constants.elementary_charge
            norm_coef = dd['cs'] * T_speak_eV / (dd['a0'] * dd['wc'])

        for is_av in range(ns_av):
            count_line += 1
            data = datas[is_av]
            prop_res, var_res, line_oper = perform_operation(
                data, x, oo_one_operation)
            props.append(prop_res * norm_coef)
            vars_res.append(var_res * norm_coef)
            desc.append(legs[is_av] + ':\ ' + line_oper + ':\ ' + norm_line + ':\ ')

    for id_res in range(len(props)):
        print(desc[id_res] + ':\ ' + '{:0.3e} +- {:0.3e}'
              .format(props[id_res], vars_res[id_res]))


def operation_along_x_old(oo):
    def perform_operation(data, x, oo_one_operation, name_x):
        def get_x_data_interval(x_domain_loc, x_loc, data_loc):
            ids_x, _, _ = mix.get_ids(x_loc, x_domain_loc)
            ids_x = range(ids_x[0], ids_x[-1] + 1)
            return  data_loc[ids_x], x_loc[ids_x]

        # operation and working domain
        sel_operation = oo_one_operation[0]
        x_domain = oo_one_operation[1]
        data_work, x_work = get_x_data_interval(x_domain, x, data)

        data_res, x_res = [], []
        line_op1, line_op2 = '', ''
        if sel_operation == 'integration-accumulation':
            line_op1, line_op2 = '\\int_0^{' + name_x + '}', 'd' + name_x
            width_domain = oo_one_operation[2]

            x_start = x_work[0]
            x_end = x_start + width_domain
            while x_end <= x_work[-1]:
                data_curr, x_curr = get_x_data_interval([x_start, x_end], x_work, data_work)
                value_curr = np.trapz(data_curr, x=x_curr)

                data_res.append(value_curr)
                x_res.append(x_end)

                x_end += width_domain
            x_res = np.array(x_res)
            data_res = np.array(data_res)
        if sel_operation == 'integration-step':
            line_op1, line_op2 = '\\int_{' + name_x + '1}^{' + name_x + '2}', 'd' + name_x
            width_domain = oo_one_operation[2]
            step_domain  = oo_one_operation[3]

            x_start = x_work[0]
            x_end   = x_start + width_domain
            while x_end <= x_work[-1]:
                data_curr, x_curr = get_x_data_interval([x_start, x_end], x_work, data_work)
                value_curr = np.trapz(data_curr, x=x_curr)

                x_point = (x_end + x_start) / 2.
                data_res.append(value_curr)
                x_res.append(x_point)

                x_start += step_domain
                x_end    = x_start + width_domain
            x_res = np.array(x_res)
            data_res = np.array(data_res)
        if sel_operation == 'none':
            data_res = np.array(data_work)
            x_res    = np.array(x_work)

        return data_res, x_res, line_op1, line_op2

    # --- INPUT parameters ---
    vvars = choose_vars(oo)
    n_vars = len(vvars)
    oo_operations = oo.get('oo_operations', None)  # for every var (not averaging)
    dds = oo.get('dds', None)

    sel_norm_x = oo.get('sel_norm_x', 'orig')
    sel_norm_ys = oo.get('sel_norm_ys', 'orig')
    labx       = oo.get('labx', None)  # x-label
    laby       = oo.get('laby', None)  # y-label
    leg_pos    = oo.get('leg_pos', None)
    tit_plot   = oo.get('tit_plot', None)  # title
    stys       = oo.get('stys', None)
    cols       = oo.get('cols', None)
    ylim       = oo.get('ylim', None)

    flag_norm     = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)
    opt_legs      = oo.get('legs', [])
    flag_abs = oo.get('flag_abs', False)

    # normalization (first stage):
    coef_x_norm, line_x_norm = 1, ''
    coef_y_norm, line_y_norm = 1, ''
    if sel_norm_x == 't-mili-seconds':
        line_x_norm = '(ms)'

    if len(sel_norm_ys) == 1:
        if sel_norm_ys[0] == 'energy-transfer':
            line_y_norm = '\ [kW/m^3]'
        if sel_norm_ys[0] == 'energy-Jpm':
            line_y_norm = '\ [J/m^3]'

    # --- Calculation ---
    data_vars, x_vars, line_start_vars, line_end_vars = [], [], [], []
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        data = vvar['data'][0]
        x     = vvar['x']
        name_x = vvar['name_x']
        oo_one_operation = oo_operations[ivar]

        data_one, x_one, line_start_one, line_end_one = \
            perform_operation(data, x, oo_one_operation, name_x)
        data_vars.append(data_one)
        x_vars.append(x_one)
        line_start_vars.append(line_start_one)
        line_end_vars.append(line_end_one)

    # --- Plotting ---
    if labx is not None: labx += line_x_norm
    if laby is not None: laby += line_y_norm
    curves = crv.Curves().xlab(labx).ylab(laby).tit(tit_plot)\
        .leg_pos(leg_pos).ylim(ylim)
    curves.flag_norm     = flag_norm
    curves.flag_semilogy = flag_semilogy
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        legs = vvar['legs']

        # normalization (second stage):
        dd_one = dds[ivar]
        if sel_norm_x == 't-mili-seconds':
            coef_x_norm = 1. / dd_one['wc'] * 1e3

        line_y_norm = ''
        if sel_norm_ys[ivar] == 'energy-transfer':
            line_y_norm = '\ [kW/m^3]'
            T_speak_eV   = dd_one['T_speak'] / constants.elementary_charge
            coef_y_norm  = dd_one['cs'] * T_speak_eV / (dd_one['a0'])
            coef_y_norm *= 1.e-3
        if sel_norm_ys[ivar] == 'energy-Jpm':
            line_y_norm = '\ [J/m^3]'
            T_speak_eV  = dd_one['T_speak'] / constants.elementary_charge
            coef_y_norm = dd_one['cs'] * T_speak_eV / (dd_one['a0'] * dd_one['wc'])

        # get data
        x    = x_vars[ivar]    * coef_x_norm
        data = data_vars[ivar] * coef_y_norm

        # styles
        sty_current = '-'
        if stys is not None:
            if ivar < len(stys):
                sty_current = stys[ivar]

        # colors
        color_current = None
        if cols is not None:
            if ivar < len(cols):
                color_current = cols[ivar]

        # legends
        one_leg = legs[0] + line_y_norm
        if ivar < len(opt_legs):
            if opt_legs[ivar] is not None:
                one_leg = opt_legs[ivar]

        curves.new()\
            .XS(x).YS(data)\
            .leg(one_leg).sty(sty_current).col(color_current)
        if flag_abs:
            curves.list_curves[-1].YS(np.abs(data))
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

    # number of dds has to be equal to n_points * n_x_fft_domains

    # choose variables and their averaging:
    vvars = get_fft_1d(oo)
    n_vars = len(vvars['datas'])

    dds = oo.get('dds', None)

    # frequency domain for plotting:
    w_domain = oo.get('w_domain', None)

    # - frequency normalization (notation) -
    sel_norm = oo.get('sel_norm', 'wc')

    line_w = None
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
    laby = oo.get('laby', None)  # y-label
    tit_plot  = oo.get('tit_plot', None)  # title

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

    line_w = None
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

    tit_plot = oo.get('tit_plot', None)  # title
    labx = oo.get('labx', 't[wci^{-1}]')
    laby = oo.get('laby', None)  # y-label
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


# NEW: calculation of the frequency and dynamic rate at one radial point:
def calc_wg(oo_var, oo_wg, oo_plot):
    # -------------------------------------------------------------------------------
    # -> oo_var - dictionary to choose a variable
    #   (input dict. for the function choose_vars(...))
    # -------------------------------------------------------------------------------
    # -> oo_wg - dictionary with parameters to calculate frequency and dynamic rate:
    # 't_work' - work time domain
    # 'sel_wg' - line, name of a method to calculate frequency and rate:
    #   'wg-adv', 'w-adv', 'g-adv', 'wg-est', 'w-est', 'g-est'
    # 'flag_two_stages' = True:
    #       Firstly, one calculates gamma,
    #       then create a signal = initial_signal * exp(-gamma*t) and
    #       calculate the frequency
    #   False:
    #       calculate gamma and frequency from the initial_signal
    # 'sel_norm': 'wc', 'vt', 'khz':
    #       output normalization
    # ---
    # 'flag_stat' = True:
    #       calculate errorbars
    # 'n_samples' - integer:
    #       number of time interval variations
    # 'min_n_peaks' - integer:
    #       minimum number of peaks in one time interval
    # 'threshold_w' - float:
    #       relative difference between estimated (linear fitting) value of frequency
    #       (w_est) and
    #       frequency value found from NL fitting (w_adv),
    #       if |(w_adv - w_est)/w_est| <= threshold_w, then we are taking w_adv as
    #       a result frequency, otherwise we don't take any value
    # 'threshold_g' - float:
    #       the same as threshold_w, but for the damping rate
    # ---
    # - FILTERING -
    #  -> If 'flag_two_stages' = True, there are three stages of the filtering:
    #       global, for gamma, for frequency;
    #  -> Globally filtered signal is a starting signal for the calculation of the
    #       both gamma and frequency;
    #  -> After that, globally filtered signal can be filtered separately
    #       before the calculation of the gamma and before the calc. of the frequency
    #  -> If 'flag_two_stages' = False, there is only global filtering
    # 'filt_global' - dict. or [dict., dict., ...]:
    #       global filtering
    # 'filt_gamma' - dict. or [dict., dict., ...]:
    #       additional filtering of the globally filtered signal
    #       before the calculation of the gamma
    # 'filt_freq' - dict. or [dict., dict., ...]:
    #       additional filtering of the globally filtered signal
    #       before the calculation of the frequency
    #  -> For the description of these dictionaries, see the function ymath.filtering
    # -------------------------------------------------------------------------------
    # -> oo_plot - dictionary for plotting:
    # 't_plot' - domain of plotting;
    # 'flag_norm' = True: normalized plots;
    # 'flag_semilogy' = True: Y-axis in logarithmic scale;
    # 'flag_plot_print' = True: plot results and print values on screen
    # -------------------------------------------------------------------------------

    # - None-filter -
    non_filt = {'sel_filt': None}
    out_res = {}

    # - FUNCTION: filtering at one stage -
    def one_stage_filtering(x, y, oo_filt_loc):
        oo_filts_res = []
        if oo_filt_loc is None or len(oo_filt_loc) is 0:
            oo_filts_res.append(non_filt)
        elif isinstance(oo_filt_loc, dict):
            oo_filts_res.append(oo_filt_loc)
        else:
            oo_filts_res = oo_filt_loc  # array of filters

        # apply one by one the filters in the array :
        dict_loc = {'x': np.array(x), 'filt': np.array(y)}
        count_filt = -1
        for one_filt in oo_filts_res:
            count_filt += 1
            dict_loc = ymath.filtering(
                dict_loc['x'], dict_loc['filt'], one_filt
            )
        return dict_loc

    # - FUNCTION: get result -
    def give_res(dict_wg_loc, name_method, name_value, coef_norm):
        value_res, err_value = None, None
        if dict_wg_loc[name_method] is not None:
            value_res = dict_wg_loc[name_method][name_value]
            err_value = dict_wg_loc[name_method][name_value + '_err']

        # normalization:
        if value_res is not None:
            value_res *= coef_norm
        if err_value is not None:
            err_value *= coef_norm

        line_value = 'None'
        if value_res is not None:
            line_value = '{:0.3e}'.format(value_res)
            if err_value is not None:
                line_value += ' +- {:0.3e}'.format(err_value)

        return value_res, line_value

    # - FUNCTION: different options of w and g measurement -
    def find_wg(t, y, sel_wg, flag_print=True):
        if sel_wg == 'wg-adv':
            res_dict = ymath.advanced_wg(t, y, flag_print=flag_print)
        elif sel_wg == 'w-adv':
            res_dict = ymath.advanced_w(t, y, flag_print=flag_print)
        elif sel_wg == 'g-adv':
            res_dict = ymath.advanced_g(t, y, flag_print=flag_print)
        elif sel_wg == 'wg-est':
            res_dict = {
                'est': ymath.estimate_wg(t, y),
                'adv': None
            }
        elif sel_wg == 'w-est':
            res_dict = {
                'est': ymath.estimate_w(t, y),
                'adv': None
            }
        elif sel_wg == 'g-est':
            res_dict = {
                'est': ymath.estimate_g(t, y),
                'adv': None
            }
        else:
            print('Error:\ Wrong wg-selector: check oo_wg[''sel_wg'']')
            sys.exit(1)

        return res_dict

    # - project structure -
    dd = oo_var.get('dds', None)[0]

    # - choose a variable -
    dict_var = choose_vars(oo_var)[0]

    # - initial data -
    data_init = dict_var['data'][0]
    t_init = dict_var['x']
    dict_fft_init = ymath.filtering(t_init, data_init, non_filt)

    # --- Frequency/rate calculation PARAMETERS ---
    t_work = oo_wg.get('t_work', [])
    sel_wg = oo_wg.get('sel_wg', 'wg-adv')

    if len(t_work) == 0 or t_work is None:
        t_work = t_init
    flag_two_stages = oo_wg.get('flag_two_stages', False)
    if sel_wg == 'w-adv' or sel_wg == 'w-est' or\
            sel_wg == 'g-adv' or sel_wg == 'g-est':
        flag_two_stages = False

    oo_filt_global = oo_wg.get('filt_global', None)
    oo_filt_gamma = oo_wg.get('filt_gamma', None)
    oo_filt_freq = oo_wg.get('filt_freq', None)

    flag_stat   = oo_wg.get('flag_stat', False)
    n_samples   = oo_wg.get('n_samples', None)
    min_n_peaks = oo_wg.get('min_n_peaks', None)
    threshold_w = oo_wg.get('threshold_w', 0.1)
    threshold_g = oo_wg.get('threshold_g', 0.2)
    n_bins = oo_wg.get('n_bins', 40)

    flag_print = False
    if n_samples <= 10:
        flag_print = True

    sel_norm = oo_wg.get('sel_norm', 'wc')
    coef_norm_global_w, coef_norm_global_g, line_norm_w, line_norm_g = None, None, '', ''
    if sel_norm == 'wc':
        line_norm_w = line_norm_g = '[wci]'
        coef_norm_global_w = coef_norm_global_g = 1
    if sel_norm == 'vt':
        line_norm_w = line_norm_g = '[sqrt(2)*vt/R0]'
        coef_norm_global_w = coef_norm_global_g = \
            dd['wc'] / (np.sqrt(2) * dd['vt'] / dd['R0'])
    if sel_norm == 'khz':
        line_norm_w = '(kHz)'
        line_norm_g = '(10^3 s)'
        coef_norm_global_w = dd['wc'] / (1e3 * 2*np.pi)
        coef_norm_global_g = dd['wc'] / 1e3

    # --- GLOBAL FILTERING ---
    dict_global = one_stage_filtering(t_init, data_init, oo_filt_global)
    data_global = dict_global['filt']
    t_global    = dict_global['x']

    # --- WORK DOMAIN ---
    ids_work, t_work, _ = mix.get_ids(t_global, t_work)
    data_work = data_global[ids_work]
    dict_fft_work_global    = ymath.filtering(t_work, data_work, non_filt)
    ids_peaks_work, _       = scipy.signal.find_peaks(np.abs(data_work))
    t_peaks_work = t_work[ids_peaks_work]

    # - PLOTTING: GLOBAL FILTERING -

    # parameters of plotting
    flag_plot_print = oo_plot.get('flag_plot_print', True)

    if flag_plot_print:
        t_plot = oo_plot.get('t_plot', [])
        if len(t_plot) == 0:
            t_plot = t_init

        ids_plot_init, _, _   = mix.get_ids(t_init,   t_plot)
        ids_plot_global, _, _ = mix.get_ids(t_global, t_plot)
        ids_plot_work, _, _   = mix.get_ids(t_work,   t_plot)

        flag_norm     = oo_plot.get('flag_norm', False)
        flag_semilogy = oo_plot.get('flag_semilogy', False)

        leg_data = dict_var['legs'][0]

        # work domain
        area_work = geom.Fill()
        area_work.xs = [t_work[0], t_work[-1], t_work[-1], t_work[0]]
        area_work.ys = ['limb', 'limb', 'limu', 'limu']
        area_work.color = 'grey'
        area_work.alpha = 0.3

        # plotting time evolution
        curves = crv.Curves() \
            .xlab('t[\omega_{ci}^{-1}]') \
            .tit(dict_var['tit']) \
            .xlim([t_plot[0], t_plot[-1]])
        curves.flag_norm = flag_norm
        curves.flag_semilogy = flag_semilogy
        if flag_norm:
            curves.ylab('norm.')
        curves.newg(area_work)
        curves.new() \
            .XS(t_init[ids_plot_init]) \
            .YS(data_init[ids_plot_init]) \
            .leg(leg_data)
        curves.new() \
            .XS(t_global[ids_plot_global]) \
            .YS(data_global[ids_plot_global]) \
            .leg(leg_data + ': globally\ filtered').sty(':')
        curves.new().norm_to(data_global[ids_plot_init])\
            .XS(t_peaks_work)\
            .YS(data_work[ids_peaks_work])\
            .leg(leg_data + ':\ peaks').sty('o')
        cpr.plot_curves(curves)

        # plotting FFT
        curves = crv.Curves() \
            .xlab('\omega[\omega_{ci}]') \
            .tit(leg_data + ':\ FFT')
        curves.flag_norm = flag_norm
        if flag_norm:
            curves.ylab('norm.')
        curves.new() \
            .XS(dict_fft_init['w2']) \
            .YS(dict_fft_init['fft_init_2']) \
            .leg('FFT:\ initial')
        curves.new() \
            .XS(dict_global['w2']) \
            .YS(dict_global['fft_filt_2']) \
            .leg('FFT:\ globally\ filtered:\ whole\ time\ domain').sty(':')
        curves.new() \
            .XS(dict_fft_work_global['w2']) \
            .YS(dict_fft_work_global['fft_init_2']) \
            .leg('FFT:\ globally\ filtered:\ work\ time\ domain').sty(':')
        cpr.plot_curves(curves)

    # --- NAIVE CALCULATION ---
    if not flag_two_stages:
        dict_wg = find_wg(t_work, data_work, sel_wg=sel_wg)

        # - results -
        w_est, line_w_est = give_res(dict_wg, 'est', 'w', coef_norm_global_w)
        g_est, line_g_est = give_res(dict_wg, 'est', 'g', coef_norm_global_g)
        w_adv, line_w_adv = give_res(dict_wg, 'adv', 'w', coef_norm_global_w)
        g_adv, line_g_adv = give_res(dict_wg, 'adv', 'g', coef_norm_global_g)

        if flag_plot_print:
            # - plotting -
            curves = crv.Curves() \
                .xlab('t[\omega_{ci}^{-1}]') \
                .tit('Freq./Gamma\ calculation')
            curves.flag_norm = flag_norm
            curves.flag_semilogy = flag_semilogy
            if flag_norm:
                curves.ylab('norm.')
            curves.new() \
                .XS(t_work) \
                .YS(data_work) \
                .leg(leg_data)
            curves.new().norm_to(data_work) \
                .XS(t_work[dict_wg['est']['ids_peaks']]) \
                .YS(data_work[dict_wg['est']['ids_peaks']]) \
                .leg('peaks').sty('o')
            curves.new().norm_to(data_work)  \
                .XS(dict_wg['est']['x_fit']) \
                .YS(dict_wg['est']['y_fit']) \
                .leg('LIN.\ REGRESSION').sty(':')
            if dict_wg['adv'] is not None:
                curves.new() \
                    .XS(dict_wg['adv']['x_fit']) \
                    .YS(dict_wg['adv']['y_fit']) \
                    .leg('NL\ FITTING').sty(':')
            cpr.plot_curves(curves)

            line_res = '--- NAIVE CALCULATION ---\n'
            line_res += '--- ESTIMATION ---\n'
            line_res += 'w' + line_norm_w + ' = ' + line_w_est + '\n'
            line_res += 'g' + line_norm_g + ' = ' + line_g_est + '\n'
            line_res += '--- NL FITTING ---\n'
            line_res += 'w' + line_norm_w + ' = ' + line_w_adv + '\n'
            line_res += 'g' + line_norm_g + ' = ' + line_g_adv
            print(line_res)
    else:
        # - FIND GAMMA -
        dict_gamma_filt = one_stage_filtering(t_work, data_work, oo_filt_gamma)
        data_gamma = dict_gamma_filt['filt']
        t_gamma    = dict_gamma_filt['x']
        dict_gamma = find_wg(t_gamma, data_gamma, sel_wg)

        w_est_prel, line_w_est_prel = give_res(dict_gamma, 'est', 'w', coef_norm_global_w)
        g_est, line_g_est           = give_res(dict_gamma, 'est', 'g', coef_norm_global_g)
        w_adv_prel, line_w_adv_prel = give_res(dict_gamma, 'adv', 'w', coef_norm_global_w)
        g_adv, line_g_adv           = give_res(dict_gamma, 'adv', 'g', coef_norm_global_g)

        # - FIND FREQUENCY -
        g_mult = g_est
        if g_adv is not None:
            g_mult = g_adv
        g_mult /= coef_norm_global_g
        data_work_exp = data_work * np.exp(- g_mult * t_work)
        dict_freq_filt = one_stage_filtering(
            t_work,
            data_work_exp,
            oo_filt_freq
        )
        data_freq = dict_freq_filt['filt']
        t_freq    = dict_freq_filt['x']
        dict_freq = find_wg(t_freq, data_freq, sel_wg)

        w_est, line_w_est           = give_res(dict_freq, 'est', 'w', coef_norm_global_w)
        g_est_zero, line_g_est_zero = give_res(dict_freq, 'est', 'g', coef_norm_global_g)
        w_adv, line_w_adv           = give_res(dict_freq, 'adv', 'w', coef_norm_global_w)
        g_adv_zero, line_g_adv_zero = give_res(dict_freq, 'adv', 'g', coef_norm_global_g)

        if flag_plot_print:
            # plotting: gamma: filtering
            curves = crv.Curves() \
                .xlab('t[\omega_{ci}^{-1}]') \
                .tit('Gamma\ calculation:\ filtering')
            curves.flag_norm        = flag_norm
            curves.flag_semilogy    = flag_semilogy
            if flag_norm:
                curves.ylab('norm.')
            curves.new() \
                .XS(t_work) \
                .YS(data_work) \
                .leg(leg_data)
            curves.new() \
                .XS(t_gamma) \
                .YS(data_gamma) \
                .leg(leg_data + ':\ filtered').sty(':')
            cpr.plot_curves(curves)

            # plotting: gamma: fft
            curves = crv.Curves() \
                .xlab('\omega[\omega_{ci}]') \
                .tit('Gamma\ calculation:\ FFT')
            curves.flag_norm = flag_norm
            if flag_norm:
                curves.ylab('norm.')
            curves.new() \
                .XS(dict_gamma_filt['w2']) \
                .YS(dict_gamma_filt['fft_init_2']) \
                .leg('FFT:\ global\ filtering')
            curves.new() \
                .XS(dict_gamma_filt['w2']) \
                .YS(dict_gamma_filt['fft_filt_2']) \
                .leg('FFT:\ gamma\ filtering').sty(':')
            cpr.plot_curves(curves)

            # plotting: gamma: fitting
            curves = crv.Curves() \
                .xlab('t[\omega_{ci}^{-1}]') \
                .tit('Gamma\ calculation:\ fitting')
            curves.flag_norm = flag_norm
            curves.flag_semilogy = flag_semilogy
            if flag_norm:
                curves.ylab('norm.')
            curves.new() \
                .XS(t_gamma) \
                .YS(data_gamma) \
                .leg(leg_data + ':\ filtered')
            curves.new().norm_to(data_work) \
                .XS(t_work[dict_gamma['est']['ids_peaks']]) \
                .YS(data_work[dict_gamma['est']['ids_peaks']]) \
                .leg('peaks').sty('o')
            curves.new().norm_to(data_work) \
                .XS(dict_gamma['est']['x_fit']) \
                .YS(dict_gamma['est']['y_fit']) \
                .leg('LIN.\ REGRESSION').sty(':')
            if dict_gamma['adv'] is not None:
                curves.new() \
                    .XS(dict_gamma['adv']['x_fit']) \
                    .YS(dict_gamma['adv']['y_fit']) \
                    .leg('NL\ FITTING').sty(':')
            cpr.plot_curves(curves)

            # plotting: frequency: filtering
            curves = crv.Curves() \
                .xlab('t[\omega_{ci}^{-1}]') \
                .tit('Freq.\ calc.:\ filtering')
            curves.flag_norm = flag_norm
            curves.flag_semilogy = flag_semilogy
            if flag_norm:
                curves.ylab('norm.')
            curves.new() \
                .XS(t_work) \
                .YS(data_work) \
                .leg(leg_data)
            curves.new() \
                .XS(t_freq) \
                .YS(data_freq) \
                .leg(leg_data + ':\ *\exp(-g*t),\ filtered').sty(':')
            cpr.plot_curves(curves)

            # plotting: frequency: fft
            curves = crv.Curves() \
                .xlab('\omega[\omega_{ci}]') \
                .tit('Freq.\ calc.:\ FFT')
            curves.flag_norm = flag_norm
            if flag_norm:
                curves.ylab('norm.')
            curves.new() \
                .XS(dict_freq_filt['w2']) \
                .YS(dict_freq_filt['fft_init_2']) \
                .leg('FFT:\ global\ filtering')
            curves.new() \
                .XS(dict_freq_filt['w2']) \
                .YS(dict_freq_filt['fft_filt_2']) \
                .leg('FFT:\ freq.\ filtering').sty(':')
            cpr.plot_curves(curves)

            # plotting: frequency: fitting
            curves = crv.Curves() \
                .xlab('t[\omega_{ci}^{-1}]') \
                .tit('Freq.\ calc.:\ fitting')
            curves.flag_norm = flag_norm
            curves.flag_semilogy = flag_semilogy
            if flag_norm:
                curves.ylab('norm.')
            curves.new() \
                .XS(t_freq) \
                .YS(data_freq) \
                .leg('filtered')
            curves.new().norm_to(data_work_exp) \
                .XS(t_work[dict_freq['est']['ids_peaks']]) \
                .YS(data_work_exp[dict_freq['est']['ids_peaks']]) \
                .leg('peaks').sty('o')
            curves.new() \
                .XS(dict_freq['est']['x_fit']) \
                .YS(dict_freq['est']['y_fit']) \
                .leg('LIN.\ REGRESSION').sty(':')
            if dict_freq['adv'] is not None:
                curves.new() \
                    .XS(dict_freq['adv']['x_fit']) \
                    .YS(dict_freq['adv']['y_fit']) \
                    .leg('NL\ FITTING').sty(':')
            cpr.plot_curves(curves)

            # - print results -
            line_res = '--- NAIVE CALCULATION ---\n'
            line_res += '- GAMMA: ESTIMATION -\n'
            line_res += 'prel. w' + line_norm_w + ' = ' + line_w_est_prel + '\n'
            line_res += 'g' + line_norm_g + ' = ' + line_g_est + '\n'
            line_res += '- GAMMA: NL FITTING -\n'
            line_res += 'prel. w' + line_norm_w + ' = ' + line_w_adv_prel + '\n'
            line_res += 'g' + line_norm_g + ' = ' + line_g_adv + '\n'
            line_res += '- FREQUENCY: ESTIMATION -\n'
            line_res += 'w' + line_norm_w + ' = ' + line_w_est + '\n'
            line_res += '(g_real - g_num)' + line_norm_g + ' = ' + line_g_est_zero + '\n'
            line_res += '- FREQUENCY: NL FITTING -\n'
            line_res += 'w' + line_norm_w + ' = ' + line_w_adv + '\n'
            line_res += '(g_real - g_num)' + line_norm_g + ' = ' + line_g_adv_zero
            print(line_res)

    # - save results of naive calculation -
    naive_res = {
        'w_est': w_est, 'g_est': g_est,
        'w_adv': w_adv, 'g_adv': g_adv,
    }
    out_res.update({'naive': naive_res})

    # --- CALCULATION WITH STATISTICS ---
    t_intervals = []
    if flag_stat:
        oo_get_intervals = {
            'nsamples': n_samples,
            't_work': t_work,
            'min_n_periods': min_n_peaks,
            't_period': np.mean(np.diff(t_peaks_work))
        }
        dict_intervals = mix.get_t_intervals(oo_get_intervals, flag_plot_print)
        t_intervals = dict_intervals['t_intervals']

        # plot time intervals
        if flag_print and flag_plot_print:
            for one_t_interval in t_intervals:
                area_calc_chosen = geom.Fill()
                area_calc_chosen.xs = [
                    one_t_interval[0], one_t_interval[-1],
                    one_t_interval[-1], one_t_interval[0]
                ]
                area_calc_chosen.ys = ['limb', 'limb', 'limu', 'limu']
                area_calc_chosen.color = 'grey'
                area_calc_chosen.alpha = 0.6

                curves = crv.Curves() \
                    .xlab('t[\omega_{ci}^{-1}]') \
                    .tit('Chosen\ time\ intervals')
                curves.flag_norm = flag_norm
                curves.flag_semilogy = flag_semilogy
                if flag_norm:
                    curves.ylab('norm.')
                curves.newg(area_work)
                curves.newg(area_calc_chosen)
                curves.new() \
                    .XS(t_global) \
                    .YS(data_global) \
                    .leg('Globally\ filtered\ data')
                curves.new().norm_to(data_global) \
                    .XS(t_work[ids_peaks_work]) \
                    .YS(data_work[ids_peaks_work]) \
                    .leg('peaks').sty('o')
                cpr.plot_curves(curves)

    # - calculation of freq/rate at one stage -
    ws, gs = [], []
    if not flag_two_stages and flag_stat:
        for i_sample in range(n_samples):
            ids_one_t_interval, t_one_interval, _ = \
                mix.get_ids(t_work, t_intervals[i_sample])

            dict_wg = find_wg(
                t_one_interval,
                data_work[ids_one_t_interval],
                sel_wg,
                flag_print=False
            )
            w_est, line_w_est = give_res(dict_wg, 'est', 'w', coef_norm_global_w)
            g_est, line_g_est = give_res(dict_wg, 'est', 'g', coef_norm_global_g)
            w_adv, line_w_adv = give_res(dict_wg, 'adv', 'w', coef_norm_global_w)
            g_adv, line_g_adv = give_res(dict_wg, 'adv', 'g', coef_norm_global_g)

            if 'est' not in sel_wg:
                if w_adv is not None:
                    if np.abs((w_adv - w_est)/w_est) <= threshold_w:
                        ws.append(w_adv)
                if g_adv is not None:
                    if np.abs((g_adv - g_est) / g_est) <= threshold_g:
                        gs.append(g_adv)
            else:
                if w_est is not None:
                    ws.append(w_est)
                if g_est is not None:
                    gs.append(g_est)
        ws = np.array(ws)
        gs = np.array(gs)

    # - calculation of freq/rate at two stages -
    elif flag_two_stages and flag_stat:
        for i_sample in range(n_samples):

            # - FIND GAMMA -
            ids_one_t_interval, t_one_interval, _ = \
                mix.get_ids(t_gamma, t_intervals[i_sample])

            dict_gamma = find_wg(
                t_one_interval,
                data_gamma[ids_one_t_interval],
                sel_wg,
                flag_print=False
            )

            g_est, line_g_est = give_res(dict_gamma, 'est', 'g', coef_norm_global_g)
            g_adv, line_g_adv = give_res(dict_gamma, 'adv', 'g', coef_norm_global_g)
            if 'est' not in sel_wg:
                if g_adv is not None:
                    if np.abs((g_adv - g_est) / g_est) <= threshold_g:
                        gs.append(g_adv)
            else:
                if g_est is not None:
                    gs.append(g_est)

            # - FIND FREQUENCY -
            g_mult = g_est
            if g_adv is not None:
                g_mult = g_adv
            g_mult /= coef_norm_global_g
            dict_freq_filt = one_stage_filtering(
                t_work,  # filtering is performed in the work domain
                data_work * np.exp(- g_mult * t_work),
                oo_filt_freq
            )
            data_freq   = dict_freq_filt['filt']
            t_freq      = dict_freq_filt['x']

            ids_one_t_interval, t_one_interval, _ = \
                mix.get_ids(t_freq, t_intervals[i_sample])

            dict_freq = find_wg(
                t_one_interval,  # calculation is performed in one of the time domains
                data_freq[ids_one_t_interval],
                sel_wg,
                flag_print=False
            )

            w_est, line_w_est = give_res(dict_freq, 'est', 'w', coef_norm_global_w)
            w_adv, line_w_adv = give_res(dict_freq, 'adv', 'w', coef_norm_global_w)
            if 'est' not in sel_wg:
                if w_adv is not None:
                    if np.abs((w_adv - w_est)/w_est) <= threshold_w:
                        ws.append(w_adv)
            else:
                if w_est is not None:
                    ws.append(w_est)
        ws = np.array(ws)
        gs = np.array(gs)

    # - statistical results -
    stat_res = {
        'w': None, 'err_w': None,
        'g': None, 'err_g': None,
    }
    if flag_stat:
        # frequency histogram
        hist_w = np.histogram(ws, bins=n_bins)
        fit_mean_w, fit_sigma_w = stat_norm.fit(ws)
        f_data_w = stat_norm.pdf(hist_w[1], fit_mean_w, fit_sigma_w)

        # Rate histogram
        hist_g = np.histogram(gs, bins=n_bins)
        fit_mean_g, fit_sigma_g = stat_norm.fit(gs)
        f_data_g = stat_norm.pdf(hist_g[1], fit_mean_g, fit_sigma_g)

        err_w = ymath.COEF_ERR * fit_sigma_w
        err_g = ymath.COEF_ERR * fit_sigma_g

        # - plotting and printing statistical results -
        if flag_plot_print:
            curves = crv.Curves() \
                .xlab('\omega[\omega_{ci}]') \
                .ylab('a.u.') \
                # .tit('Histogram:\ Frequency')
            curves.flag_maxlocator = True
            curves.new() \
                .XS(n_bins) \
                .YS(ws) \
                .set_hist().alpha(1).leg('histogram')
            curves.new() \
                .XS(hist_w[1]) \
                .YS(f_data_w) \
                .col('red').leg('normal\ distribution')
            cpr.plot_curves(curves)

            curves = crv.Curves() \
                .xlab('\gamma[\omega_{ci}]') \
                .ylab('a.u.') \
                # .tit('Histogram:\ Damping Rate')
            curves.flag_maxlocator = True
            curves.new() \
                .XS(n_bins) \
                .YS(gs) \
                .set_hist().alpha(1).leg('histogram')
            curves.new() \
                .XS(hist_g[1]) \
                .YS(f_data_g) \
                .col('red').leg('normal\ distribution')
            cpr.plot_curves(curves)

            # Print results
            line_stat = '--- STATISTICS ---\n'
            line_stat += 'number of frequency samples = ' + '{:d}'.format(len(ws)) + '\n'
            line_stat += 'number of rate samples = ' + '{:d}'.format(len(gs)) + '\n'
            line_stat += 'w' + line_norm_w + ' = ' + '{:0.3e}'.format(fit_mean_w) \
                        + '+-' + '{:0.3e}'.format(err_w) + '\n'
            line_stat += 'g' + line_norm_g + ' = ' + '{:0.3e}'.format(fit_mean_g) \
                        + '+-' + '{:0.3e}'.format(err_g)
            print(line_stat)

        # - save results of naive calculation -
        stat_res = {
            'w': fit_mean_w, 'err_w': err_w,
            'g': fit_mean_g, 'err_g': err_g,
        }
    out_res.update({'stat': stat_res})

    return out_res


# NEW: calculation of the frequency and dynamic rate in a radial interval:
def calc_wg_s(oo_s, oo_var, oo_wg, oo_plot):
    # Calculate frequency (w) and damping/growth rate (g) at radial points in the some
    #   radial interval;
    # oo_s - dict. to describe radial domain where (w,g) will be calculated.
    #  'cond_res_w' - function:
    #    defines conditions to which result w has to satisfy
    #  'cond_res_g' - function:
    #    defines conditions to which result g has to satisfy
    smax = oo_s.get('smax', 0.0)
    smin = oo_s.get('smin', 1.0)
    ds = oo_s.get('ds', None)
    s  = oo_s.get('s_array', None)

    cond_res_w = oo_s.get('cond_res_w', None)
    cond_res_g = oo_s.get('cond_res_g', None)

    ns = int((smax - smin) / ds + 1)
    if s is None:
        s = np.linspace(smin, smax, ns)

    oo_plot.update({
        'flag_plot_print': False
    })

    s_res, w_est, g_est, w_adv, g_adv = [], [], [], [], []
    s_res_stat, w_stat, w_stat_err, g_stat, g_stat_err = [], [], [], [], []
    for s1 in s:
        print('--- s = {:0.3f} ---'.format(s1))

        oo_var.update({
            'avrs': [['ts', 'point-s', [s1]]],
        })
        try:
            out_res = calc_wg(oo_var, oo_wg, oo_plot)
        except:
            print('No results')
            print('---')
            continue

        # Naive results
        w1_est = out_res['naive']['w_est']
        g1_est = out_res['naive']['g_est']
        w1_adv = out_res['naive']['w_adv']
        g1_adv = out_res['naive']['g_adv']

        flag_check = True
        if w1_est is None or g1_est is None or \
                w1_adv is None or g1_adv is None:
            flag_check = False
        if cond_res_w is not None:
            if not cond_res_w(w1_est) or not cond_res_w(w1_adv):
                flag_check = False
        if cond_res_g is not None:
            if not cond_res_g(g1_est) or not cond_res_g(g1_adv):
                flag_check = False

        if flag_check:
            s_res.append(s1)
            w_est.append(w1_est)
            g_est.append(g1_est)
            w_adv.append(w1_adv)
            g_adv.append(g1_adv)
        else:
            print('No results for naive calculation')

        # Statistical results
        w1 = out_res['stat']['w']
        w1_err = out_res['stat']['err_w']
        g1 = out_res['stat']['g']
        g1_err = out_res['stat']['err_g']

        flag_check = True
        if w1 is None or g1 is None:
            flag_check = False
        if cond_res_w is not None:
            if not cond_res_w(w1):
                flag_check = False
        if cond_res_g is not None:
            if not cond_res_g(g1):
                flag_check = False

        if flag_check:
            s_res_stat.append(s1)
            w_stat.append(w1)
            g_stat.append(g1)
            w_stat_err.append(w1_err)
            g_stat_err.append(g1_err)
        else:
            print('No results for statistical calculation')
        print('---')

    curves = crv.Curves().xlab('s').ylab('\omega[\omega_{ci}]')\
        .tit('Naive\ frequency')
    curves.new().XS(s_res).YS(w_est).sty('o--')
    cpr.plot_curves(curves)

    curves = crv.Curves().xlab('s').ylab('\omega[\omega_{ci}]') \
        .tit('Statistical\ frequency')
    curves.new().XS(s_res_stat).YS(w_stat).sty('o--')\
        .set_errorbar(True, ys=w_stat_err)
    cpr.plot_curves(curves)

    return


# NEW: calculation of the w/g in different time intervals:
def calc_wg_t_simplified(oo_wg_t, oo_plot):
    # configurations
    oo_cong_names = oo_wg_t['conf_names']
    oo_ts   = oo_wg_t['oo_ts']
    oo_vars = oo_wg_t['oo_vars']
    oo_wgs  = oo_wg_t['oo_wgs']

    # number of configurations
    n_conf = len(oo_cong_names)

    # output properties
    tit_plot = oo_plot.get('tit_plot', None)
    stys_plot = oo_plot.get('stys_plot', [':'])
    sel_norm_w = oo_plot.get('sel_norm_w', 'wc')
    sel_norm_t = oo_plot.get('sel_norm_t', 'wc')  # 'wc', 's'
    oo_texts  = oo_plot.get('text', None)

    # calculate data for different configurations
    conf_res = []
    for i_conf in range(n_conf):
        conf_name = oo_cong_names[i_conf]

        print('Calculation: ' + conf_name)

        oo_t   = dict(oo_ts[i_conf])
        oo_var = dict(oo_vars[i_conf])
        oo_wg  = dict(oo_wgs[i_conf])

        dd = oo_var['dds'][0]

        tmin = oo_t['tmin']
        tmax = oo_t['tmax']
        width_t = oo_t['width_t']
        step_t = oo_t['step_t']
        flag_inc_boundary = oo_t.get('flag_inc_boundary', False)
        flag_fv = oo_t.get('flag_fv', False)
        flag_je = oo_t.get('flag_je', False)

        half_width_t = width_t/2

        flag_stat = oo_wg.get('flag_stat', False)
        sel_wg = oo_wg.get('sel_wg', 'wg-adv')

        # normalization:
        coef_norm_w, coef_norm_g, line_norm_w, line_norm_g = choose_wg_normalization(dd, sel_norm_w)

        # time normalization for plotting
        coef_norm_t, line_norm_t = None, ''
        if sel_norm_t == 'wc':
            coef_norm_t, line_norm_t = 1, '[\omega_{ci}]'
        if sel_norm_t == 's':
            coef_norm_t, line_norm_t = 1.e3 / dd['wc'], '(ms)'

        # create time intervals:
        t_right = tmin + width_t
        t_ints  = [[tmin, t_right]]
        t_center = tmin + half_width_t
        while t_right < tmax:
            t_center += step_t

            t_left  = t_center - half_width_t
            t_right = t_center + half_width_t
            if t_right > tmax:
                if not flag_inc_boundary:
                    break
                else:
                    t_right = tmax
            t_ints.append([t_left, t_right])

        # result number of time intervals:
        nt = len(t_ints)

        # precise method of w,g calculation
        if 'adv' in  sel_wg:
            line_res_method = '_adv'
        else:
            line_res_method = '_est'

        line_res = 'naive'
        if flag_stat:
            line_res = 'stat'
            line_res_method = ''

        # --- CALCULATION of w(t) ---
        t = np.zeros(nt)
        ws, ws_err = np.zeros(nt), np.zeros(nt)
        gs, gs_err = np.zeros(nt), np.zeros(nt)
        t[:], ws[:], ws_err[:] = [np.nan]*3
        for i_int in range(nt):
            t_int = t_ints[i_int]
            t[i_int] = t_int[0] + (t_int[1] - t_int[0])/2.

            oo_wg.update({
                't_work': t_int,
            })
            res_wg =  calc_wg(oo_var, oo_wg, {'flag_plot_print': False})

            if flag_stat:
                ws_err[i_int] = res_wg[line_res]['err_w']
                gs_err[i_int] = res_wg[line_res]['err_g']
            ws[i_int] = res_wg[line_res]['w' + line_res_method]

        # --- initial variable ---
        dict_var = choose_vars(oo_var)[0]
        data_var = dict_var['data'][0]

        t_data   = dict_var['x']
        ids_t_data, t_data, _ = mix.get_ids(t_data, [t[0], t[-1]])
        data_var = data_var[ids_t_data[0]:ids_t_data[-1] + 1]

        conf_res_one = {
            'name': conf_name,
            't': t,
            'ws': ws,
            'ws_err': ws_err,
            'flag_fv': flag_fv,
            'flag_je': flag_je,
            'data_var': data_var,
            't_data': t_data,
        }

        # --- RESULTS ---
        conf_res.append(conf_res_one)

    # Exp. NLED AUG:
    exp_w_kHZ = [48, 50.5, 52.25, 53.75, 55, 56, 56.8, 57.4]
    exp_t_s = [0.8410, 0.8415, 0.8420, 0.8425, 0.8430, 0.8435, 0.8440, 0.8445]
    exp_t_ms = (np.array(exp_t_s) - exp_t_s[0]) * 1e3
    exp_w_rel = np.array(exp_w_kHZ) / exp_w_kHZ[0]

    exp_t_ms_cont = np.linspace(exp_t_ms[0], exp_t_ms[-1], 101)
    exp_w_rel_cont = np.interp(exp_t_ms_cont, exp_t_ms, exp_w_rel)

    # --- PLOT RESULTS ---
    if len(stys_plot) < n_conf:
        last_sty = stys_plot[-1]
        stys_plot += [last_sty] * (n_conf - len(stys_plot))

    # plotting
    curves_orig = crv.Curves().xlab('t' + line_norm_t).ylab('signal')
    curves_orig.flag_semilogy = True

    # curves = crv.Curves().xlab('t' + line_norm_t).ylab('\omega' + line_norm_w).tit(tit_plot)
    curves = crv.Curves().xlab('t' + line_norm_t).ylab('\omega/\omega_{0}')
    # additional text:
    for oo_text in oo_texts:
        oText = crv.PlText(oo_text)
        curves.newt(oText)

    curves_exp = crv.Curves().xlab('(t - t_{0})(ms)')\
        .ylab('\omega/\omega_{0}').tit('Experimental\ spectrogram')
    for i_conf in range(n_conf):
        one_conf = conf_res[i_conf]
        curves_orig.new() \
            .XS(one_conf['t_data'] * coef_norm_t) \
            .YS(one_conf['data_var']) \
            .leg(one_conf['name'])
        # curves.new()\
        #     .XS(one_conf['t']  * coef_norm_t)\
        #     .YS(one_conf['ws'] * coef_norm_w) \
        #     .set_errorbar(True, ys=one_conf['ws_err'] * coef_norm_w)\
        #     .leg(one_conf['name'])\
        #     .sty(stys_plot[i_conf])
        curves.new() \
            .XS(one_conf['t'] * coef_norm_t) \
            .YS(one_conf['ws']/one_conf['ws'][0]) \
            .set_errorbar(True, ys=one_conf['ws_err']/one_conf['ws'][0]) \
            .sty(stys_plot[i_conf])

        curves.new() \
            .XS(0.15 * exp_t_ms_cont) \
            .YS(exp_w_rel_cont) \
            .sty(':')
        curves_exp.new()\
            .XS(exp_t_ms)\
            .YS(exp_w_rel)\
            .leg('Long\ branch')\
            .sty('o')
        curves_exp.new() \
            .XS(exp_t_ms_cont) \
            .YS(exp_w_rel_cont) \
            .leg('interpolation') \
            .sty(':')
    cpr.plot_curves(curves_orig)
    cpr.plot_curves(curves)
    cpr.plot_curves(curves_exp)


# NEW: calculation of the w/g in different time intervals:
def calc_wg_t(oo_wg_t, oo_plot):
    # configurations
    oo_cong_names = oo_wg_t['conf_names']
    oo_ts   = oo_wg_t['oo_ts']
    oo_vars = oo_wg_t['oo_vars']
    oo_wgs  = oo_wg_t['oo_wgs']

    # number of configurations
    n_conf = len(oo_cong_names)

    # output properties
    tit_plot = oo_plot.get('tit_plot', None)
    stys_plot = oo_plot.get('stys_plot', [':'])
    sel_norm_w = oo_plot.get('sel_norm_w', 'wc')
    sel_norm_t = oo_plot.get('sel_norm_t', 'wc')  # 'wc', 's'

    # calculate data for different configurations
    conf_res = []
    for i_conf in range(n_conf):
        conf_name = oo_cong_names[i_conf]

        print('Calculation: ' + conf_name)

        oo_t   = dict(oo_ts[i_conf])
        oo_var = dict(oo_vars[i_conf])
        oo_wg  = dict(oo_wgs[i_conf])

        dd = oo_var['dds'][0]

        tmin = oo_t['tmin']
        tmax = oo_t['tmax']
        width_t = oo_t['width_t']
        step_t = oo_t['step_t']
        flag_inc_boundary = oo_t.get('flag_inc_boundary', False)
        species_name = oo_t.get('species_name', None)
        flag_fv = oo_t.get('flag_fv', False)
        flag_je = oo_t.get('flag_je', False)
        s1_q    = oo_t.get('s1_q', None)
        n_analyt_vres = oo_t.get('n_analyt_vres', 1)

        half_width_t = width_t/2

        flag_stat = oo_wg.get('flag_stat', False)
        sel_wg = oo_wg.get('sel_wg', 'wg-adv')

        # normalization:
        coef_norm_w, coef_norm_g, line_norm_w, line_norm_g = choose_wg_normalization(dd, sel_norm_w)

        # time normalization for plotting
        coef_norm_t, line_norm_t = None, ''
        if sel_norm_t == 'wc':
            coef_norm_t, line_norm_t = 1, '[\omega_{ci}]'
        if sel_norm_t == 's':
            coef_norm_t, line_norm_t = 1.e3 / dd['wc'], '(ms)'

        # create time intervals:
        t_right = tmin + width_t
        t_ints  = [[tmin, t_right]]
        t_center = tmin + half_width_t
        while t_right < tmax:
            t_center += step_t

            t_left  = t_center - half_width_t
            t_right = t_center + half_width_t
            if t_right > tmax:
                if not flag_inc_boundary:
                    break
                else:
                    t_right = tmax
            t_ints.append([t_left, t_right])

        # result number of time intervals:
        nt = len(t_ints)

        # precise method of w,g calculation
        if 'adv' in  sel_wg:
            line_res_method = '_adv'
        else:
            line_res_method = '_est'

        line_res = 'naive'
        if flag_stat:
            line_res = 'stat'
            line_res_method = ''

        # --- CALCULATION of w(t) ---
        t = np.zeros(nt)
        ws, ws_err = np.zeros(nt), np.zeros(nt)
        gs, gs_err = np.zeros(nt), np.zeros(nt)
        t[:], ws[:], ws_err[:], gs[:], gs_err[:] = [np.nan]*5
        for i_int in range(nt):
            t_int = t_ints[i_int]
            t[i_int] = t_int[0] + (t_int[1] - t_int[0])/2.

            oo_wg.update({
                't_work': t_int,
            })
            res_wg =  calc_wg(oo_var, oo_wg, {'flag_plot_print': False})

            if flag_stat:
                ws_err[i_int] = res_wg[line_res]['err_w']
                gs_err[i_int] = res_wg[line_res]['err_g']
            ws[i_int] = res_wg[line_res]['w' + line_res_method]
            gs[i_int] = res_wg[line_res]['g' + line_res_method]

        # analytical resonance velocity
        vres = get_v_res(dd, s1_q, ws, species_name)
        curves_draw_vres = geom.Curve()
        for i_vres in range(n_analyt_vres):
            vres_data_1 = vres / (i_vres + 1)
            curves_draw_vres.add_curve(t * coef_norm_t, + vres_data_1)
            curves_draw_vres.add_curve(t * coef_norm_t, - vres_data_1)
        curves_draw_vres.color = 'grey'
        curves_draw_vres.style = ':'
        curves_draw_vres.width = 4

        conf_res_one = {
            'name': conf_name,
            't': t,
            'ws': ws,
            'gs': gs,
            'ws_err': ws_err,
            'gs_err': gs_err,
            'vres': vres,
            'curves_draw_vres': curves_draw_vres,
            'flag_fv': flag_fv,
            'flag_je': flag_je,
        }

        # --- clumps and holes and zeros: parallel velocities, result frequencies ---
        if flag_fv:
            # delta f
            oo_get_f = {'avrs':  [['tvpar', 'none-']], 'dds': [dd],
                        'ovars': [['distribution', 'df_vel_1d', species_name]]}
            vvar_delta_f = choose_vars(oo_get_f)[0]
            deltaf = vvar_delta_f['data']
            t_f = vvar_delta_f['t']
            vpar_f = vvar_delta_f['vpar']

            # d \delta f / dv
            oo_get_f['ovars'] = [['distribution', 'ddeltaf_vel_1d-dv', species_name]]
            vvar_d_delta_f = choose_vars(oo_get_f)[0]
            ddeltaf = vvar_d_delta_f['data']

            # clumps, holes, zeros:
            ids_pos, _, _ = mix.get_ids(vpar_f, [0.0, np.max(vpar_f[-1])])
            ids_pos = range(ids_pos[0], ids_pos[-1]+1)
            vpar_f_pos = vpar_f[ids_pos]

            id_f_max = np.argmax(vvar_delta_f['data'][:, ids_pos], axis=1)
            id_f_min = np.argmin(vvar_delta_f['data'][:, ids_pos], axis=1)
            id_f_zero = np.argmin(vvar_d_delta_f['data'][:, ids_pos], axis=1)

            vpar_max, vpar_min, vpar_zero = \
                np.zeros(np.size(id_f_max)), np.zeros(np.size(id_f_min)), np.zeros(np.size(id_f_zero))
            for id_t in range(np.size(id_f_max)):
                vpar_max[id_t]  = vpar_f_pos[id_f_max[id_t]]
                vpar_min[id_t]  = vpar_f_pos[id_f_min[id_t]]
                vpar_zero[id_t] = vpar_f_pos[id_f_zero[id_t]]
            w_from_max  = get_w_from_v_res(dd, s1_q, vpar_max, species_name)
            w_from_min  = get_w_from_v_res(dd, s1_q, vpar_min, species_name)
            w_from_zero = get_w_from_v_res(dd, s1_q, vpar_zero, species_name)

            curves_draw_vres_max  = geom.Curve()
            curves_draw_vres_min  = geom.Curve()
            curves_draw_vres_zero = geom.Curve()
            curves_draw_vres_max.add_curve(t_f * coef_norm_t, vpar_max)
            curves_draw_vres_min.add_curve(t_f * coef_norm_t, vpar_min)
            curves_draw_vres_zero.add_curve(t_f * coef_norm_t, vpar_zero)
            curves_draw_vres_max.color, curves_draw_vres_min.color, curves_draw_vres_zero.color = \
                'black', 'black', 'black'
            curves_draw_vres_max.style, curves_draw_vres_min.style, curves_draw_vres_zero.style = \
                ':', ':', ':'
            curves_draw_vres_max.width, curves_draw_vres_min.width, curves_draw_vres_zero.width = \
                4, 4, 4

            label_deltaf  = species_name + ':\ \delta f(t, v_{res})'
            label_ddeltaf = species_name + ':\ \partial \delta f(t, v_{res})/\partial v'

            conf_res_one.update({
                'vpar_max':  vpar_max,
                'vpar_min':  vpar_min,
                'vpar_zero': vpar_zero,
                'w_from_max':  w_from_max,
                'w_from_min':  w_from_min,
                'w_from_zero': w_from_zero,
                't_f': t_f,
                'vpar_f': vpar_f,
                'deltaf': deltaf,
                'ddeltaf': ddeltaf,
                'curves_draw_vres_max': curves_draw_vres_max,
                'curves_draw_vres_min': curves_draw_vres_min,
                'curves_draw_vres_zero': curves_draw_vres_zero,
                'label_deltaf': label_deltaf,
                'label_ddeltaf': label_ddeltaf,
            })

        # --- Investigate time evolution of a species J*E signal: <JE(t, vpar)>_{nT_gam} ---
        if flag_je:
            oo_get_je = {
                'ovars': [['mpr', 'je', species_name]],
                'avrs': [['tvpar', 'none-']],
                'dds': [dd],
            }
            vvar_je = choose_vars(oo_get_je)[0]
            t_je = vvar_je['t']
            vpar = vvar_je['vpar']

            je_data = np.zeros([nt, len(vpar)])
            for i_int in range(nt):
                tmin1 = t_ints[i_int][0]
                tmax1 = t_ints[i_int][-1]
                T_loc = 2*np.pi / ws[i_int]
                n_periods = np.floor((tmax1 - tmin1)/T_loc)
                ids_min1, tmin1, _ = mix.get_ids(vvar_je['t'], tmin1)
                ids_max1, tmax1, _ = mix.get_ids(vvar_je['t'],
                                                 tmin1 + T_loc * n_periods
                                                 )
                je_data[i_int, :] = np.squeeze(np.nansum(
                    vvar_je['data'][ids_min1:ids_max1+1, :], axis=0
                ))

            conf_res_one.update({
                'je_tvpar': je_data,
                'je_tje_vpar': vvar_je['data'],
                't_je': t_je,
                'vpar_je': vpar,
                'tit_je': vvar_je['tit'],
            })

        # --- RESULTS ---
        conf_res.append(conf_res_one)

    # --- PLOT RESULTS ---
    if len(stys_plot) < n_conf:
        last_sty = stys_plot[-1]
        stys_plot += [last_sty] * (n_conf - len(stys_plot))

    # frequency
    curves = crv.Curves().xlab('t' + line_norm_t).ylab('\omega' + line_norm_w).tit(tit_plot)
    for i_conf in range(n_conf):
        one_conf = conf_res[i_conf]
        curves.new()\
            .XS(one_conf['t'] * coef_norm_t)\
            .YS(one_conf['ws'] * coef_norm_w)\
            .set_errorbar(True, ys=one_conf['ws_err'] * coef_norm_w)\
            .leg(one_conf['name'])\
            .sty(stys_plot[i_conf])
    cpr.plot_curves(curves)

    # distribution function (t,vpar):
    for i_conf in range(n_conf):
        one_conf = conf_res[i_conf]
        if not one_conf['flag_fv']:
            continue

        vpar_plot, ids_vpar_plot = mix.get_array_oo(
            oo_plot, one_conf['vpar_f'], 'vpar')
        ids_vpar_plot = range(ids_vpar_plot[0], ids_vpar_plot[-1] + 1)

        curves = crv.Curves() \
            .xlab('t' + line_norm_t).ylab('v_{\parallel}') \
            .tit(one_conf['name']).titn(one_conf['label_deltaf'])
        curves.new() \
            .XS(one_conf['t_f'] * coef_norm_t) \
            .YS(vpar_plot) \
            .ZS(one_conf['deltaf'][:, ids_vpar_plot]) \
            .lev(60)
        curves.newg(one_conf['curves_draw_vres'])
        curves.newg(one_conf['curves_draw_vres_max'])
        curves.newg(one_conf['curves_draw_vres_min'])
        cpr.plot_curves_3d(curves)

        curves = crv.Curves() \
            .xlab('t' + line_norm_t).ylab('v_{\parallel}') \
            .tit(one_conf['name']).titn(one_conf['label_ddeltaf'])
        curves.new() \
            .XS(one_conf['t_f'] * coef_norm_t) \
            .YS(vpar_plot) \
            .ZS(one_conf['ddeltaf'][:, ids_vpar_plot]) \
            .lev(60)
        curves.newg(one_conf['curves_draw_vres'])
        curves.newg(one_conf['curves_draw_vres_zero'])
        cpr.plot_curves_3d(curves)

        curves = crv.Curves() \
            .xlab('t' + line_norm_t).ylab('\omega' + line_norm_w) \
            .tit(one_conf['name'])
        curves.new() \
            .XS(one_conf['t_f'] * coef_norm_t) \
            .YS(one_conf['w_from_max'] * coef_norm_w) \
            .leg('clumps').sty('-')
        curves.new() \
            .XS(one_conf['t_f'] * coef_norm_t) \
            .YS(one_conf['w_from_min'] * coef_norm_w) \
            .leg('holes').sty('-')
        curves.new() \
            .XS(one_conf['t_f'] * coef_norm_t) \
            .YS(one_conf['w_from_zero'] * coef_norm_w) \
            .leg('zeros').sty(':')
        curves.new() \
            .XS(one_conf['t'] * coef_norm_t) \
            .YS(one_conf['ws'] * coef_norm_w) \
            .leg('fitting').sty(':o').col('black')
        cpr.plot_curves(curves)

    # <je(t,vpar)>_nT
    for i_conf in range(n_conf):
        one_conf = conf_res[i_conf]
        if one_conf['flag_je']:
            vpar_plot, ids_vpar_plot = mix.get_array_oo(
                oo_plot, one_conf['vpar_je'], 'vpar')

            curves = crv.Curves()\
                .xlab('t' + line_norm_t).ylab('v_{\parallel}')\
                .tit(one_conf['name'] + ':\ <' + one_conf['tit_je'] + '>_{nT_{egam}}')
            curves.new() \
                .XS(one_conf['t'] * coef_norm_t) \
                .YS(vpar_plot) \
                .ZS(one_conf['je_tvpar'][:, ids_vpar_plot[0]:ids_vpar_plot[-1]+1]) \
                .lev(60)
            curves.newg(one_conf['curves_draw_vres'])
            cpr.plot_curves_3d(curves)

    # je(t,vpar)
    for i_conf in range(n_conf):
        one_conf = conf_res[i_conf]
        if one_conf['flag_je']:
            vpar_plot, ids_vpar_plot = mix.get_array_oo(
                oo_plot, one_conf['vpar_je'], 'vpar')

            curves = crv.Curves() \
                .xlab('t' + line_norm_t).ylab('v_{\parallel}') \
                .tit(one_conf['name'] + ':\ ' + one_conf['tit_je'])
            curves.new() \
                .XS(one_conf['t_je'] * coef_norm_t) \
                .YS(vpar_plot) \
                .ZS(one_conf['je_tje_vpar'][:, ids_vpar_plot[0]:ids_vpar_plot[-1] + 1]) \
                .lev(60)
            curves.newg(one_conf['curves_draw_vres'])
            cpr.plot_curves_3d(curves)


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

    oo_avr = dict(oo)
    vvars, xs, tit_vars, lines_avr = [], [], [], []
    for i_var in range(n_vars):
        dd = dds[i_var]
        line_id_var = '{:d}'.format(i_var + 1)

        if i_var == 0:
            oo_avr.update({'flag_var_first': True})
        else:
            oo_avr.update({'flag_var_first': False})
        vvar = average_one_var(avrs[i_var], ovars[i_var], dds[i_var], oo_avr)
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
def average_one_var(avr, ovar, dd, oo):
    # oo.avrs = [[sel_coords type opt1 opt2 ...], [], ...]

    sel_coords  = avr[0]  # chosen coordinate system
    oavr = avr[1:len(avr)]

    # additional parameteres:
    flag_var_first   = oo.get('flag_var_first', False)

    # choose coordinate system, where the system will be considered:
    if sel_coords == 'ts':
        vvar_res = choose_one_var_ts(ovar, dd)
    elif sel_coords == 'tchi':
        vvar_res = choose_one_var_tchi(ovar, dd)
    elif sel_coords == 'tvpar':
        vvar_res = choose_one_var_tvpar(ovar, dd)
    elif sel_coords == 'vparmu':
        vvar_res = choose_one_var_vparmu(ovar, dd)
    elif sel_coords == 't':
        vvar_res = choose_one_var_t(ovar, dd)
        oavr = ['none-']
    elif sel_coords == 'rz':
        vvar_res = choose_one_var_rz(ovar, dd)
    elif sel_coords == 'schi':
        vvar_res = choose_one_var_schi(ovar, dd)
    else:
        print('Wrong name of plane in average_one_var')
        sys.exit(-1)

    if flag_var_first and 'x2' in vvar_res:
        x1_ref = vvar_res[vvar_res['x1']]
        x2_ref = vvar_res[vvar_res['x2']]
        oo.update({vvar_res['x1'] + '_ref': x1_ref,
                   vvar_res['x2'] + '_ref': x2_ref})

    # averaging of the chosen system
    vvar_avr = ymath.avr_x1x2(vvar_res, oavr, oo)

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
    if opt_type == 'fields3d':
        vvar = fields3d.choose_one_var_ts(oovar, dd)

    # save information about coordinate system:
    vvar['x1'], vvar['fx1'], vvar['labx'] = 't', '{:0.3e}', 't'
    vvar['x2'], vvar['fx2'], vvar['laby'] = 's', '{:0.3f}', 's'

    return vvar


# NEW: take variable(r,z)
def choose_one_var_rz(ovar, dd):
    opt_type = ovar[0]
    oovar = ovar[1:len(ovar)]

    vvar = {}
    if opt_type == 'nonzonal':
        vvar = itg.choose_one_var_rz(oovar, dd)
    if opt_type == 'fields3d':
        vvar = fields3d.choose_one_var_rz(oovar, dd)

    # save information about coordinate system:
    vvar['x1'], vvar['fx1'], vvar['labx'] = 'r', '{:0.3f}', 'R'
    vvar['x2'], vvar['fx2'], vvar['laby'] = 'z', '{:0.3f}', 'Z'

    return vvar


# NEW: take variable(r,z)
def choose_one_var_schi(ovar, dd):
    opt_type = ovar[0]
    oovar = ovar[1:len(ovar)]

    vvar = {}
    if opt_type == 'nonzonal':
        vvar = itg.choose_one_var_rz(oovar, dd)
    if opt_type == 'fields3d':
        vvar = fields3d.choose_one_var_rz(oovar, dd)

    # save information about coordinate system:
    vvar['x1'], vvar['fx1'], vvar['labx'] = 's', '{:0.3f}', 's'
    vvar['x2'], vvar['fx2'], vvar['laby'] = 'chi', '{:0.3f}', '\chi'

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
    if opt_type.lower() == 'distribution':
        vvar = distribution.choose_one_var_tvpar(oovar, dd)
    if opt_type.lower() == 'mpr':
        vvar = MPR.choose_one_var_tvpar(oovar, dd)

    # save information about coordinate system:
    vvar['x1'], vvar['fx1'], vvar['labx'] = 't',    '{:0.3e}', 't'
    vvar['x2'], vvar['fx2'], vvar['laby'] = 'vpar', '{:0.3f}', 'v_{\parallel}'

    return vvar


# NEW: take variable(vpar,mu)
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


# NEW: take variable(t)
def choose_one_var_t(ovar, dd):
    opt_type = ovar[0]
    oovar = ovar[1:len(ovar)]

    vvar = {}
    if opt_type == 'mpr':
        vvar = MPR.choose_one_var_t(oovar, dd)

    # save information about coordinate system:
    vvar['fx'], vvar['labx'] = '{:0.3e}', 't'
    vvar['laby'] = ''
    vvar['name_x'] = 't',
    vvar['format_x'] =  vvar['fx'],

    vvar['data'] = [vvar['data']]  # to be consistent with other signals after averaging

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
def MPR_gamma_velocity_domains(dd, oo_vel, oo_vars, oo_wg, oo_plot):
    # --------------------------------------------------------------
    # For description of the dictionaries oo_vars, oo_wg, oo_plot
    # see the function MPR_gamma()
    # --------------------------------------------------------------
    # Dictionary oo_vars does not have fields: 'mu_domain', 'vpar_domain'.
    # Dictionary oo_vars also has the following fields:
    #   't_int' = [t_start, t_end]:
    #       time interval, where je(vpar, mu, t) will be averaged
    # --------------------------------------------------------------
    # Dictionary oo_wg also has the following fields:
    #  's1' - integer:
    #       radial point, where safety factor is taken and the gam/egam
    #       frequency is taken
    # ------------------------------------------------------------------
    # oo_vel - [dict., dict., ...], see code, section Velocity domains
    # ------------------------------------------------------------------
    # oo_plot - dict. for plotting:
    #   'mu_plot' - domain of mu to plot
    #   'vpar_plot' - domain of vpar to plot
    #   'texts_plot' - [dict., dict., ...]:
    #       lines of text to add to a plot;
    #       for description of every dict. see
    #       file curve.py, class PlText, function init_from_oo(dict.)
    #   'geoms_plot' - [dict., dict., ...]:
    #       geometrical figures to add to a plot;
    #       for description of every dict. see
    #       file Geom.py, class Geom and classed, derived from it, function

    # --- Parameters ---
    je_name = 'jdote_es-mean-t'
    sel_species = oo_vars.get('sel_species', 'total')
    t_int = oo_vars.get('t_int', [])

    s1 = oo_wg.get('s1', None)
    w0_wc = oo_wg.get('gam-w', None)

    mu_plot = oo_plot.get('mu_plot', [])
    vpar_plot = oo_plot.get('vpar_plot', [])
    n_vpar_res = oo_plot.get('n_vpar_res', 0)
    oo_texts = oo_plot.get('texts_plot', [])
    oo_geoms = oo_plot.get('geoms_plot', [])
    tit_vmu = oo_plot.get('tit_vmu', None)
    flag_pass_trap_cone = oo_plot.get('flag_pass_trap_cone', True)

    # velocity normalization:
    coef_max = np.max(dd[sel_species].nT_equil['T'])
    Te_speak = dd['electrons'].T_speak(dd)
    cs_speak = np.sqrt(Te_speak / dd['pf'].mass)
    norm_v = coef_max / cs_speak

    # Te_max = np.max(dd['electrons'].nT_equil['T_J'])
    # cs_max = np.sqrt(Te_max / dd['pf'].mass)
    # norm_v = 1. / cs_max

    # --- Species energy transfer signal ---
    oo_vvar = {
        'avrs': [['vparmu', 'none-']],
        'dds': [dd],
    }
    ovar_array = [['mpr', je_name, sel_species, t_int]]
    oo_vvar.update({'ovars': ovar_array})
    je_dict = choose_vars(oo_vvar)[0]
    je      = np.array(je_dict['data'])
    mu      = np.array(je_dict['mu'])
    vpar    = np.array(je_dict['vpar'])
    nmu = np.size(mu)

    # --- Analytical resonances ---
    w0 = w0_wc * dd['wc']

    q_s1, _, line_s1 = equil_profiles.q_s1(dd, s1)
    vres = (q_s1 * dd['R0'] * w0) * norm_v

    gHLines = None
    if n_vpar_res is not 0:
        gHLines = geom.HLine()
        gHLines.ys = []
        for i_res in range(n_vpar_res):
            # gHLines.ys = [-vres, -vres/2, vres/2, vres]
            gHLines.ys += [-vres/(i_res+1), vres/(i_res+1)]
        gHLines.color = 'black'
        gHLines.style = '--'
        gHLines.width = 4

    # --- Passing-Trapped boundary ---
    cone_pt, mu_pt, pt_bottom, pt_up = geom.pass_trap_boundary(dd, mu)

    # --- Velocity domains ---
    n_areas = len(oo_vel)
    names_areas = []
    je_areas = []
    ids_je_areas = []
    separate_boundaries_areas = []
    vertical_lines_areas = []
    for i_area in range(n_areas):
        oo_area = oo_vel[i_area]
        name_type = oo_area.get('type', None)

        sep_boundary_area = None

        vertical_lines_areas.append(None)

        name_area = oo_area.get('name', '')

        # - strap between the passing-trapped boundary and a chosen parabola -
        if name_type == 'strap_pr_pb':
            # build a parabola mu = a * vpar^2 + mu0
            mu0, mu1, vpar1 = oo_area['mu0'], oo_area['mu1'], oo_area['vpar1']
            _, _, _, pb_bottom, pb_up = \
                geom.parabola_x(mu, mu0, mu1, vpar1)

            # boundaries
            sep_boundary_area = [pt_bottom, pb_bottom, pb_up, pt_up]
        # ---
        # - area inside a parabola =
        if name_type == 'parabola':
            # build a parabola mu = a * vpar^2 + mu0
            mu0, mu1, vpar1 = oo_area['mu0'], oo_area['mu1'], oo_area['vpar1']
            _, _, _, pb_bottom, pb_up = \
                geom.parabola_x(mu, mu0, mu1, vpar1)

            # boundaries
            sep_boundary_area = [pb_bottom, pb_up]
        # ---
        # - area inside a vertical parallelogram =
        if name_type == 'parallelogram_v':
            flag_upper_points = oo_area.get('flag-upper-points', True)
            point_left, point_right = oo_area['point-left'], oo_area['point-right']
            vlength = oo_area['vlength']

            pav_bottom, pav_up, pav_corners = geom.parallelogram_v(
                mu, point_left, point_right, vlength, flag_upper_points
            )

            # boundaries
            sep_boundary_area = [pav_bottom, pav_up]

            # vertical lines
            vertical_line_left  = [pav_corners[2], pav_corners[0]]
            vertical_line_right = [pav_corners[3], pav_corners[1]]
            vertical_lines_areas[-1] = [vertical_line_left, vertical_line_right]
        if name_type == 'outside-lines':
            lines = oo_area.get('lines', [])

            # lines
            lines_area = geom.lines(mu, lines)

            # bottom and upper vpar-boundaries
            line_low_vpar = geom.lines(mu, [ [[mu[0],  vpar[0]], [mu[-1],  vpar[0]]]  ])
            line_up_vpar  = geom.lines(mu, [ [[mu[0], vpar[-1]], [mu[-1], vpar[-1]]] ])

            # boundaries

            sep_boundary_area = line_low_vpar + lines_area + line_up_vpar

        # ---
        # - build boundaries -
        boundaries_area = geom.create_boundaries(
            nmu,
            sep_boundary_area
        )

        # - build J*E inside the strap -
        je_area, ids_area_je = geom.build_area(mu, vpar, je, boundaries_area)

        # - results -
        separate_boundaries_areas.append(sep_boundary_area)
        je_areas.append(je_area)
        ids_je_areas.append(ids_area_je)
        names_areas.append(name_area)

    # --- Plot chosen velocity domains ---
    if len(mu_plot) is 0 or mu_plot is None:
        mu_plot = mu
    else:
        _, mu_plot, _ = mix.get_ids(mu, mu_plot)
    if len(vpar_plot) is 0 or vpar_plot is None:
        vpar_plot = vpar
    else:
        _, vpar_plot, _ = mix.get_ids(vpar, vpar_plot)

    tit_vmu_res = tit_vmu
    if tit_vmu_res is not None:
        if len(tit_vmu_res) is 0:
            tit_vmu_res = je_dict['tit']

    colMap, colMap_center = 'seismic', 0

    color_area = 'white'
    curves = crv.Curves()\
        .xlab('\mu').ylab('v_\parallel') \
        .tit(tit_vmu_res)
    curves.flag_legend = False
    curves.newg(gHLines)
    curves.xlim([mu_plot[0],   mu_plot[-1]])
    curves.ylim([vpar_plot[0], vpar_plot[-1]])
    curves.new() \
        .XS(mu)\
        .YS(vpar)\
        .ZS(je).lev(60) \
        .cmp(colMap, colMap_center)
    for i_area in range(n_areas):
        boundaries_area = separate_boundaries_areas[i_area]
        vertical_lines = vertical_lines_areas[i_area]
        for sel_boundary in boundaries_area:
            curves.new() \
                .XS(mu) \
                .YS(sel_boundary) \
                .col(color_area).sty(':')
        if vertical_lines is not None:
            gVerticalLinesArea = geom.Line()
            for vertical_line in vertical_lines:
                gVerticalLinesArea.lines.append(vertical_line)
                gVerticalLinesArea.color = color_area
                gVerticalLinesArea.style = ':'
                gVerticalLinesArea.width = 3
            curves.newg(gVerticalLinesArea)
    for oo_text in oo_texts:
        oText = crv.PlText(oo_text)
        curves.newt(oText)
    for oo_geom in oo_geoms:
        oGeom = None
        if oo_geom['type'] == 'annotate':
            oGeom = geom.Annotate(oo_geom)
        if oGeom is not None:
            curves.newg(oGeom)
    if flag_pass_trap_cone:
        curves.new() \
            .XS(mu_pt) \
            .YS(cone_pt) \
            .col('black').sty(':').w(20)
    cpr.plot_curves_3d(curves)

    # --- Plot areas ---
    if gHLines is not None:
        gHLines.color = 'black'
    for i_area in range(n_areas):
        curves = crv.Curves() \
            .xlab('\mu').ylab('v_\parallel') \
            .tit(je_dict['tit'] + ':\ ' + names_areas[i_area])
        curves.flag_legend = False
        curves.newg(gHLines)
        curves.xlim([mu_plot[0],   mu_plot[-1]])
        curves.ylim([vpar_plot[0], vpar_plot[-1]])
        curves.newg(gHLines)
        curves.new() \
            .XS(mu) \
            .YS(vpar) \
            .ZS(je_areas[i_area]).lev(60) \
            .cmp('seismic')
        curves.new() \
            .XS(mu_pt) \
            .YS(cone_pt) \
            .col('black').sty('-')
        cpr.plot_curves_3d(curves)

    # --- Calculate gamma in chosen velocity domains ---
    oo_vars['flag_given_signal_je'] = True
    for i_area in range(n_areas):
        ovar_loc = ['jdote_es', sel_species, ids_je_areas[i_area], names_areas[i_area]]
        oo_vars['signal_je'] = MPR.choose_one_var_t(ovar_loc, dd, flag_vpar_boundaries=True)
        oo_vars['signal_je']['legs'] = ['_']
        MPR_gamma(dd, oo_vars, oo_wg, oo_plot)

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
    #   'flag_given_signal_je' = True:
    #       take a given J*E signal oo_vars['signal_je']
    # -> oo_plot - dict.:
    #   't_plot' - 1d array [t_start, t_end]:
    #       time domain, where data will be plotted
    # -> oo_wg - dict.: parameters for MPR.GAM_calc_wg() diagnostic:
    #   't_work' - 1d array [t_start, t_end]:
    #       time interval, where gamma will be calculated
    #   'filt_je' - dict. or [dict., dict., ...]:
    #       filter of the J*E signal
    #   'filt_ef' - dict. or [dict., dict., ...]:
    #       filter of the field energy signal
    #   REMARK: final time interval is defined by the last filtering domain.
    #           Work domain is taken from this final time domain.
    # ---
    #   'flag_naive_t_peaks' = True: only for naive calculation
    #       take a time interval between peaks to calculate the damping/growth rate,
    #       otherwise, take an arbitrary time interval
    #   'naive_n_periods' - integer:
    #       number of GAM periods that will define a length of the time interval
    #   'sel_norm': 'wc', 'vt':
    #       output normalization
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

    # - None-filter -
    non_filt = {'sel_filt': None}

    # - FUNCTION: filtering at one stage -
    def several_filters(x, y, oo_filt_loc):
        oo_filts_res = []
        if oo_filt_loc is None or len(oo_filt_loc) is 0:
            oo_filts_res.append(non_filt)
        elif isinstance(oo_filt_loc, dict):
            oo_filts_res.append(oo_filt_loc)
        else:
            oo_filts_res = oo_filt_loc  # array of filters

        # apply one by one the filters in the array :
        dict_loc = {'x': np.array(x), 'filt': np.array(y)}
        count_filt = -1
        for one_filt in oo_filts_res:
            count_filt += 1
            dict_loc = ymath.filtering(
                dict_loc['x'], dict_loc['filt'], one_filt
            )
        return dict_loc

    # -FUNCTION: Filter a signal -
    def filter_signal(init_dict, oo_filt, oo_plot):
        # initial fft
        dict_fft_ini = ymath.filtering(init_dict['x'], init_dict['data'], non_filt)
        leg_data = init_dict['legs'][0]

        # filtering
        res_dict = dict(init_dict)
        filt = several_filters(init_dict['x'], init_dict['data'], oo_filt)
        res_dict['data'] = np.array(filt['filt'])
        res_dict['x']    = np.array(filt['x'])

        # find peaks
        if not flag_inv_peaks:
            ids_peaks, _ = scipy.signal.find_peaks(res_dict['data'])
        else:
            ids_peaks, _ = scipy.signal.find_peaks(-res_dict['data'])
        t_peaks = res_dict['x'][ids_peaks]

        # plotting time evolution
        t_plot = oo_plot.get('t_plot', [])
        if t_plot is None or len(t_plot) is 0:
            t_plot = init_dict['x']
        else:
            _, t_plot, _ = mix.get_ids(init_dict['x'], t_plot)

        line_add_to_title = oo_plot.get('line_add_to_title', '')

        curves = crv.Curves() \
            .xlab('t[\omega_{ci}^{-1}]') \
            .tit(init_dict['tit'] + ':\ ' + line_add_to_title)
        curves.flag_norm        = flag_norm
        curves.flag_semilogy    = flag_semilogy
        if flag_norm:
            curves.ylab('norm.')
        curves.xlim([t_plot[0], t_plot[-1]])
        curves.new() \
            .XS(init_dict['x']) \
            .YS(init_dict['data']) \
            .leg('initial\ signal')
        curves.new().norm_to(init_dict['data']) \
            .XS(res_dict['x']) \
            .YS(res_dict['data']) \
            .leg('filtered\ signal').sty(':')
        # curves.new().norm_to(init_dict['data']) \
        #     .XS(t_peaks) \
        #     .YS(res_dict['data'][ids_peaks]) \
        #     .leg('peaks').sty('o')
        cpr.plot_curves(curves)

        # plotting FFT
        curves = crv.Curves() \
            .xlab('\omega[\omega_{ci}]') \
            .tit(leg_data + ':\ FFT'  + ':\ ' + line_add_to_title)
        curves.flag_norm = flag_norm
        if flag_norm:
            curves.ylab('norm.')
        curves.new() \
            .XS(dict_fft_ini['w2']) \
            .YS(dict_fft_ini['fft_init_2']) \
            .leg('FFT:\ initial')
        curves.new() \
            .XS(filt['w2']) \
            .YS(filt['fft_filt_2']) \
            .leg('FFT:\ filtered').sty(':')
        cpr.plot_curves(curves)

        return res_dict

    # --- Parameters ---
    sel_species = oo_vars.get('sel_species', 'total')
    mu_domain   = oo_vars.get('mu_domain', [])
    vpar_domain = oo_vars.get('vpar_domain', [])

    flag_ZFZF = oo_wg.get('flag_ZFZF', False)

    flag_inv_peaks = oo_wg.get('flag_inv_peaks', False)

    oo_filt_je = oo_wg.get('filt_je', None)
    oo_filt_ef = oo_wg.get('filt_ef', None)
    oo_filt_er = oo_wg.get('filt_er', None)  # !!! don't filter out residue component !!!

    flag_norm = oo_plot.get('flag_norm', False)
    flag_semilogy = oo_plot.get('flag_semilogy', False)

    flag_given_signal_je = oo_vars.get('flag_given_signal_je', False)
    given_signal = oo_vars.get('signal_je', None)

    # --- Names ---
    je_name = 'jdote_es'
    ef_name = 'efield'

    # --- GET JE and Efield (t) ---
    oo_vvar = {
        'avrs': [['t']],
        'dds': [dd],
    }
    ovar_array = [['mpr', je_name, sel_species, mu_domain, vpar_domain]]

    if not flag_given_signal_je:
        oo_vvar.update({'ovars': ovar_array})
        je_dict = choose_vars(oo_vvar)[0]
    else:
        je_dict = given_signal

    ovar_array[0][1] = ef_name
    ovar_array[0][2] = 'total'
    oo_vvar.update({'ovars': ovar_array})
    ef_dict = choose_vars(oo_vvar)[0]

    del oo_vvar, ovar_array

    # --- FILTERING ---
    if oo_filt_je is not None:
        je_dict = filter_signal(je_dict, oo_filt_je, oo_plot)
    if oo_filt_ef is not None:
        ef_dict = filter_signal(ef_dict, oo_filt_ef, oo_plot)

    # --- INTERPOLATE Efield to the t-axis of JE ---
    ef_dict['data'] = np.interp(je_dict['x'], ef_dict['x'], ef_dict['data'])
    ef_dict['x'] = np.array(je_dict['x'])

    # - update dictionary oo_wg -
    oo_wg.update({
        'je_dict': je_dict,
        'ef_dict': ef_dict
    })

    # --- EXCLUDE ZFZF from signals ---
    if flag_ZFZF:
        flag_norm     = oo_plot.get('flag_norm', False)
        flag_semilogy = oo_plot.get('flag_semilogy', False)

        t  = je_dict['x']

        t_plot_boundaries = oo_plot.get('t_plot', [])
        if len(t_plot_boundaries) is 0:
            t_plot_boundaries = t
        ids_t_plot, _, _ = mix.get_ids(t, t_plot_boundaries)

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
            'dds': [dd],
            'sel_legs1': 'woPr',
        }
        erbar_dict = choose_vars(oo_erbar)[0]

        erbar_dict_to_filt = dict(erbar_dict)
        erbar_dict_to_filt['data'] = np.array(erbar_dict['data'][0])
        erbar_dict = filter_signal(erbar_dict_to_filt, oo_filt_er, oo_plot)
        erbar = np.interp(t, erbar_dict['x'], erbar_dict['data'])

        # erbar = np.interp(t, erbar_dict['x'], erbar_dict['data'][0])

        # - estimate GAM, ZFZF amplitdes in Erbar -
        res_e = estimate_GAM_ZFZF_amplitudes(gam_g, erbar, t, t_point_inf)

        # - Exclude ZFZF from Ef and JE -
        eta = res_e['A_zfzf'] / res_e['A_gam']
        res_f = MPR_exclude_ZFZF(gam_w, gam_g, ef, je, t, id_peak, eta)
        id_peak = res_f['id_peak']

        # - PLOTTING: electric field -
        curves = crv.Curves().xlab('t[\omega_{ci}^{-1}]').tit(erbar_dict['tit'])
        curves.flag_norm = flag_norm
        curves.flag_semilogy = flag_semilogy
        curves.new() \
            .XS(t).YS(erbar) \
            .leg('ZF + GAM')
        curves.new() \
            .XS(t).YS(res_e['e_gam_cos']) \
            .leg('GAM:\ cos')
        curves.new().norm_to(res_e['e_gam_cos']) \
            .XS(res_e['t_peaks']).YS(res_e['e_gam_peaks']) \
            .leg('GAM:\ cos:\ peaks').sty('o')
        curves.new() \
            .XS(t).YS(res_e['e_gam']) \
            .leg('GAM')
        cpr.plot_curves(curves)

        # - PLOTTING: field energy -
        curves = crv.Curves()\
            .xlab('t[\omega_{ci}^{-1}]')\
            .tit(ef_dict['tit'])
        if flag_norm:
            curves.ylab('norm.')
        curves.flag_norm = flag_norm
        curves.flag_semilogy = flag_semilogy
        curves.new() \
            .XS(t[ids_t_plot[0]:ids_t_plot[-1]+1]) \
            .YS(res_e['e_gam'][ids_t_plot[0]:ids_t_plot[-1]+1])\
            .leg(erbar_dict['legs'][0]).sty(':')
        curves.new()\
            .XS(t[ids_t_plot[0]:ids_t_plot[-1]+1])\
            .YS(ef[ids_t_plot[0]:ids_t_plot[-1]+1])\
            .leg('\mathcal{E}:\ ZFZF+GAM')
        curves.new().norm_to(ef[ids_t_plot[0]:ids_t_plot[-1]+1]) \
            .XS(t[res_f['ids_peaks']]) \
            .YS(ef[res_f['ids_peaks']]) \
            .leg('peaks').sty('o')
        curves.new().norm_to(ef[ids_t_plot[0]:ids_t_plot[-1]+1]) \
            .XS(t[res_f['ids_peaks'][id_peak]]) \
            .YS(ef[res_f['ids_peaks'][id_peak]]) \
            .leg('chosen\ peak').sty('s')
        curves.new()\
            .XS(t[ids_t_plot[0]:ids_t_plot[-1]+1])\
            .YS(res_f['Egam'][ids_t_plot[0]:ids_t_plot[-1]+1])\
            .leg('\mathcal{E}:\ GAM').sty(':')
        cpr.plot_curves(curves)

        # - PLOTTING: energy transfer signal -
        curves = crv.Curves()\
            .xlab('t[\omega_{ci}^{-1}]')\
            .tit(je_dict['tit'])
        if flag_norm:
            curves.ylab('norm.')
        curves.flag_norm = flag_norm
        curves.flag_semilogy = flag_semilogy
        curves.new()\
            .XS(t[ids_t_plot[0]:ids_t_plot[-1]+1])\
            .YS(je[ids_t_plot[0]:ids_t_plot[-1]+1])\
            .leg('ZF+GAM')
        curves.new()\
            .XS(t[ids_t_plot[0]:ids_t_plot[-1]+1])\
            .YS(res_f['Pgam'][ids_t_plot[0]:ids_t_plot[-1]+1])\
            .leg('GAM').sty(':')
        cpr.plot_curves(curves)

        # - Comparison with theory -
        th_gam = (gam_w * np.sin(2*gam_w*t) - gam_g - gam_g * np.cos(2*gam_w*t)) * \
                 np.exp(2 * gam_g * t)
        th_zf = th_gam + \
                (gam_w * np.sin(gam_w * t) - gam_g * np.cos(gam_w * t)) * np.exp(gam_g * t)

        curves = crv.Curves() \
            .xlab('t[\omega_{ci}^{-1}]') \
            .ylab('norm.') \
            .tit('\mathcal{P}\ with\ only\ GAM\ component')
        curves.flag_norm = True
        curves.new()\
            .XS(t[ids_t_plot])\
            .YS(th_gam[ids_t_plot])\
            .leg('Analytical\ theory')
        curves.new()\
            .XS(t[ids_t_plot])\
            .YS(res_f['Pgam'][ids_t_plot])\
            .leg('ORB5').sty(':')
        cpr.plot_curves(curves)

        curves = crv.Curves() \
            .xlab('t[\omega_{ci}^{-1}]') \
            .ylab('norm.') \
            .tit('\mathcal{P}\ with\ GAM\ and\ ZFZF\ components')
        curves.flag_norm = True
        curves.new() \
            .XS(t[ids_t_plot]) \
            .YS(th_zf[ids_t_plot]) \
            .leg('Analytical\ theory')
        curves.new() \
            .XS(t[ids_t_plot]) \
            .YS(je[ids_t_plot]) \
            .leg('ORB5').sty(':')
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
    MPR.calc_gamma(dd, oo_wg, oo_plot, oo_desc)

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

    # if eta > 0:
    #     idd = np.array([E_peak1, E_peak2]).argmax()
    # else:
    #     # in this case you should take the peak at the beginning where
    #     # there are still alternating peaks
    #     idd = np.array([E_peak1, E_peak2]).argmin()
    # if idd == 1:
    #     id_peak += 1

    if eta < 0:
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


# NEW: Test animation:
def animation_1d(oo):
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
    labx = oo.get('labx', None)  # x-label
    laby = oo.get('laby', None)  # y-label
    tit_plot = oo.get('tit_plot', None)  # title

    xlims = oo.get('xlims', None)
    ylims = oo.get('ylims', None)

    flag_norm = oo.get('flag_norm', False)
    flag_semilogy = oo.get('flag_semilogy', False)

    # plotting:
    curves = crv.Curves().xlab(labx).ylab(laby).tit(tit_plot)
    curves.flag_norm     = flag_norm
    curves.flag_semilogy = flag_semilogy
    curves.xlim(xlims)
    curves.ylim(ylims)
    for ivar in range(n_vars):
        vvar = vvars[ivar]
        leg = vvar['legs'][0]
        if vvar['x1'] is 'r' and vvar['x2'] is 'z':
            _, ids_s = mix.get_array_oo(oo, vvar['s'], 's')
            _, ids_chi = mix.get_array_oo(oo, vvar['chi'], 'chi')
            x1 = mix.get_slice(vvar['r'], ids_chi, ids_s)
            x2 = mix.get_slice(vvar['z'], ids_chi, ids_s)
            data = mix.get_slice(vvars[ivar]['data'], ids_s, ids_chi)
        else:
            x1, ids_x1 = mix.get_array_oo(oo, vvar[vvar['x1']], vvar['x1'])
            x2, ids_x2 = mix.get_array_oo(oo, vvar[vvar['x2']], vvar['x2'])
            data = mix.get_slice(vvars[ivar]['data'], ids_x1, ids_x2)

        curves.new().XS(x1).YS(x2).ZS(data).lev(60).leg(leg)
    cpr.animation_1d(curves)


# NEW: resonant velocity:
def get_v_res(dd, s1, w_wci, sel_species):
    coef_max = np.max(dd[sel_species].nT_equil['T'])
    Te_speak = dd['electrons'].T_speak(dd)
    cs_speak = np.sqrt(Te_speak / dd['pf'].mass)
    norm_v = coef_max / cs_speak

    # Te_max = np.max(dd['electrons'].nT_equil['T_J'])
    # cs_max = np.sqrt(Te_max / dd['pf'].mass)
    # norm_v = 1 / cs_max

    w0 = w_wci * dd['wc']

    q_s1, _, line_s1 = equil_profiles.q_s1(dd, s1)
    vres = (q_s1 * dd['R0'] * w0) * norm_v

    return vres


# NEW: get frequency from a resonant velocity:
def get_w_from_v_res(dd, s1, vres, sel_species):
    coef_max = np.max(dd[sel_species].nT_equil['T'])
    Te_speak = dd['electrons'].T_speak(dd)
    cs_speak = np.sqrt(Te_speak / dd['pf'].mass)
    norm_v = coef_max / cs_speak

    # Te_max = np.max(dd['electrons'].nT_equil['T_J'])
    # cs_max = np.sqrt(Te_max / dd['pf'].mass)
    # norm_v = 1 / cs_max

    q_s1, _, line_s1 = equil_profiles.q_s1(dd, s1)
    w0 = vres / (q_s1 * dd['R0'] * norm_v)
    w_wci = w0 / dd['wc']

    return w_wci


# get normalization:
def choose_wg_normalization(dd, sel_norm):
    coef_norm_w, coef_norm_g, line_norm_w, line_norm_g = \
        None, None, '', ''
    if sel_norm.lower() == 'wc':
        line_norm_w = line_norm_g = '[\omega_{ci}]'
        coef_norm_w = coef_norm_g = 1
    if sel_norm.lower() == 'vt':
        line_norm_w = line_norm_g = '[sqrt(2)*v_{th,i}/R_0]'
        coef_norm_w = coef_norm_g = \
            dd['wc'] / (np.sqrt(2) * dd['vt'] / dd['R0'])
    if sel_norm.lower() == 'khz':
        line_norm_w = '(kHz)'
        line_norm_g = '(10^3 s)'
        coef_norm_w = dd['wc'] / (1e3 * 2 * np.pi)
        coef_norm_g = dd['wc'] / 1e3
    return coef_norm_w, coef_norm_g, line_norm_w, line_norm_g


# Test continuous wavelet transform
def test_cwt():
    # # Static frequencies
    # t = np.linspace(0, 143, 401)
    #
    # T1, T2 = 5., 15
    # w1, w2 = 1. / T1, 1./T2
    # dt = np.max(t) / (len(t)-1)
    # y = np.cos(2*np.pi * w1 * t) + np.cos(2*np.pi * w2 * t)
    #
    # print('w1(Hz) = {:0.3e}'.format(w1))
    # print('w2(Hz) = {:0.3e}'.format(w2))
    #
    # curves = crv.Curves().xlab('t')
    # curves.new().XS(t).YS(y)
    # cpr.plot_curves(curves)
    #
    # widths = np.arange(1, 16)
    # cwt_res, freqs = pywt.cwt(y, widths, 'mexh', sampling_period=dt)
    # curves = crv.Curves().xlab('t').ylab('freq(Hz)')
    # curves.new().XS(t).YS(freqs).ZS(cwt_res.T)
    # cpr.plot_curves_3d(curves)

    # Evolving frequencies
    t = np.linspace(0, 143, 401)
    dt = np.max(t) / (len(t) - 1)

    T1 = 15
    w01 = 1. / T1
    rate_w1 = w01 / 40.
    w1 = w01 + rate_w1 * t
    y = np.cos(2 * np.pi * w1 * t)

    print('w1(Hz) = {:0.3e}'.format(w01))

    curves = crv.Curves().xlab('t')
    curves.new().XS(t).YS(y)
    cpr.plot_curves(curves)

    curves = crv.Curves().xlab('t').ylab('w(Hz)')
    curves.new().XS(t).YS(w1)
    cpr.plot_curves(curves)

    # opt 1
    widths = np.arange(1, 25)
    cwt_res, freqs = pywt.cwt(y, widths, 'mexh', sampling_period=dt)
    curves = crv.Curves().xlab('t').ylab('freq(Hz)')
    curves.new().XS(t).YS(freqs).ZS(cwt_res.T)
    cpr.plot_curves_3d(curves)


# - FUNCTION: filtering at one stage -
def global_several_filters(x, y, oo_filt_loc):
    oo_filts_res = []
    if oo_filt_loc is None or len(oo_filt_loc) is 0:
        oo_filts_res.append(GLOBAL_NONE_FILTER)
    elif isinstance(oo_filt_loc, dict):
        oo_filts_res.append(oo_filt_loc)
    else:
        oo_filts_res = oo_filt_loc  # array of filters

    # apply one by one the filters in the array :
    dict_loc = {'x': np.array(x), 'filt': np.array(y)}
    count_filt = -1
    for one_filt in oo_filts_res:
        count_filt += 1
        dict_loc = ymath.filtering(
            dict_loc['x'], dict_loc['filt'], one_filt
        )
    return dict_loc


# -FUNCTION: Filter a signal -
def global_filter_signal(t, y, oo_filt):
    # initial fft
    dict_fft_ini = ymath.filtering(t, y, GLOBAL_NONE_FILTER)
    init_dict = {'x': t, 'data': y}

    # filtering
    res_dict = dict(init_dict)
    filt = global_several_filters(init_dict['x'], init_dict['data'], oo_filt)
    res_dict['data'] = np.array(filt['filt'])
    res_dict['x']    = np.array(filt['x'])

    return res_dict['data'], res_dict['x']



