import numpy as np
import importlib as imp
from scipy import signal


def reload():
    # Important: put here all modules that you want to reload
    return


def reload_module(obj_module):
    imp.reload(obj_module)
    obj_module.reload()


def get_attribute(ff, path):
    list_attrs = ff[path].attrs
    ids_attr = list(list_attrs)
    list_attrs = [list_attrs[name].decode("utf-8") for name in ids_attr]
    return list_attrs


def find(x, x1):
    # works only for monotonic increasing arrays!!!
    id_x1 = np.where(x >= x1)[0]
    if id_x1.size != 0:
        id_x1 = id_x1[0]
        x1 = x[id_x1]
    else:
        id_x1 = None
        x1 = None
    return id_x1, x1


def get_array(x, x_lower, x_upper):
    id_start = np.where(x >= x_lower)[0]
    if id_start.size != 0:
        id_start = id_start[0]
    else:
        id_start = 0

    id_end = np.where(x <= x_upper)[0]
    if id_end.size != 0:
        id_end = id_end[-1]
    else:
        id_end = x.size
    ids = np.array([id_start, id_end])
    return x[ids[0]:ids[-1]+1], ids


def get_array_oo(oo, x, name_x):
    x_start = oo.get(name_x + '_start', x[0])
    x_end = oo.get(name_x + '_end', x[-1])
    return get_array(x, x_start, x_end)  # x_new, ids_x


def get_slice(x, ids1, ids2=None, ids3=None):
    # 1D
    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, type(None)) and \
            isinstance(ids3, type(None)):
        return x[ids1]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, type(None)) and \
            isinstance(ids3, type(None)):
        return x[ids1[0]:ids1[-1]+1]

    # 2D
    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (int, np.int64)) \
            and isinstance(ids3, type(None)):
        return x[ids1, ids2]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (int, np.int64)) and \
            isinstance(ids3, type(None)):
        return x[ids1[0]:ids1[-1]+1, ids2]

    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, type(None)):
        return x[ids1, ids2[0]:ids2[-1]+1]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, type(None)):
        return x[ids1[0]:ids1[-1]+1, ids2[0]:ids2[-1]+1]

    # 3D
    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (int, np.int64)) \
            and isinstance(ids3, (int, np.int64)):
        return x[ids1, ids2, ids3]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (int, np.int64)) and \
            isinstance(ids3, (int, np.int64)):
        return x[ids1[0]:ids1[-1]+1, ids2, ids3]

    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, (int, np.int64)):
        return x[ids1, ids2[0]:ids2[-1]+1, ids3]

    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (int, np.int64)) and \
            isinstance(ids3, (np.ndarray, list)):
        return x[ids1, ids2, ids3[0]:ids3[-1]+1]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, (int, np.int64)):
        return x[ids1[0]:ids1[-1]+1, ids2[0]:ids2[-1]+1, ids3]

    if isinstance(ids1, (int, np.int64)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, (np.ndarray, list)):
        return x[ids1, ids2[0]:ids2[-1]+1, ids3[0]:ids3[-1]+1]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (int, np.int64)) and \
            isinstance(ids3, (np.ndarray, list)):
        return x[ids1[0]:ids1[-1]+1, ids2, ids3[0]:ids3[-1]+1]

    if isinstance(ids1, (np.ndarray, list)) and \
            isinstance(ids2, (np.ndarray, list)) and \
            isinstance(ids3, (np.ndarray, list)):
        return x[ids1[0]:ids1[-1]+1, ids2[0]:ids2[-1]+1, ids3[0]:ids3[-1]+1]


def test_array(x, name_axis, format_axis=':0.3f'):
    form = '{' + format_axis + '}'
    line_x1 = form.format(x[0])
    line_x2 = form.format(x[-1])
    desc_line = name_axis + ' = [' + line_x1 + ', ' + line_x2 + ']'
    return desc_line


def get_interval(x, x_intervals, name_x, format_x):
    nx_intervals = np.shape(x_intervals)[0]
    ids_x_intervals, x_final_intervals, lines_x = [], [], []
    for count_x in range(nx_intervals):
        if x_intervals[count_x] is None:
            x_1, x_2 = x[0], x[-1]
        else:
            x_1, x_2 = x_intervals[count_x][0], x_intervals[count_x][-1]
        one_x, one_ids_x = get_array(x, x_1, x_2)
        x_final_intervals.append(one_x)
        ids_x_intervals.append(one_ids_x)

        line_format = '{:' + format_x + '}'
        line_format_x = '[' + line_format + ',' + line_format + ']'
        line_x = name_x + ' = ' + line_format_x.format(one_x[0], one_x[-1])
        lines_x.append(line_x)
    return ids_x_intervals, x_final_intervals, lines_x


# averaging in space:
def find_avs(oo):
    # for several signals oo.vars
    # every signal is averaged in several s-intervals
    # s-intervals are the same for every signal

    vvars = oo.get('vars', [])  # signals (t, s) to average
    ts = oo.get('ts', [])  # time grids
    ss = oo.get('ss', [])  # space grids
    ns_av = oo.get('ns_av', 0)  # number of intervals to average in space
    label_s = oo.get('label_s', 's')  # describe space coordinate

    n_vars = len(vvars)  # number of signals

    opts_av = oo.get('opts_av', [])  # mean, rms
    def_opt_av = 'mean'
    if len(opts_av) == 0:
        for ivar in range(n_vars):
            opts_av.append(def_opt_av)

    # intervals
    vars_avs = []
    line_av = ''
    for ivar in range(n_vars):
        var, t, s = vvars[ivar], ts[ivar], ss[ivar]
        t_int, ids_t_int = get_array_oo(oo, t, 't')
        opt_av = opts_av[ivar]

        data_s_av, lines_s_av = {}, {}
        for is_av in range(ns_av):
            s_av, ids_s_av = get_array_oo(oo, s, 's_av{:d}'.format(is_av + 1))
            temp = get_slice(var, ids_t_int, ids_s_av)
            if opt_av == 'mean':
                data_s_av[is_av]  = np.mean(temp, axis=1)
                line_av = 'mean_s:\ '
            if opt_av == 'rms':
                data_s_av[is_av]  = np.sqrt(np.mean(temp**2, axis=1))
                line_av = 'rms_s:\ '
            line_av = line_av + test_array(s_av, label_s, ':0.3f')
            lines_s_av[is_av] = line_av

        res = {'data': data_s_av, 't': t_int, 'lines_avs': lines_s_av}
        vars_avs.append(res)
    return vars_avs


# averaging in time:
def find_avt(oo):
    vvars = oo.get('vars', [])  # signals (t, s) to average
    ts = oo.get('ts', [])  # time grids
    ss = oo.get('ss', [])  # space grids
    nt_av = oo.get('nt_av', 0)  # number of intervals to average in time
    label_t = oo.get('label_t', 't[wci^{-1}]')  # describe time normalization

    n_vars = len(vvars)  # number of signals

    opts_av = oo.get('opts_av', [])  # mean, rms
    def_opt_av = 'mean'
    if len(opts_av) == 0:
        for ivar in range(n_vars):
            opts_av.append(def_opt_av)

    # intervals
    vars_avs = []
    data_av, line_av = np.nan, ''
    for ivar in range(n_vars):
        var, t, s = vvars[ivar], ts[ivar], ss[ivar]
        s_int, ids_s_int = get_array_oo(oo, s, 's')
        opt_av = opts_av[ivar]

        data_t_av, lines_t_av = {}, {}
        for it_av in range(nt_av):
            t_av, ids_t_av = get_array_oo(oo, t, 't_av{:d}'.format(it_av + 1))
            temp = get_slice(var, ids_t_av, ids_s_int)
            if opt_av == 'mean':
                data_av  = np.mean(temp, axis=0)
                line_av = 'mean_t:\ '
            if opt_av == 'rms':
                data_av  = np.sqrt(np.mean(temp**2, axis=0))
                line_av = 'rms_t:\ '
            line_av = line_av + test_array(t_av, label_t, ':0.2e')
            lines_t_av[it_av] = line_av

            # ## FOR NORMALIZATION
            # data_av[np.isnan(data_av)] = -np.inf
            # data_av[np.isinf(data_av)] = -np.inf
            # data_av = data_av / np.max(np.abs(data_av))
            # data_av[np.isinf(data_av)] = np.nan

            data_t_av[it_av] = data_av

        res = {'data': data_t_av, 's': s_int, 'lines_avt': lines_t_av}
        vars_avs.append(res)
    return vars_avs


# NEW: get indices of a corresponding domain:
def get_ids(x, x_domain, format_x='{:0.3e}'):
    # x = [x1, x2, x3, ...], where must be x[i] < x[i+1]
    # x_domain = some value or an array from two numbers
    x_domain = np.array([x_domain])
    if len(np.shape(x_domain)) == 2:
        x_domain = x_domain[0]

    id_start = np.where(x >= x_domain[0])[0]
    if id_start.size != 0:
        id_start = id_start[0]
    else:
        id_start = len(x) - 1

    id_end = np.where(x <= x_domain[-1])[0]
    if id_end.size != 0:
        id_end = id_end[-1]
    else:
        id_end = 0

    if id_end < id_start:
        id_end = id_start

    ids = [i for i in range(id_start, id_end + 1)]
    x_res = np.array(x[ids[0]:ids[-1] + 1])

    if len(ids) == 1:
        x_res = x_res[0]
        ids = ids[0]
        line_x = format_x.format(x_res)
    else:
        line_temp = '[' + format_x + ', ' + format_x + ']'
        line_x = line_temp.format(x_res[0], x_res[-1])

    return ids, x_res, line_x


# Check if all list-elements are unique in a list y
def is_unique(y):
    # y - list: [[...], [...], [...], ...]
    ny = len(y)
    for id_y in range(ny):
        if y[id_y] in y[id_y+1:ny]:
            return False
    return True


# Create an array with different time intervals
def get_t_intervals(oo, flag_print):
    # ---------------------------------------------------------
    # Create an array with several random time intervals within
    # a particular working time domain. Every time interval
    # is equal to an integer number of the GAM intervals.
    # ---------------------------------------------------------
    # -> nsamples - number of time intervals to choose
    # -> t_work - working time domain
    # -> min_n_periods - minimum number of a period that
    #   should be inside of every time interval
    # -> t_period - length of the period

    # parameters:
    nsamples = oo.get('nsamples', None)
    t_work = oo.get('t_work', None)

    min_n_periods = oo.get('min_n_periods', None)
    t_period = oo.get('t_period', None)

    # maximum begin time point
    id_max_start_point, _, _ = get_ids(
        t_work, t_work[-1] - min_n_periods * t_period
    )
    id_max_start_point -= 1
    if id_max_start_point == 0:
        print('ERROR: work time interval is too narrow.')
        return None

    # array with random begin time points
    ids_points_begin = np.random.randint(
        id_max_start_point + 1, size=nsamples
    ).astype('uint64')

    # For every random start point, define a length of a time interval
    ids_chosen_points = [[None, None]]
    for id_point_begin in ids_points_begin:

        chosen_comb = [None, None]
        count_n_begin_points = -1
        while chosen_comb in ids_chosen_points:

            # check if number of samples is too high
            count_n_begin_points += 1
            if count_n_begin_points > len(t_work):
                print('Error: number of samples is too high, '
                      'or work time domain is too narrow')
                return None

            # number of the GAM periods inside of the domain
            # from the current start point till
            # the right boundary of the working time domain
            max_n_periods = np.int(
                (t_work[-1] - t_work[id_point_begin]) / t_period
            )

            # array of available number of GAM periods
            possible_n_periods = \
                np.array([n_one for n_one in
                          range(min_n_periods, max_n_periods + 1)
                          ])

            rand_n_period = None
            while chosen_comb in ids_chosen_points:
                # remove already used number of periods
                possible_n_periods = \
                    possible_n_periods[(possible_n_periods != rand_n_period)]

                # if there are not more options of period numbers, then
                # change a begin time point
                if len(possible_n_periods) is 0:
                    id_point_begin_new = id_point_begin
                    while id_point_begin_new == id_point_begin:
                        id_point_begin_new = np.random.randint(id_max_start_point + 1)
                    id_point_begin = id_point_begin_new
                    break

                # choose a length of a time interval
                rand_n_period = np.random.choice(possible_n_periods, replace=True)

                # end time point
                id_point_end, _, _ = get_ids(
                    t_work,
                    t_work[id_point_begin] + rand_n_period * t_period
                )
                chosen_comb = [id_point_begin, id_point_end]

        ids_chosen_points.append(chosen_comb)

    del ids_chosen_points[0]

    res_t_intervals = []
    for one_ids_time_interval in ids_chosen_points:
        # noinspection PyTypeChecker
        res_t_intervals.append(
            t_work[one_ids_time_interval[0]:one_ids_time_interval[-1] + 1]
        )

    if is_unique(ids_chosen_points) and flag_print:
        print('All chosen time intervals are unique.')

    res = {
        't_intervals': res_t_intervals,
        'ids_intervals': ids_chosen_points,
    }

    return res


def get_t_intervals_old(nsamples, t_work, min_n_gam_periods, gam_t_period):
    # ---------------------------------------------------------
    # Create an array with several random time intervals within
    # a particular working time domain. Every time interval
    # is equal to an integer number of the GAM intervals.
    # ---------------------------------------------------------
    # -> nsamples - number of time intervals to choose
    # -> t_work - working time domain
    # -> min_n_gam_periods - minimum number of the GAM periods that
    #   should be inside of every time interval
    # -> gam_t_period - GAM period

    # maximum begin time point
    id_max_start_point, _, _ = get_ids(
        t_work, t_work[-1] - min_n_gam_periods * gam_t_period
    )
    id_max_start_point -= 1
    if id_max_start_point == 0:
        print('ERROR: work time interval is too narrow.')
        return None

    # array with random begin time points
    ids_points_begin = np.random.randint(
        id_max_start_point + 1, size=nsamples
    ).astype('uint64')

    # For every random start point, define a length of a time interval
    ids_chosen_points = [[None, None]]
    for id_point_begin in ids_points_begin:

        chosen_comb = [None, None]
        count_n_begin_points = -1
        while chosen_comb in ids_chosen_points:

            # check if number of samples is too high
            count_n_begin_points += 1
            if count_n_begin_points > len(t_work):
                print('Error: number of samples is too high, '
                      'or work time domain is too narrow')
                return None

            # number of the GAM periods inside of the domain
            # from the current start point till
            # the right boundary of the working time domain
            max_n_gam_periods = np.int(
                (t_work[-1] - t_work[id_point_begin]) / gam_t_period
            )

            # array of available number of GAM periods
            possible_n_gam_periods = \
                np.array([n_one for n_one in
                          range(min_n_gam_periods, max_n_gam_periods + 1)
                          ])

            rand_n_period = None
            while chosen_comb in ids_chosen_points:
                # remove already used number of periods
                possible_n_gam_periods = \
                    possible_n_gam_periods[(possible_n_gam_periods != rand_n_period)]

                # if there are not more options of period numbers, then
                # change a begin time point
                if len(possible_n_gam_periods) is 0:
                    id_point_begin_new = id_point_begin
                    while id_point_begin_new == id_point_begin:
                        id_point_begin_new = np.random.randint(id_max_start_point + 1)
                    id_point_begin = id_point_begin_new
                    break

                # choose a length of a time interval
                rand_n_period = np.random.choice(possible_n_gam_periods, replace=True)

                # end time point
                id_point_end, _, _ = get_ids(
                    t_work,
                    t_work[id_point_begin] + rand_n_period * gam_t_period
                )
                chosen_comb = [id_point_begin, id_point_end]

        ids_chosen_points.append(chosen_comb)

    del ids_chosen_points[0]

    res_t_intervals = []
    for one_ids_time_interval in ids_chosen_points:
        # noinspection PyTypeChecker
        res_t_intervals.append(
            t_work[one_ids_time_interval[0]:one_ids_time_interval[-1] + 1]
        )

    if is_unique(ids_chosen_points):
        print('All chosen time intervals are unique.')

    res = {
        't_intervals': res_t_intervals,
        'ids_intervals': ids_chosen_points,
    }

    return res


# Should be recheck before use it
def get_t_intervals_zeros(nsamples, t_work, min_n_gam_periods, gam_t_period, t_je_abs_peaks):
    # time array where areas around J*E = 0 are excluded
    je_quarter_period = gam_t_period / 8
    part_peak = 0.0
    ids_t_work_without_zero = []
    for t_one_peak in t_je_abs_peaks:
        id_left, t_left, _ = get_ids(
            t_work, t_one_peak - part_peak * je_quarter_period
        )
        id_right, t_right, _ = get_ids(
            t_work, t_one_peak + part_peak * je_quarter_period
        )
        ids_current = [i for i in range(id_left, id_right+1)]
        ids_t_work_without_zero += ids_current
    t_work_woZeros = t_work[ids_t_work_without_zero]

    ids_t_work_without_zero = np.array(ids_t_work_without_zero)

    # maximum begin time point
    id_max_start_point, _, _ = get_ids(
        t_work, t_work_woZeros[-1] - min_n_gam_periods * gam_t_period
    )
    id_max_start_point -= 1
    if id_max_start_point == 0:
        print('ERROR: work time interval is too narrow.')
        return None

    # time interval without zeros excluding very right area of the several gam periods
    ids_t_work_without_zero_max_begin = ids_t_work_without_zero[
        (ids_t_work_without_zero <= id_max_start_point)
    ]

    # random begin time points
    ids_points_begin = np.random.choice(
        ids_t_work_without_zero_max_begin, size=nsamples, replace=True
    )  # check size of ids_t_work_without_zero_max_begin after choice

    ids_chosen_points = [[None, None]]
    for id_point_begin in ids_points_begin:

        chosen_comb = [None, None]
        # count_n_begin_points = -1
        while chosen_comb in ids_chosen_points:

            # # check if number of samples is too high
            # count_n_begin_points += 1
            # if count_n_begin_points > len(t_work_woZeros):
            #     print('Error: number of samples is too high, '
            #           'or work time domain is too narrow')
            #     break

            # maximum J*E periods from the current start point till
            # the right boundary of the working time domain
            max_n_gam_periods = np.int(
                (t_work_woZeros[-1] - t_work[id_point_begin]) / gam_t_period
            )

            # array of available number of J*E periods
            possible_n_gam_periods = \
                np.array([n_one for n_one in range(min_n_gam_periods, max_n_gam_periods+1)])

            rand_n_period = None
            while chosen_comb in ids_chosen_points:
                # remove already used number of periods
                possible_n_gam_periods = \
                    possible_n_gam_periods[(possible_n_gam_periods != rand_n_period)]

                # if there are not more options of period numbers, then
                # change a begin time point
                if len(possible_n_gam_periods) is 0:
                    id_point_begin_new = id_point_begin
                    while id_point_begin_new == id_point_begin:
                        id_point_begin_new = np.random.choice(
                            ids_t_work_without_zero_max_begin, replace=True
                        )
                    id_point_begin = id_point_begin_new
                    break

                rand_n_period = np.random.choice(possible_n_gam_periods, replace=True)

                # end time point
                id_point_end, _, _ = get_ids(
                    t_work,
                    t_work[id_point_begin] + rand_n_period * gam_t_period
                )
                chosen_comb = [id_point_begin, id_point_end]

        ids_chosen_points.append(chosen_comb)

    del ids_chosen_points[0]

    res_t_intervals = []
    for one_ids_time_interval in ids_chosen_points:
        # noinspection PyTypeChecker
        res_t_intervals.append(
            t_work[one_ids_time_interval[0]:one_ids_time_interval[-1]+1]
        )

    if is_unique(ids_chosen_points):
        print('All chosen time intervals are unique.')

    res = {
        't_intervals':   res_t_intervals,
        'ids_intervals': ids_chosen_points,
        'ids_wo_zeros': ids_t_work_without_zero,
        'ids_wo_zeros_max_begin': ids_t_work_without_zero_max_begin,
    }

    return res







