import numpy as np
import importlib as imp
import Global_variables as GLO
import sys


def reload():
    # Important: put here all modules that you want to reload
    reload_module(GLO)
    return


def reload_module(obj_module):
    imp.reload(obj_module)
    obj_module.reload()


def normalization(sel_norm, dd=None):
    line_norm, coef_norm = '', 1

    if sel_norm == 't-ms':
        line_norm = '\ (ms)'
        if dd is not None:
            coef_norm = 1. / dd['wc'] * 1e3
    if sel_norm == 't-csr':
        line_norm = '\ [\omega_s^{-1}]'
        if dd is not None:
            coef_norm = (dd['cs']/dd['R0']) / dd['wc']
    if sel_norm == 'energy-transfer-W':
        line_norm = '\ [W]'
        if dd is not None:
            coef_norm = dd['T_speak'] * dd['wc']
    if sel_norm == 'energy-J':
        line_norm = '\ [J]'
        if dd is not None:
            coef_norm = dd['T_speak']
    if sel_norm == 'n-m3':
        line_norm = '\ [m^{-3}]'
        if dd is not None:
            coef_norm = dd['ele-nbar-m3']

    if sel_norm is not None:
        if len(sel_norm.split('-')) == 2:
            selector_type, sel_norm_new = sel_norm.split('-')
            if selector_type == 'frequency' or selector_type == 'gamma':
                coef_norm_w, coef_norm_g, line_norm_w, line_norm_g = \
                    choose_wg_normalization(sel_norm_new, dd)
                if selector_type == 'frequency':
                    line_norm, coef_norm = line_norm_w, coef_norm_w
                if selector_type == 'gamma':
                    line_norm, coef_norm = line_norm_g, coef_norm_g

    res_data = {
        'line_norm': line_norm,
        'coef_norm': coef_norm,
    }
    return res_data


def choose_wg_normalization(sel_norm, dd=None):
    coef_norm_w, coef_norm_g, line_norm_w, line_norm_g = \
        None, None, '', ''

    if sel_norm is None:
        sel_norm = 'wc'

    if sel_norm.lower() == 'wc':
        line_norm_w = line_norm_g = '\ [\omega_{ci}]'
        coef_norm_w = coef_norm_g = 1
    if sel_norm.lower() == 'vt':
        line_norm_w = line_norm_g = '\ [sqrt(2)*v_{th,i}/R_0]'
        if dd is not None:
            coef_norm_w = coef_norm_g = \
                dd['wc'] / (np.sqrt(2) * dd['vt'] / dd['R0'])
    if sel_norm.lower() == 'khz':
        line_norm_w = '\ [kHz]'
        line_norm_g = '\ [10^3\ s^{-1}]'
        if dd is not None:
            coef_norm_w = dd['wc'] / (1e3 * 2 * np.pi)
            coef_norm_g = dd['wc'] / 1e3
    if sel_norm == 'csa':
        line_norm_w = line_norm_g = '\ [c_s/a_0]'
        if dd is not None:
            coef_norm_w = coef_norm_g =\
                dd['wc'] / (dd['cs'] / dd['a0'])
    if sel_norm == 'csr':
        line_norm_w = line_norm_g = '\ [\omega_s]'  # which is cs/R0
        if dd is not None:
            coef_norm_w = coef_norm_g =\
                dd['wc'] / (dd['cs'] / dd['R0'])
    return coef_norm_w, coef_norm_g, line_norm_w, line_norm_g


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


def get_x_data_oo(oo, x, name_x, data):
    x_start = oo.get(name_x + '_start', x[0])
    x_end = oo.get(name_x + '_end', x[-1])
    return get_x_data_interval([x_start, x_end], x, data) # data_new, x_new


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


def get_id_int_strict(x, x_point):
    # x = [x1, x2, x3, ...], where must be x[i] < x[i+1]
    # x_point = some value
    # if x_point != x anywhere, return np.nan
    id = np.where(x == x_point)[0]
    if id.size != 0:
        id = id[0]
    else:
        id = np.nan

    return id


def get_x_data_interval(x_domain, x, data):
    ids_x, _, _ = get_ids(x, x_domain)
    ids_x = range(ids_x[0], ids_x[-1] + 1)
    return  data[ids_x], x[ids_x]


def get_data_interval_2d(data, x, x_domain, y, y_domain):
    ids_x, _, _ = get_ids(x, x_domain)
    ids_x = range(ids_x[0], ids_x[-1] + 1)
    ids_y, _, _ = get_ids(y, y_domain)
    ids_y = range(ids_y[0], ids_y[-1] + 1)
    data_res = data[ids_x[0]:ids_x[-1] + 1, ids_y[0]:ids_y[-1] + 1]
    x_res = x[ids_x]
    y_res = y[ids_y]
    return  data_res, x_res, y_res


def get_data_at_x1(x1, x, data):
    data_x1 = np.interp(x1, x, data)
    return data_x1


def create_consequent_time_intervals(oo_t):
    tmin = oo_t['tmin']
    tmax = oo_t['tmax']
    width_t = oo_t['width_t']
    step_t = oo_t['step_t']
    flag_inc_boundary = oo_t.get('flag_inc_boundary', False)
    half_width_t = width_t / 2

    t_right = tmin + width_t
    t_ints = [[tmin, t_right]]
    t_center = tmin + half_width_t
    while t_right < tmax:
        t_center += step_t

        t_left = t_center - half_width_t
        t_right = t_center + half_width_t
        if t_right > tmax:
            if not flag_inc_boundary:
                break
            else:
                t_right = tmax
        t_ints.append([t_left, t_right])

    return t_ints


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


# Create a line from a list of lines:
def create_line_from_list(list_lines):
    if list_lines is None:
        return None

    format_begin = r'\boldmath $'
    format_middle = '$\n \\boldmath $'
    format_end = '$'

    res_line = ''
    if isinstance(list_lines, list):
        for one_line in list_lines:
            if one_line is not None:
                res_line += one_line if res_line == '' else \
                    format_middle + one_line
    else:
        res_line = list_lines
    res_line = format_begin + res_line + format_end
    return res_line


# Print error message and exit from the program:
def error_mes(message):
    print('Error: ' + message)
    sys.exit(-1)


# translates an rgb tuple of int to a tkinter friendly color code
def to_rgb(rgb):
    return "#%02x%02x%02x" % rgb







