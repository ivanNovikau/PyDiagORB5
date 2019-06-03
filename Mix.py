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
        line_x = format_x.format(x_res)
    else:
        line_temp = '[' + format_x + ', ' + format_x + ']'
        line_x = line_temp.format(x_res[0], x_res[-1])

    return ids, x_res, line_x









