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
    return out
