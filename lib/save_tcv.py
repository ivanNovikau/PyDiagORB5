import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import zf_gam as zf
import ITG_gamma as itg
import write_data as wr
import common
import numpy as np
from scipy import interpolate
import h5py as h5


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(zf)
    mix.reload_module(itg)
    mix.reload_module(wr)
    mix.reload_module(common)


def save_signals(dd):
    # file to write to:
    out_file = 'data_tcv.h5'
    path_to_write = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/TCV/NL/muT12-mun5/n80-potsc/'
    path_to_write += out_file

    # radial point
    s1 = 0.953

    # get 1d signals to save:
    oo_t_s1 = {
        'ovars': [
                    ['zonal', 'erbar'],
                    ['nonzonal', 'er_r', 0.0],
                 ],
        'avrs': [
                    ['ts', 'point-s', [s1]],
                    ['ts', 'point-s', [s1]],
                ],
        'dds': [dd] * 2,
    }
    vvars = common.choose_vars(oo_t_s1)
    n_vars = len(vvars)

    # data and grids
    t     = vvars[0]['x']
    erbar = vvars[0]['data'][0]
    er    = np.interp(t, vvars[1]['x'], vvars[1]['data'][0])

    wc = dd['wc']

    # save data:
    wr.create_file(path_to_write)
    wr.save_data(path_to_write, 's_point', s1,
                 desc=u'radial point, where the signals are taken')
    wr.save_data(path_to_write, 'chi_point', 0,
                 desc=u'poloidal point, where Er is taken')
    wr.save_data(path_to_write, 'wc', wc, desc='cyclotron frequency')
    wr.save_data(path_to_write, 't', t, desc='time grid, [wc^{-1}]')
    wr.save_data(path_to_write, 'Er_zonal', erbar, desc='zonal radial electric field (GAM) at s_point')
    wr.save_data(path_to_write, 'Er',          er, desc='full radial electric field at s_point and chi_point')

    f = h5.File(path_to_write, 'r')
    t = np.array(f['t'])
    erbar = np.array(f['Er_zonal'])
    er = np.array(f['Er'])
    f.close()

    curves = crv.Curves().xlab('t')
    curves.flag_semilogy = True
    curves.new().XS(t).YS(erbar)
    cpr.plot_curves(curves)

    curves = crv.Curves().xlab('t')
    curves.flag_semilogy = True
    curves.new().XS(t).YS(er)
    cpr.plot_curves(curves)
