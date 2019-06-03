import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import zf_gam as zf
import general as gn
import write_data as wr
import transport
from scipy.fftpack import next_fast_len
from polycoherence import _plot_signal, polycoherence, plot_polycoherence
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
    mix.reload_module(zf)
    mix.reload_module(gn)
    mix.reload_module(wr)
    mix.reload_module(transport)


# NEW: take signal (vpar,mu)
def choose_one_var_vparmu(ovar, dd):
    opt_var      = ovar[0]
    species_name = ovar[1]

    vvar, vpar, mu, tit_var = None, None, None, ''
    if opt_var == 'jdote_es-mean-t':
        t_int = ovar[2]

        jdote_es_dict = dd[species_name].jdote_es(dd)
        t, vpar, mu = jdote_es_dict['t'], jdote_es_dict['vpar'], jdote_es_dict['mu']
        ids_t, t_int, line_t = mix.get_ids(t, t_int)

        vvar = np.mean(jdote_es_dict['data'][ids_t, :, :], axis=0)

        tit_var  = species_name + ':\ <J*E>_{t}'
        tit_var = 't = ' + line_t + ':\ ' + tit_var

    res = {
        'data': vvar,
        'vpar': vpar,
        'mu': mu,
        'tit': tit_var
    }

    return res