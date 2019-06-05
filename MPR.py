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
import scipy.io
from scipy import constants


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
        'data': vvar,  # (t, mu, vpar)
        'vpar': vpar,
        'mu': mu,
        'tit': tit_var
    }

    return res


def data_GENE_CPS2019():
    path_data = 'd:/Work-Projects/MyProgs/ORB_data/MPR/nled/kin/linear/GENE-CPS2019-DF/data.mat'
    dGENE = scipy.io.loadmat(path_data)

    Te_div_B0 = dGENE['Te_B0'][0, 0]
    cs = dGENE['cs'][0, 0]
    Lref = dGENE['lref'][0, 0]
    vpar = dGENE['out_Ivan_v'][0]
    mu   = dGENE['out_Ivan_mu'][0]
    data_t_jet_jed_jef = dGENE['out_Ivan']

    t    = data_t_jet_jed_jef[:, 0]
    Ptot = data_t_jet_jed_jef[:, 1]  # total J*E
    Pth  = data_t_jet_jed_jef[:, 2]  # J*E of thermal species
    Pf   = data_t_jet_jed_jef[:, 3]  # J*E of fast species

    gamma_vmu = dGENE['out_Ivan_gamma_vm']  # (vpar, mu)

    res = {
        'Te_div_B0': Te_div_B0, 'cs': cs, 'Lref': Lref,
        'vpar': vpar, 'mu': mu, 't': t,
        'Ptot': Ptot, 'Pth': Pth, 'Pf': Pf,
        'gamma_vmu': gamma_vmu
    }

    return res


# Comparison with GENE for the paper CPS-2019:
def comparison_with_GENE_CPS2019(dd, oo):
    # Parameters
    dd_f0   = oo.get('dd_f0', None)

    labx_ms = 'v_{\parallel}(m/s)'
    labx_norm = 'v_{\parallel}'

    leg_gene = 'GENE:\ ' + '-\langle \mathcal{P} \\rangle_{\mu}'
    leg_orb  = 'ORB5:\ ' + '\langle \mathcal{P} \\rangle_{\mu}'

    # Find iniitial distribution functions of fast species
    rd.distribution_1d(dd_f0, 'fast')

    # velocity normalization:
    Te_max = np.max(dd['electrons'].nT_equil['T_J'])
    cs_max = np.sqrt(Te_max / dd['pf'].mass)
    norm_v = 1. / cs_max

    # --- GENE data ---
    gd = data_GENE_CPS2019()
    gvpar_norm = gd['vpar'] * norm_v
    gamma_vmu_aver_mu = np.mean(gd['gamma_vmu'], axis=1)
    gamma_vmu_aver_mu = - gamma_vmu_aver_mu

    # --- ORB5 data ---
    # Initial distribution function:
    f0       = dd_f0['fast'].f_1d['f_vel_1d'][0]
    ovpar_f0 = dd_f0['fast'].f_1d['vpar']

    # <J*E>_mu
    t_point_csR0 = 73.5  # time moment where GENE gamma is taken
    renorm_t_coef = dd['wc'] / (gd['cs']/gd['Lref'])
    t_point = t_point_csR0 * renorm_t_coef

    jdote_es_dict = dd['fast'].jdote_es(dd)
    id_t, _, line_t = mix.get_ids(jdote_es_dict['t'], t_point)

    oje_t1 = np.squeeze(jdote_es_dict['data'][id_t, :, :])
    oje_aver_mu = np.mean(oje_t1, axis=0)
    ovpar = jdote_es_dict['vpar']

    # --- Plot <J*E>_mu for fast species ---
    curves = crv.Curves().xlab(labx_norm).ylab('norm.\ values')\
        .tit('fast\ deuterium')
    curves.flag_norm = True
    curves.new()\
        .XS(gvpar_norm)\
        .YS(gamma_vmu_aver_mu)\
        .leg(leg_gene)
    curves.new() \
        .XS(ovpar) \
        .YS(oje_aver_mu) \
        .leg(leg_orb)
    curves.new() \
        .XS(ovpar_f0) \
        .YS(f0)\
        .sty('--').col('grey').leg('f_0')
    cpr.plot_curves(curves)
