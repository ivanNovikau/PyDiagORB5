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
import gam_theory
import gam_exp
import Geom as geom
import common

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
    mix.reload_module(gam_theory)
    mix.reload_module(gam_exp)
    mix.reload_module(geom)
    mix.reload_module(common)


def GENE_ES_global(oo):
    dd = oo.get('dds', None)[0]
    flag_semilogy = oo.get('flag_semilogy', False)

    # --- GENE data ---
    path_gene = 'd:/Work-Projects/MyProgs/ORB_data/turbulence/ITPA/GENE/ES-global/'

    x_gene = np.loadtxt(path_gene + '/x.dat')
    t_gene = np.loadtxt(path_gene + '/t.dat')  # in a/cs units
    phi_gene = np.loadtxt(path_gene + '/phi.dat')  # (x,t), in T/(e*rho_star) units

    # - renormalize GENE data -
    rho_star = 2 / dd['Lx']
    inv_rho_star = 1. / rho_star
    renorm_phi = dd['wc'] * dd['R0'] / dd['cs']

    curves = crv.Curves().xlab('x').ylab('t').tit('Phi')
    curves.new().XS(x_gene).YS(t_gene * renorm_phi).ZS(phi_gene).lev(60)
    cpr.plot_curves_3d(curves)

    max_phi_gene = np.max(phi_gene, axis=0)  # -> function on t

    # t_gene       = t_gene * inv_rho_star  # in 1/wc units
    t_gene = t_gene * renorm_phi  # in 1/wc units
    max_phi_gene = max_phi_gene / renorm_phi  # in T/e units

    # --- ORB5 DATA ---
    vvar = common.choose_vars(oo)[0]
    data_orb = vvar['data'][0]
    t_orb = vvar['x']
    leg_orb = vvar['legs'][0]

    curves = crv.Curves().xlab('t').tit('ORB5:\ ' + leg_orb)
    curves.new().XS(t_orb).YS(data_orb)
    cpr.plot_curves(curves)

    # # --- COMPARE ---
    curves = crv.Curves().xlab('t').tit('max_s:\ \Phi')
    curves.flag_semilogy = flag_semilogy
    curves.new().XS(t_gene).YS(max_phi_gene).leg('GENE')
    curves.new().XS(t_orb).YS(data_orb).leg('ORB5')
    cpr.plot_curves(curves)



































