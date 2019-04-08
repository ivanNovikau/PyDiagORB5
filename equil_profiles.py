import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np
from scipy import interpolate


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)


def q_prof(dd):
    rd.q(dd)

    curves = crv.Curves().xlab('s').ylab('q').tit('safety\ factor\ profile')
    curves.new('q').XS(dd['q']['s']).YS(dd['q']['data']).leg('q')
    cpr.plot_curves(curves)


def nT_profs(dd):
    # output equilibrium temperature from ORB5 as it is
    curves_T = crv.Curves().xlab('s').ylab('T')\
        .tit('species\ temperature\ (ORB5\ output)').set_diff_styles()
    for sp_name in dd['species_names']:
        curves_T.new(sp_name)\
            .XS(dd[sp_name].nT_equil['s'])\
            .YS(dd[sp_name].nT_equil['T'])\
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # normalized temperature
    curves_T.new_tit('species\ norm.\ temperature').ylab('norm.\ T')
    curves_T.flag_norm = True
    cpr.plot_curves(curves_T)

    # temperature in keV
    curves_T = crv.Curves().xlab('s').ylab('T(keV)')\
        .tit('species\ temperature\ in\ keV').set_diff_styles()
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['T_keV']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # output equilibrium density from ORB5 as it is
    curves_n = crv.Curves().xlab('s').ylab('norm.\ n') \
        .tit('species\ norm.\ density').set_diff_styles()
    curves_n.flag_norm = True
    for sp_name in dd['species_names']:
        curves_n.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['n']) \
            .leg(sp_name)
    cpr.plot_curves(curves_n)

    # normalized temperature gradient
    curves_T = crv.Curves().xlab('s').ylab('norm.\ \\nabla(T)') \
        .tit('species\ norm.\ temperature\ gradient')\
        .set_diff_styles()
    curves_T.flag_norm = True
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['gradT']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # normalized density gradient
    curves_T = crv.Curves().xlab('s').ylab('norm.\ \\nabla(n)') \
        .tit('species\ norm.\ density\ gradient') \
        .set_diff_styles()
    curves_T.flag_norm = True
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['gradn']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # normalized log temperature gradient
    curves_T = crv.Curves().xlab('s').ylab('norm.\ \\nabla(log(T))') \
        .tit('species\ norm.\ log.\ temperature\ gradient') \
        .set_diff_styles()
    curves_T.flag_norm = True
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['grad_logT']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)

    # normalized log density gradient
    curves_T = crv.Curves().xlab('s').ylab('norm.\ \\nabla(log(n))') \
        .tit('species\ norm.\ log.\ density\ gradient') \
        .set_diff_styles()
    curves_T.flag_norm = True
    for sp_name in dd['species_names']:
        curves_T.new(sp_name) \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['grad_logn']) \
            .leg(sp_name)
    cpr.plot_curves(curves_T)
