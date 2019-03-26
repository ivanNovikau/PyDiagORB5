import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)


def plot_t(dd, s1, s2, chi1, oo={}):
    rd.phibar(dd)
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    _, ids_s = mix.get_array(s, s1, s2)
    t, ids_t = mix.get_array_oo(oo, t, 't')
    id_chi, _ = mix.find(chi, chi1)

    # non-zonal (at phi = 0) Phi in chosen intervals
    pot_nz_chi = mix.get_slice(dd['potsc']['data'], ids_t, id_chi, ids_s)

    # averaging of the Phi on s
    pot_nz = np.mean(pot_nz_chi, axis=1)

    # plot data:
    cpr.plot_x1(t, pot_nz)


def plot_schi(dd, t1, oo={}):
    rd.phibar(dd)
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    s, ids_s     = mix.get_array_oo(oo, s, 's')
    chi, ids_chi = mix.get_array_oo(oo, chi, 'chi')
    id_t1, _ = mix.find(t, t1)

    # signal in the chosen intervals
    pot_nz = mix.get_slice(dd['potsc']['data'], id_t1, ids_chi, ids_s)

    # form 3d curve:
    curves = crv.Curves().xlab('s').ylab('\chi').tit('\Phi')
    curves.new('Phi_schi').XS(s).YS(chi).ZS(pot_nz) \
        .leg('\Phi').cmp('hot')
    cpr.plot_curves_3d(curves, 'mat')


def plot_rz(dd, t1, oo={}):
    rd.phibar(dd)
    t = dd['potsc']['t']  # create a new reference
    r = dd['potsc']['r']  # create a new reference
    z = dd['potsc']['z']  # create a new reference

    # intervals
    id_t1, _ = mix.find(t, t1)

    # signal in the chosen intervals
    pot_nz = np.array(dd['potsc']['data'][id_t1, :, :])  # actually copy data

    # form 3d curve:
    curves = crv.Curves().xlab('r').ylab('z').tit('\Phi')
    curves.new('Phi_schi').XS(r).YS(z).ZS(pot_nz) \
        .leg('\Phi').cmp('jet')
    cpr.plot_curves_3d(curves, 'mat')


def calc_gamma_chi0(dd, s1, s2, chi1, oo={}):
    rd.phibar(dd)
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    _, ids_s = mix.get_array(s, s1, s2)
    t, ids_t = mix.get_array_oo(oo, t, 't')
    id_chi, _ = mix.find(chi, chi1)

    # non-zonal (at phi = 0) Phi in chosen intervals
    pot_nz_chi = mix.get_slice(dd['potsc']['data'], ids_t, id_chi, ids_s)

    # averaging of the Phi on s
    pot_nz = np.mean(pot_nz_chi, axis=1)

    # estimation:
    wg_est = ymath.estimate_wg(t, pot_nz)

    # plotting: estimation: time evolution and peaks
    curves_est = crv.Curves().xlab('t').ylab('\Phi')
    curves_est.flag_semilogy = True
    curves_est.new('init')\
        .XS(t).YS(pot_nz).leg('init')
    curves_est.new('peaks')\
        .XS(wg_est['x_peaks']).YS(wg_est['y_peaks'])\
        .leg('peaks').sty('o').col('green')
    curves_est.new('fitting')\
        .XS(t).YS(wg_est['y_fit'])\
        .leg('fitting').col('red').sty('--')
    cpr.plot_curves(curves_est)

    # print results of estimation
    print('--- Estimation ---')
    print('w = {:0.3e}'.format(wg_est['w']))
    print('g = {:0.3e}'.format(wg_est['g']))

    # advanced w,g calculation
    ainf = {'est': wg_est,
            'x_start': wg_est['x_peaks'][0],
            'x_end': wg_est['x_peaks'][-1]
            }
    wg_adv = ymath.advanced_wg(t, pot_nz, ainf)
    curves_adv = crv.Curves().xlab('t').ylab('\Phi')
    curves_adv.flag_semilogy = True
    curves_adv.new('init') \
        .XS(t).YS(pot_nz).leg('init')
    curves_adv.new('fitting') \
        .XS(wg_adv['x_fit']).YS(wg_adv['y_fit'])\
        .leg('adv. fitting').col('red').sty('--')
    cpr.plot_curves(curves_adv)

    # print results of the advanced plotting
    print('--- Advanced ---')
    print('w = {:0.3e}'.format(wg_adv['w']))
    print('g = {:0.3e}'.format(wg_adv['g']))


def calc_gamma_chimax(dd, s1, s2, oo={}):
    rd.phibar(dd)
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    _, ids_s = mix.get_array(s, s1, s2)
    t, ids_t = mix.get_array_oo(oo, t, 't')
    chi, ids_chi = mix.get_array_oo(oo, chi, 'chi')

    # non-zonal Phi in chosen intervals
    pot_chi = mix.get_slice(dd['potsc']['data'], ids_t, ids_chi, ids_s)

    # averaging on s
    pot_aver_s = np.mean(pot_chi, axis=2)

    # maximums along chi:
    pot_max = np.amax(np.abs(pot_aver_s), axis=1)

    # gamma calculation
    wg_est = ymath.estimate_g(t, pot_max)
    wg_adv = ymath.advanced_g(t, pot_max, {'est': wg_est})

    # plotting: estimation: time evolution and peaks
    curves = crv.Curves().xlab('t').ylab('\Phi')\
        .tit('Estimation:\                g = {:0.3e}'.format(wg_est['g']))\
        .tit('\\\\Full\ signal\ fitting:\ g = {:0.3e}'.format(wg_adv['g']))
    curves.flag_semilogy = True
    curves.new('max') \
        .XS(t).YS(pot_max).leg('max.\ along\ \chi')
    if 'x_peaks' in wg_est:
        curves.new('peaks') \
            .XS(wg_est['x_peaks']).YS(wg_est['y_peaks']) \
            .leg('peaks').sty('o').col('grey')
    curves.new('fitting') \
        .XS(t).YS(wg_est['y_fit']) \
        .leg('estimation').col('red').sty('--')
    curves.new('fitting') \
        .XS(wg_adv['x_fit']).YS(wg_adv['y_fit']) \
        .leg('full\ signal\ fitting').col('green').sty(':')
    cpr.plot_curves(curves)

    # print results (keep it to have from screen)
    print('--- Estimation ---')
    print('g = {:0.3e}'.format(wg_est['g']))
    print('--- Advanced ---')
    print('g = {:0.3e}'.format(wg_adv['g']))


def plot_schi_max_along_chi(dd, t1, oo={}):
    # plot max of phi on (s, chi) at a particular time step
    rd.phibar(dd)
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    s, ids_s     = mix.get_array_oo(oo, s, 's')
    chi, ids_chi = mix.get_array_oo(oo, chi, 'chi')
    id_t1, _ = mix.find(t, t1)

    # non-zonal Phi in chosen intervals
    pot_chi = mix.get_slice(dd['potsc']['data'], id_t1, ids_chi, ids_s)

    # maximums along chi:
    pot_max = np.amax(np.abs(pot_chi), axis=0)

    # build array for max:
    chi_max = np.zeros(np.size(s))
    for i in range(np.size(pot_max)):
        max1 = pot_max[i]
        pot_s1 = pot_chi[:, i]
        id_chi1 = np.where(pot_s1 == max1)
        if len(id_chi1[0]) is 0:
            id_chi1 = np.where(pot_s1 == -max1)
        id_chi1 = id_chi1[0][0]
        chi_max[i] = chi[id_chi1]

    # plotting:
    curves = crv.Curves().xlab('s').ylab('\chi').tit('\Phi')
    curves.new('schi') \
        .XS(s).YS(chi).ZS(pot_chi)\
        .leg('\Phi').cmp('hot')
    curves.new('max').XS(s).YS(chi_max).leg('max')
    cpr.plot_curves_3d(curves)


def plot_rz_max_along_chi(dd, t1, oo={}):
    rd.phibar(dd)
    t = dd['potsc']['t']  # create a new reference
    r = dd['potsc']['r']  # create a new reference
    z = dd['potsc']['z']  # create a new reference

    # intervals
    id_t1, _ = mix.find(t, t1)

    # signal in the chosen intervals
    pot_nz = np.array(dd['potsc']['data'][id_t1, :, :])  # actually copy data

    # maximums along chi:
    pot_max = np.amax(np.abs(pot_nz), axis=0)

    # build array for max:
    # r_max = np.zeros(np.shape(pot_nz))
    # z_max = np.zeros(np.shape(pot_nz))
    r_max = np.zeros(np.size(pot_max))
    z_max = np.zeros(np.size(pot_max))
    for id_j in range(np.size(pot_max)):
        max1 = pot_max[id_j]
        pot_s1 = pot_nz[:, id_j]
        id_i = np.where(pot_s1 == max1)
        if len(id_i[0]) is 0:
            id_i = np.where(pot_s1 == -max1)
        id_i = id_i[0][0]
        r_max[id_j] = r[id_i, id_j]
        z_max[id_j] = z[id_i, id_j]

    # plotting:
    curves = crv.Curves().xlab('r').ylab('z').tit('\Phi')
    curves.new('rz') \
        .XS(r).YS(z).ZS(pot_nz)\
        .leg('\Phi').cmp('hot')
    curves.new('max').XS(r_max).YS(z_max).leg('max')
    cpr.plot_curves_3d(curves)


def anim_st_chi0(dd, chi0, oo={}):
    rd.phibar(dd)
    t = dd['potsc']['t']  # create a new reference
    s = dd['potsc']['s']  # create a new reference
    chi = dd['potsc']['chi']  # create a new reference

    # intervals
    s, ids_s = mix.get_array_oo(oo, s, 's')
    t, ids_t = mix.get_array_oo(oo, t, 't')
    id_chi1, _ = mix.find(chi, chi0)

    # signal in the chosen intervals
    pot_nz = mix.get_slice(dd['potsc']['data'], ids_t, id_chi1, ids_s)

    # form 2d curve:
    curves = crv.Curves().xlab('s').wlab('t').zlab('\Phi').tit('\Phi')
    curves.new('anim_Phi_st').XS(s).WS(t).ZS(pot_nz).leg('\Phi')
    curves.set_limits()
    curves.flag_norm = True

    # animation:
    cpr.animation_curves_2d(curves)




