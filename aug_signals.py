import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import zf_gam as zf
import ITG_gamma as itg
import write_data as wr
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


def save_aug_signals(dd, oo):
    # path to a file, where signals are to be written to
    path_to_write = dd['path_to_write']

    # radial points, where the signals will be considered
    s_points  = oo.get('s_points', [0.5])
    ns_points = np.size(s_points)

    # --- zonal Er, Phi ---
    zf.erbar(dd)
    phibar = dd['phibar']['data']
    erbar = dd['erbar']['data']

    # space, time, frequency grids
    s = dd['erbar']['s']
    t = dd['erbar']['t']
    nt = np.size(t)

    ffres = ymath.fft_y(t)
    w2 = ffres['w2']
    nw2 = np.size(w2)

    # time and frequency normalization
    coef_norm_t = 1 / dd['wc']
    coef_norm_w = dd['wc'] / 1e3

    # poloidal angles, where a full signal will be considered
    chi_points = oo.get('chi_points', [0.0])
    nchi_points = np.size(chi_points)

    # --- full Phi and Er ---
    oo_phi = {'chi_s': chi_points}
    names_potsc = rd.potsc_chi(dd, oo_phi)
    names_ersc  = itg.ersc(dd, oo_phi)

    # --- Non-zonal Phi and Er
    oo_phi = {'chi_s': chi_points}
    names_phinz = itg.phinz(dd, oo_phi)
    names_ernz  = itg.ernz(dd, oo_phi)

    chi_points_res = []
    for one_name in names_ersc:
        chi_points_res.append(dd[one_name]['chi_1'])

    s_points_res = np.zeros(ns_points)
    phibar_t = np.zeros([ns_points, nt])
    phibar_w = np.zeros([ns_points, nw2])
    erbar_t = np.zeros([ns_points, nt])
    erbar_w = np.zeros([ns_points, nw2])
    phi_t = np.zeros([ns_points, nchi_points, nt])
    phi_w = np.zeros([ns_points, nchi_points, nw2])
    er_t = np.zeros([ns_points, nchi_points, nt])
    er_w = np.zeros([ns_points, nchi_points, nw2])
    phinz_t = np.zeros([ns_points, nchi_points, nt])
    phinz_w = np.zeros([ns_points, nchi_points, nw2])
    ernz_t = np.zeros([ns_points, nchi_points, nt])
    ernz_w = np.zeros([ns_points, nchi_points, nw2])

    def find_signal(dd, t, one_name):
        f_interp = interpolate.interp2d(
            dd['potsc_grids']['t'], dd['potsc_grids']['s'], dd[one_name]['data'].T)
        y = f_interp(t, s).T[:, id_s1]
        yw = ymath.fft_y(t, y, oo={'flag_f2_arranged': True})['f2_arranged']
        return y, yw

    for count_s1 in range(ns_points):
        s1 = s_points[count_s1]
        id_s1, s1 = mix.find(s, s1)
        s_points_res[count_s1] = s1

        # --- Zonal Er, Phi ---
        phibar_t[count_s1] = phibar[:, id_s1]
        erbar_t[count_s1]  = erbar[:, id_s1]

        phibar_w[count_s1] = ymath.fft_y(
            t, phibar_t[count_s1], oo={'flag_f2_arranged': True})['f2_arranged']
        erbar_w[count_s1] = ymath.fft_y(
            t, erbar_t[count_s1],  oo={'flag_f2_arranged': True})['f2_arranged']

        # --- Full and non-zonal Er and Phi ---
        for count_chi1 in range(nchi_points):
            # - Full Phi -
            phi_t[count_s1, count_chi1], phi_w[count_s1, count_chi1] = \
                find_signal(dd, t, names_potsc[count_chi1])
            
            # - Full Er -
            er_t[count_s1, count_chi1], er_w[count_s1, count_chi1] = \
                find_signal(dd, t, names_ersc[count_chi1])

            # - Nonzonal Phi -
            phinz_t[count_s1, count_chi1], phinz_w[count_s1, count_chi1] = \
                find_signal(dd, t, names_phinz[count_chi1])

            # - Nonzonal Er -
            ernz_t[count_s1, count_chi1], ernz_w[count_s1, count_chi1] = \
                find_signal(dd, t, names_ernz[count_chi1])

    # change normalization (the normalization should be performed just before the saving)
    t  = t  * coef_norm_t
    w2 = w2 * coef_norm_w

    # save results:
    wr.create_file(path_to_write)

    wr.save_data(path_to_write, 's_points', s_points_res,
                 desc=u'radial points, where the signals are taken')
    wr.save_data(path_to_write, 'chi_points', chi_points_res,
                 desc=u'poloidal points, where Full_Er is taken')
    wr.save_data(path_to_write, 't', t,
                 desc=u'time grid, seconds')
    wr.save_data(path_to_write, 'w', w2,
                 desc=u'frequency grid, kHz')

    wr.save_data(path_to_write, 'Zonal_Phi', phibar_t,
                 desc=u'electric potential (s_points, t), averaged on a flux surface')
    wr.save_data(path_to_write, 'FFT_Zonal_Phi', phibar_w,
                 desc=u'FFT of Zonal_Er (s_points, w)')
    wr.save_data(path_to_write, 'Zonal_Er', erbar_t,
                 desc=u'radial electric field (s_points, t), averaged on a flux surface')
    wr.save_data(path_to_write, 'FFT_Zonal_Er', erbar_w,
                 desc=u'FFT of Zonal_Er (s_points, w)')

    wr.save_data(path_to_write, 'Phi', phi_t,
                 desc=u'full electric potential (s_points, chi_point, t)')
    wr.save_data(path_to_write, 'FFT_Phi', phi_w,
                 desc=u'FFT of Phi (s_points, chi_point, w)')

    wr.save_data(path_to_write, 'Er', er_t,
                 desc=u'full radial electric field (s_points, chi_point, t)')
    wr.save_data(path_to_write, 'FFT_Er', er_w,
                 desc=u'FFT of Er (s_points, chi_point, w)')

    wr.save_data(path_to_write, 'Nonzonal_Phi', phinz_t,
                 desc=u'Phi - Zonal_Phi (s_points, chi_point, t)')
    wr.save_data(path_to_write, 'FFT_Nonzonal_Phi', phinz_w,
                 desc=u'FFT of Nonzonal_Phi (s_points, chi_point, w)')

    wr.save_data(path_to_write, 'Nonzonal_Er', ernz_t,
                 desc=u' - d(Phi - Zonal_Phi)/ds (s_points, chi_point, t)')
    wr.save_data(path_to_write, 'FFT_Nonzonal_Er', ernz_w,
                 desc=u'FFT of Nonzonal_Er (s_points, chi_point, w)')


def check_2d_st_interpolation(dd, oo):
    coef_norm_t = 1 / dd['wc']

    # --- Zonal Er ---
    zf.erbar(dd)
    erbar = dd['erbar']['data']
    s = dd['erbar']['s']
    t = dd['erbar']['t']

    # intervals for plotting
    t, ids_t = mix.get_array_oo(oo, t, 't')
    s, ids_s = mix.get_array_oo(oo, s, 's')
    erbar = mix.get_slice(erbar, ids_t, ids_s)

    # --- Full Phi and Er ---
    chi_points = oo.get('chi_points', [0.0])
    oo_phi = {'chi_s': chi_points}
    names_potsc = rd.potsc_chi(dd, oo_phi)

    phi_st = {}
    for one_name in names_potsc:
        f_interp = interpolate.interp2d(
            dd['potsc_grids']['t'], dd['potsc_grids']['s'], dd[one_name]['data'].T)
        phi_st[one_name] = f_interp(t, s).T

    # renormalize the time:
    t = t * coef_norm_t

    # Plot: erbar(s,t)
    curves = crv.Curves().xlab('t(seconds)').ylab('s').tit('Zonal\ Er')
    curves.new('f').XS(t).YS(s).ZS(erbar).lev(60)
    cpr.plot_curves_3d(curves)

    # Plot: original VS interpolated
    for one_name in names_potsc:
        s_grid, ids_s = mix.get_array_oo(oo, dd['potsc_grids']['s'], 's')
        t_grid, ids_t = mix.get_array_oo(oo, dd['potsc_grids']['t'], 't')
        phi_current = mix.get_slice(dd[one_name]['data'], ids_t, ids_s)
        curves = crv.Curves().xlab('t(seconds)').ylab('s').tit('original\ Phi')
        curves.new('f')\
            .XS(t_grid * coef_norm_t)\
            .YS(s_grid)\
            .ZS(phi_current)\
            .lev(60)
        cpr.plot_curves_3d(curves)

        curves = crv.Curves().xlab('t(seconds)').ylab('s').tit('interp.\ Phi')
        curves.new('f') \
            .XS(t) \
            .YS(s) \
            .ZS(phi_st[one_name]) \
            .lev(60)
        cpr.plot_curves_3d(curves)


def check_radial_structure_2d(dd, oo):
    coef_norm_t = 1 / dd['wc']

    # --- Full Phi(r,z) ---
    t_points = oo.get('t_points', [0.0])  # normalized to wc
    names_potsc_t = rd.potsc_t(dd, {'t_points': t_points})

    # --- Non-zonal Phi(r,z) ---
    rd.phibar(dd)
    phibar = dd['phibar']['data']
    s = dd['phibar']['s']
    t = dd['phibar']['t']
    f_interp = interpolate.interp2d(t, s, phibar.T)
    phibar_interp = f_interp(dd['potsc_grids']['t'], dd['potsc_grids']['s']).T

    count_t = -1
    phi_nz = np.zeros([np.size(t_points),
                       np.size(dd['potsc_grids']['chi']),
                       np.size(dd['potsc_grids']['s'])
                       ])
    for one_name in names_potsc_t:
        count_t += 1
        id_t_point = dd[one_name]['id_t_point']
        for id_chi in range(np.size(dd['potsc_grids']['chi'])):
            phi_nz[count_t, id_chi, :] = dd[one_name]['data'][id_chi, :] \
                                         - phibar_interp[id_t_point, :]

    # Plot Phi(r,z):
    count_t = -1
    for one_name in names_potsc_t:
        count_t += 1
        one_t_point = t_points[count_t]
        line_t = 't = {:0.3e}'.format(one_t_point * coef_norm_t)

        curves = crv.Curves().xlab('r').ylab('z').tit('Phi:\ ' + line_t)
        curves.new('f') \
            .XS(dd['potsc_grids']['r']) \
            .YS(dd['potsc_grids']['z']) \
            .ZS(dd[one_name]['data'].T) \
            .lev(60)
        cpr.plot_curves_3d(curves)

        curves = crv.Curves().xlab('r').ylab('z').tit('Non-zonal\ Phi:\ ' + line_t)
        curves.new('f') \
            .XS(dd['potsc_grids']['r']) \
            .YS(dd['potsc_grids']['z']) \
            .ZS(phi_nz[count_t, :, :].T) \
            .lev(60)
        cpr.plot_curves_3d(curves)


def check_interpolation_1d(dd, oo):
    s_point = oo.get('s_point', 0.0)
    chi_point = oo.get('chi_point', 0.0)

    coef_norm_t = 1 / dd['wc']

    f = h5.File(dd['path_to_write'], 'r')

    chi_points = np.array(f['chi_points'])
    s_points = np.array(f['s_points'])
    t = np.array(f['t'])  # already in seconds

    id_s_point,   s_point   = mix.find(s_points, s_point)
    id_chi_point, chi_point = mix.find(chi_points, chi_point)
    line_s1   = 's = {:0.3f}'.format(s_point)
    line_chi1 = 'chi = {:0.3f}'.format(chi_point)

    phi_interp = np.array(f['Phi'][id_s_point, id_chi_point, :])

    # --- Original Phi ---
    oo_phi = {'chi_s': [chi_point]}
    name_potsc = rd.potsc_chi(dd, oo_phi)[0]

    id_s_potsc, s_potsc = mix.find(dd['potsc_grids']['s'], s_point)
    line_s_potsc = 's = {:0.3f}'.format(s_potsc)

    # Plot: Phi and Er:
    curves_phi = crv.Curves().xlab('t(seconds)').tit('Phi:\ ' + line_chi1)
    curves_phi.new('orig')\
        .XS(dd['potsc_grids']['t'] * coef_norm_t)\
        .YS(dd[name_potsc]['data'][:, id_s_potsc])\
        .leg('original:\ ' + line_s_potsc)
    curves_phi.new('interp')\
        .XS(t)\
        .YS(phi_interp)\
        .leg('interp.:\ ' + line_s1).sty(':')

    cpr.plot_curves(curves_phi)

    # close the file:
    f.close()


def check_signals_1d(dd, oo):
    s_point = oo.get('s_point', 0.0)

    f = h5.File(dd['path_to_write'], 'r')

    chi_points = np.array(f['chi_points'])
    s_points   = np.array(f['s_points'])
    t = np.array(f['t'])

    nchi_points = np.size(chi_points)

    id_s_point, s_point = mix.find(s_points, s_point)
    line_s1 = 's = {:0.3f}'.format(s_point)
    phibar = np.array(f['Zonal_Phi'][id_s_point, :])
    erbar = np.array(f['Zonal_Er'][id_s_point, :])
    phi = np.array(f['Phi'][id_s_point, :, :])
    er = np.array(f['Er'][id_s_point, :, :])
    phinz = np.array(f['Nonzonal_Phi'][id_s_point, :, :])
    ernz = np.array(f['Nonzonal_Er'][id_s_point, :, :])

    # Plot: Zonal Phi and Er:
    curves = crv.Curves().xlab('w(kHz)').tit('Zonal\ Phi:\ ' + line_s1)
    curves.new('f').XS(t).YS(phibar)
    cpr.plot_curves(curves)

    curves = crv.Curves().xlab('w(kHz)').tit('Zonal\ Er:\ ' + line_s1)
    curves.new('f').XS(t).YS(erbar)
    cpr.plot_curves(curves)

    # Plot: Phi and Er:
    curves_phi = crv.Curves().xlab('w(kHz)') \
        .tit('Phi:\ ' + line_s1)
    curves_er = crv.Curves().xlab('w(kHz)') \
        .tit('Er:\ ' + line_s1)
    curves_phinz = crv.Curves().xlab('w(kHz)') \
        .tit('Non-zonal\ Phi:\ ' + line_s1)
    curves_ernz = crv.Curves().xlab('w(kHz)') \
        .tit('Non-zonal\ Er:\ ' + line_s1)
    for count_chi in range(nchi_points):
        one_chi  = chi_points[count_chi]
        line_chi = 'chi = {:0.3f}'.format(one_chi)

        # - Phi -
        curves_phi.new('f').XS(t).YS(phi[count_chi, :]).leg(line_chi)

        # - Er -
        curves_er.new('f').XS(t).YS(er[count_chi, :]).leg(line_chi)

        # - Non-zonal Phi -
        curves_phinz.new('f').XS(t).YS(phinz[count_chi, :]).leg(line_chi)

        # - Non-zonal Er -
        curves_ernz.new('f').XS(t).YS(ernz[count_chi, :]).leg(line_chi)
    cpr.plot_curves(curves_phi)
    cpr.plot_curves(curves_er)
    cpr.plot_curves(curves_phinz)
    cpr.plot_curves(curves_ernz)

    # close the file:
    f.close()


def check_fft_1d(dd, oo):
    s_point = oo.get('s_point', 0.0)

    f = h5.File(dd['path_to_write'], 'r')

    chi_points = np.array(f['chi_points'])
    s_points   = np.array(f['s_points'])
    w = np.array(f['w'])
    w, ids_w = mix.get_array_oo(oo, w, 'w')

    nchi_points = np.size(chi_points)

    id_s_point, s_point = mix.find(s_points, s_point)
    line_s1 = 's = {:0.3f}'.format(s_point)

    phibar = np.array(f['FFT_Zonal_Phi'][id_s_point, ids_w[0]:ids_w[-1]+1])
    erbar = np.array(f['FFT_Zonal_Er'][id_s_point, ids_w[0]:ids_w[-1]+1])
    phi = np.array(f['FFT_Phi'][id_s_point, :, ids_w[0]:ids_w[-1]+1])
    er = np.array(f['FFT_Er'][id_s_point, :, ids_w[0]:ids_w[-1]+1])
    phinz = np.array(f['FFT_Nonzonal_Phi'][id_s_point, :, ids_w[0]:ids_w[-1] + 1])
    ernz = np.array(f['FFT_Nonzonal_Er'][id_s_point, :, ids_w[0]:ids_w[-1] + 1])

    # Plot: Zonal Er and Phi:
    curves = crv.Curves()\
        .xlab('\omega(kHz)').tit('FFT:\ Zonal\ Phi:\ ' + line_s1)
    curves.new('f').XS(w).YS(phibar).leg(line_s1)
    cpr.plot_curves(curves)

    curves = crv.Curves().xlab('\omega(kHz)').tit('FFT:\ Zonal\ Er:\ ' + line_s1)
    curves.new('f').XS(w).YS(erbar).leg(line_s1)
    cpr.plot_curves(curves)

    # Plot: Phi and Er:
    curves_phi = crv.Curves().xlab('\omega(kHz)') \
        .tit('FFT:\ Phi:\ ' + line_s1)
    curves_er = crv.Curves().xlab('\omega(kHz)') \
        .tit('FFT:\ Er:\ ' + line_s1)
    curves_phinz = crv.Curves().xlab('\omega(kHz)') \
        .tit('FFT:\ Nonzonal\ Phi:\ ' + line_s1)
    curves_ernz = crv.Curves().xlab('\omega(kHz)') \
        .tit('FFT:\ Nonzonal\ Er:\ ' + line_s1)
    for count_chi in range(nchi_points):
        one_chi = chi_points[count_chi]
        line_chi = 'chi = {:0.3f}'.format(one_chi)

        # - Phi -
        curves_phi.new().XS(w).YS(phi[count_chi, :]).leg(line_chi)

        # - Er -
        curves_er.new().XS(w).YS(er[count_chi, :]).leg(line_chi)

        # - Nonzonal Phi -
        curves_phinz.new().XS(w).YS(phinz[count_chi, :]).leg(line_chi)

        # - Nonzonal Er -
        curves_ernz.new().XS(w).YS(ernz[count_chi, :]).leg(line_chi)
    cpr.plot_curves(curves_phi)
    cpr.plot_curves(curves_er)
    cpr.plot_curves(curves_phinz)
    cpr.plot_curves(curves_ernz)

    # close the file:
    f.close()

























