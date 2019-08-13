import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import zf_gam as zf
import general as gn
import write_data as wr
import transport
import write_data
import numpy as np
import h5py as h5
from scipy.interpolate import BSpline
from scipy.interpolate import PPoly
import sys


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
    mix.reload_module(write_data)


def init(dd):
    if 'fields' in dd['3d']:
        return

    # --- read necessary additional parameters ---
    read_add(dd)

    # --- create 3d-fields field ---
    d3 = dd['3d']
    d3['fields'] = {}
    ff = dd['3d']['fields']

    # --- read 3d-field data ---
    path_to_file = dd['path'] + '/orb5_res_parallel.h5'
    f = h5.File(path_to_file, 'r')

    ff['t']             = np.array(f['/data/var3d/generic/pot3d/time'])
    ff['m_mins']        = np.array(f['/data/var3d/generic/pot3d/mmin'])  # [n, s]
    coef_bspl_dft = np.array(f['/data/var3d/generic/pot3d/data'])  # [t, n, s, m]

    ff['coef_bspl_dft'] = coef_bspl_dft[:]['real'] + 1j * coef_bspl_dft[:]['imaginary']

    f.close()

    # --- create poloidal and toroidal Fourier transform of the field ---
    form_grids(d3, ff)


def init_antenna(dd, name_antenna_file='orb5_res_parallel.h5',
                 structure_name='antenna', flag_orb_shape=False):
    # if structure_name in dd['3d']:
    #     return

    # --- read necessary additional parameters ---
    read_add(dd)

    # --- create antenna field ---
    d3 = dd['3d']
    d3[structure_name] = {}
    ffa = dd['3d'][structure_name]

    if structure_name == 'antenna-init':
        ffa['nsel_type'] = 'phi'

    # --- define antenna parameters ---
    frequency = None

    # # --- read input antenna parameters ---
    # path_to_file = dd['path'] + '/orb5_res.h5'
    # f = h5.File(path_to_file, 'r')
    # try:
    #     ffa['amplitudes']   = np.array(f['/parameters/antenna/amplitudes'])
    #     ffa['n']            = np.array(f['/parameters/antenna/n'])
    #     ffa['m']            = np.array(f['/parameters/antenna/m'])
    #
    #     # antenna's type
    #     antenna_type = f['/parameters/antenna/nsel_type'].attrs
    #     ids_attr = list(antenna_type)
    #     ffa['nsel_type'] = [antenna_type[name].decode("utf-8")
    #                            for name in ids_attr[0:len(ids_attr) - 1]][0]
    # except:
    #     f.close()
    # f.close()

    # --- read initial antenna structure ---
    path_to_file = dd['path'] + '/' + name_antenna_file
    f = h5.File(path_to_file, 'r')
    try:
        path_parameter = '/parameters/w'
        if path_parameter in f:
            frequency = np.array(f[path_parameter])[0]
            ffa['w']  = frequency

        path_antenna = '/data/var3d/generic/' + ffa['nsel_type'] + '_antenna/run.1/'
        ffa['m_mins']   = np.array(f[path_antenna + 'mmin'])      # [n, s]

        nn = np.shape(ffa['m_mins'])[0]
        ns = np.shape(ffa['m_mins'])[1]
        nm = 2 * d3['deltam'] + 1

        if not flag_orb_shape:
            if frequency == 0:
                path_str = path_antenna + 'bspl_dft'
                ffa['coef_bspl_dft'] = np.array(f[path_str])[None, :, :, :]  # [dummy_t, n, s, m]
                ffa['t'] = np.zeros(1)
            else:
                ffa['coef_bspl_dft'] = np.zeros([2, nn, ns, nm], dtype=np.complex64)

                path_str = path_antenna + 'bspl_dft_t0'
                ffa['coef_bspl_dft'][0, :, :, :] = np.array(f[path_str])[:, :, :]  # [dummy_t, n, s, m]

                path_str = path_antenna + 'bspl_dft_tquarter'
                ffa['coef_bspl_dft'][1, :, :, :] = np.array(f[path_str])[:, :, :]  # [dummy_t, n, s, m]

                T = 2*np.pi/frequency
                ffa['t'] = np.array([0.0, 0.25 * T])
        else:
            ffa['coef_bspl_dft'] = np.zeros([2, nn, ns, nm], dtype=np.complex64)

            path_str = path_antenna + 'orb_bspl_dft_t0'
            orb_shape_bspl_dft_t0 = np.array(f[path_str])[:, :]  # [ns * nm, nn]

            path_str = path_antenna + 'orb_bspl_dft_tquarter'
            orb_shape_bspl_dft_tquarter = np.array(f[path_str])[:, :]  # [ns * nm, nn]

            orb_shape_bspl_dft_t0 = \
                orb_shape_bspl_dft_t0[:]['real'] + 1j * orb_shape_bspl_dft_t0[:]['imaginary']
            orb_shape_bspl_dft_tquarter = \
                orb_shape_bspl_dft_tquarter[:]['real'] + 1j * orb_shape_bspl_dft_tquarter[:]['imaginary']

            ffa['coef_bspl_dft'][0, :, :, :] = orb_shape_bspl_dft_t0.reshape((nn, ns, nm))
            ffa['coef_bspl_dft'][1, :, :, :] = orb_shape_bspl_dft_tquarter.reshape((nn, ns, nm))

            T = 2 * np.pi / frequency
            ffa['t'] = np.array([0.0, 0.25 * T])
    except:
        f.close()
        sys.exit(-1)
    f.close()

    # --- create poloidal and toroidal Fourier transform of the field ---
    form_grids(d3, ffa)


def read_add(dd):
    if 'nidbas' in dd['3d']:
        return

    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')

    # order of splines
    nidbas = np.array(f['/parameters/fields/nidbas'])[0].astype(np.int)

    # number of toroidal grid points
    nphi = np.array(f['/parameters/fields/nphi'])[0].astype(np.int)

    # number of poloidal grid points
    nchi = np.array(f['/parameters/fields/nchi'])[0].astype(np.int)

    # left and right boundaries of toroidal filter
    nfilt1 = np.array(f['/parameters/fields/nfilt1'])[0].astype(np.int)
    nfilt2 = np.array(f['/parameters/fields/nfilt2'])[0].astype(np.int)

    # width of the poloidal band:
    deltam = np.array(f['/parameters/fields/deltam'])[0].astype(np.int)

    # grid of toroidal grids
    phi_grid = np.array(f['/data/var2d/generic/potsp/coord2'])

    # step of the flux-tube filter
    n_flux_tube = np.array(f['/parameters/fields/n_flux_tube'])[0].astype(np.int)

    # radial boundaries
    sfmin = np.array(f['/parameters/fields/sfmin'])[0]
    sfmax = np.array(f['/parameters/fields/sfmax'])[0]

    # close file
    f.close()

    # grids:
    rd.potsc_grids(dd)
    s   = dd['potsc_grids']['s']
    chi = dd['potsc_grids']['chi']
    r = dd['potsc_grids']['r']
    z = dd['potsc_grids']['z']

    # results
    dd['3d'].update({
        'nidbas': nidbas,
        'nphi': nphi, 'nchi': nchi,
        'nfilt1': nfilt1, 'nfilt2': nfilt2,
        'deltam': deltam,
        'n_flux_tube': n_flux_tube,
        'sfmin': sfmin, 'sfmax': sfmax,
        's': s, 'chi': chi, 'phi': phi_grid,
        'r': r, 'z': z,
    })


def form_grids(d3, ff, sel_opt_bspl_s='opt2'):
    coef_bspl_dft = ff['coef_bspl_dft']

    n_nmodes = np.shape(coef_bspl_dft)[1]  # number of toroidal modes
    ns       = np.shape(coef_bspl_dft)[2]  # number of radial points
    m_width  = np.shape(coef_bspl_dft)[3]  # number of poloidal modes for every n-mode

    s = d3['s']
    ds = (d3['sfmax'] - d3['sfmin']) / (np.size(s) - 1)
    nidbas = d3['nidbas']

    # result signal: poloidally and toroidally Fourier transformed field:
    data_ft = None

    # --- find absolute values of poloidal modes ---
    m_offsets = np.meshgrid(np.arange(m_width), np.arange(ns))[0]
    m_modes   = np.zeros([n_nmodes, ns, m_width])
    for id_n_mode in np.arange(n_nmodes):
        loc_m_mins = np.array([ff['m_mins'][id_n_mode]]).transpose()
        m_modes[id_n_mode, :, :] = loc_m_mins + m_offsets

    # --- get toroidal mode numbers ---
    n_modes      = [i for i in range(d3['nfilt1'], d3['nfilt2']+1)]
    n_modes_flux = [
        i for i in range(d3['nfilt1'], d3['nfilt2']+1) if i % d3['n_flux_tube'] == 0
    ]

    # --- get poloidal and toroidal Fourier transformation of signal, but still Bspline in radial direction ---
    # s_bspl_ft = coef_bspl_dft[:]['real'] + 1j * coef_bspl_dft[:]['imaginary']
    s_bspl_ft = np.array(coef_bspl_dft)
    for i_n_mode, ntor in enumerate(n_modes):
        if ntor not in n_modes_flux:
            continue
        m_modes_n1 = m_modes[i_n_mode, :, :]  # [s, m_width]

        s_bspl_ft[:, i_n_mode, :, :] *= \
            ( np.sinc(ntor/d3['nphi']) * np.sinc(m_modes_n1/d3['nchi']) )**(d3['nidbas']+1) * \
                np.exp(1j * np.pi * (d3['nidbas'] - 1) * (m_modes_n1 / d3['nchi'] + ntor / d3['nphi']))

    # Include negative n-modes:
    for i_n_mode, ntor in enumerate(n_modes):
        if ntor not in n_modes_flux:
            continue
        if ntor > 0:
            s_bspl_ft[:, i_n_mode, :, :] = 2 * s_bspl_ft[:, i_n_mode, :, :]

    # data_ft = s_bspl_ft

    # --- sum radial Bspline in radial direction ---
    # OPTION 1:
    if sel_opt_bspl_s == 'opt1':
        basic_bspl_obj = BSpline.basis_element(np.array([i for i in range(-nidbas, 1 + 1)]))
        basic_bspl_pp = PPoly.from_spline(basic_bspl_obj)

        wx = (s - ff['sfmin']) / ds  # we use s-grid from potsc, so here we always get integer numbers
        id_grid_point = np.floor(wx)

        data_ft = np.zeros_like(s_bspl_ft)
        for id_knot in range(nidbas):
            value_bspl = basic_bspl_pp(wx - id_grid_point - id_knot)
            data_ft[:, :, :, :] += value_bspl[None, None, :, None] * s_bspl_ft[:, :, :, :]

    # OPTION 2:
    if sel_opt_bspl_s == 'opt2':
        add_knots_1 = np.array([- i * ds / 2. for i in range(nidbas, 0, -1)])
        s_knots = np.concatenate((add_knots_1, s + ds / 2))
        obj_bspl = BSpline(s_knots, s_bspl_ft, nidbas, axis=2)
        data_ft = obj_bspl(s)

    # results:
    ff.update({
        'n_modes': n_modes,
        'n_modes_flux': n_modes_flux,
        'm_modes': m_modes,  # [n, s, m_width]
        's': s, 'chi': d3['chi'], 'phi': d3['phi'],
        'data_ft': data_ft,    # [t, n, s, m]
    })


def save_rhsf(dd, freq_n1_wc, A_scaling=1, name_file_antenna='user_antenna.h5'):
    # freq_n1 - frequency of a toroidal mode n1
    # A_scaling - scaling of the result structures
    # name_file_antenna - name of a file where one is going to save
    #                       the antenna's radial structure

    # --- mode period ---
    T_period = 2 * np.pi / freq_n1_wc

    # --- read 3d-field data ---
    init(dd)
    d3 = dd['3d']
    ff = d3['fields']

    mmins = ff['m_mins']

    # time point of the quarter of the n1 mode's period:
    id_t_quarter, t_quarter, _ = mix.get_ids(ff['t'], ff['t'][-1] - 0.75 * T_period)

    # radial structure at last time point
    coef_bspl_dft_BEGIN   = ff['coef_bspl_dft'][-1, :, :, :] * A_scaling

    # radial structure in a quarter of the n1 mode's period
    coef_bspl_dft_QUARTER = ff['coef_bspl_dft'][id_t_quarter, :, :, :] * A_scaling

    # --- SAVE ANTENNA'S STRUCTURE ---
    path_to_file = dd['path'] + name_file_antenna
    ff = h5.File(path_to_file, 'w')

    try:
        # save frequency and scaling
        grp = ff.create_group("parameters")

        ddata = grp.create_dataset('w', data=np.array([freq_n1_wc]))
        ddata.attrs[u'descr'] = 'Frequency of a chosen toroidal mode.'

        ddata = grp.create_dataset('A', data=np.array([A_scaling]))
        ddata.attrs[u'descr'] = 'Scale of a chosen toroidal mode.'

        ddata = grp.create_dataset('number_n_modes', data=1)
        ddata.attrs[u'descr'] = 'Number of n modes.'

        # save radial structure
        grp = ff.create_group("data/var3d/generic/phi_antenna/run.1")

        ddata = grp.create_dataset('mmin', data=mmins)
        ddata.attrs[u'descr'] = 'Minimal poloidal mode numbers.'

        if freq_n1_wc > 0:
            ddata = grp.create_dataset('bspl_dft_t0', data=coef_bspl_dft_BEGIN)
            ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients: t = 0.'

            ddata = grp.create_dataset('bspl_dft_tquarter', data=coef_bspl_dft_QUARTER)
            ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients: t = 0.25 T.'
        else:
            ddata = grp.create_dataset('bspl_dft', data=coef_bspl_dft_BEGIN)
            ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients.'

        # change the shape:
        nn_bspl = np.shape(coef_bspl_dft_BEGIN)[0]
        ns_bspl = np.shape(coef_bspl_dft_BEGIN)[1]
        wm_bspl = np.shape(coef_bspl_dft_BEGIN)[2]
        neq_f = ns_bspl*wm_bspl

        # adjust shape of the matrix to the shape used in ORB5
        new_bspl_dft_BEGIN      = coef_bspl_dft_BEGIN.reshape((nn_bspl, neq_f))
        new_bspl_dft_QUARTER    = coef_bspl_dft_QUARTER.reshape((nn_bspl, neq_f))

        # type compatible with the complex-like type that ORB5 understands
        comp_datatype = np.dtype([
            ('real', np.float),
            ('imaginary', np.float)
        ])

        temp_t0 = np.zeros(np.shape(new_bspl_dft_BEGIN), dtype=comp_datatype)
        temp_tq = np.zeros(np.shape(new_bspl_dft_QUARTER), dtype=comp_datatype)
        for id_n in range(nn_bspl):
            for id_eq in range(neq_f):
                temp_t0[id_n, id_eq] = (
                    new_bspl_dft_BEGIN[id_n, id_eq].real,
                    new_bspl_dft_BEGIN[id_n, id_eq].imag
                )
                temp_tq[id_n, id_eq] = (
                    new_bspl_dft_QUARTER[id_n, id_eq].real,
                    new_bspl_dft_QUARTER[id_n, id_eq].imag
                )

        ddata = grp.create_dataset('orb_bspl_dft_t0', data=temp_t0)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients: t = 0. ' \
                                'Shape is [(ns+nidbas)*width_m, number_of_n_modes]'

        ddata = grp.create_dataset('orb_bspl_dft_tquarter', data=temp_tq)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients: t = 0.25 T.' \
                                'Shape is [(ns+nidbas)*width_m, number_of_n_modes]'

        # test_array = [(10, 0.2), (20, 0.4)]
        # test_array = np.array(test_array, dtype=comp_datatype)
        # ddata = grp.create_dataset('test_array', data=test_array, dtype=comp_datatype)
        # ddata.attrs[u'descr'] = 'Test array'
        #
        # test_array2 = np.zeros(2)
        # test_array2[0] = 10.2
        # test_array2[1] = 20.4
        # ddata = grp.create_dataset('test_array2', data=test_array2)
        # ddata.attrs[u'descr'] = 'Test array2'

    except:
        ff.close()
        sys.exit(-1)

    ff.close()

    return


def signal_chi1_phi1(ff, chi_point, phi_point=0.0):
    data_ft = ff['data_ft']  # [t, n, s, m]

    t = ff['t']
    s = ff['s']

    n_modes      = ff['n_modes']
    n_modes_flux = ff['n_modes_flux']
    m_modes      = ff['m_modes']  # [n, s, m_width]

    # Get a signal in real poloidal and toroidal coordinates:
    signal_ts = np.zeros([np.size(t), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        for id_n, ntor in enumerate(n_modes):
            if ntor not in n_modes_flux:
                continue
            for id_m, m_mode in enumerate(m_modes[id_n, id_s1, :]):
                signal_ts[:, id_s1] += data_ft[:, id_n, id_s1, id_m] \
                            * np.exp(1j * m_mode * chi_point + 1j * ntor * phi_point)

    return np.real(signal_ts)


def signal_t1_phi1(ff, t_point, phi_point=0.0):
    data_ft = ff['data_ft']  # [t, n, s, m]

    t = ff['t']
    s = ff['s']
    chi = ff['chi']

    n_modes      = ff['n_modes']
    n_modes_flux = ff['n_modes_flux']
    m_modes      = ff['m_modes']  # [n, s, m_width]

    id_t1, t1, _ = mix.get_ids(t, t_point)

    # Get a signal in real poloidal and toroidal coordinates:
    signal_chis = np.zeros([np.size(chi), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        for id_n, ntor in enumerate(n_modes):
            if ntor not in n_modes_flux:
                continue
            for id_m, m_mode in enumerate(m_modes[id_n, id_s1, :]):
                signal_chis[:, id_s1] += data_ft[id_t1, id_n, id_s1, id_m] \
                        * np.exp(1j * m_mode * chi) * np.exp(1j * ntor * phi_point)

    return np.real(signal_chis), t1


def choose_one_var_ts(ovar, dd):
    opt_var = ovar[0]
    vvar, tit_var = None, ''
    s, t = None, None
    res = {}

    if opt_var == 'potsc':
        chi_point = ovar[1]

        phi_point = 0.0
        if len(ovar) > 2:
            phi_point = ovar[2]

        init(dd)
        ff = dd['3d']['fields']
        s = ff['s']
        t = ff['t']

        vvar = signal_chi1_phi1(ff, chi_point, phi_point)
        line_chi = '\chi = {:0.1f},\ '.format(chi_point)
        line_phi = '\\varphi = {:0.1f}'.format(phi_point)
        tit_var = '3D:\ \Phi' + '_{' + line_chi + line_phi + '}'
    if opt_var == 'potsc-antenna-init':
        chi_point = ovar[1]
        file_name = ovar[2]
        phi_point = 0.0
        flag_orb_shape = False
        if len(ovar) > 3:
            phi_point = ovar[3]
        if len(ovar) > 4:
            flag_orb_shape = ovar[4]

        structure_name = 'antenna-init'

        init_antenna(dd, name_antenna_file=file_name, structure_name=structure_name,
                     flag_orb_shape=flag_orb_shape)
        ff = dd['3d'][structure_name]
        s = ff['s']
        t = ff['t']

        vvar = signal_chi1_phi1(ff, chi_point, phi_point)
        line_chi = '\chi = {:0.1f},\ '.format(chi_point)
        line_phi = '\\varphi = {:0.1f}'.format(phi_point)
        tit_var = '3D:\ Init.\ Antenna:\ \Phi' + '_{' + line_chi + line_phi + '}'

    # result:
    res.update({
        'data': vvar,
        's': s,
        't': t,
        'tit': tit_var
    })

    return res


def choose_one_var_rz(ovar, dd):
    opt_var   = ovar[0]
    vvar, s, chi, tit_var = None, None, None, ''
    res = {}

    if opt_var == 'potsc':
        t_point = ovar[1]

        phi_point = 0.0
        if len(ovar) > 2:
            phi_point = ovar[2]

        init(dd)
        ff = dd['3d']['fields']
        s   = ff['s']
        chi = ff['chi']

        vvar, t1 = signal_t1_phi1(ff, t_point, phi_point)
        line_time = 't = {:0.3e},\ '.format(t1)
        line_phi = '\\varphi = {:0.1f}'.format(phi_point)
        tit_var = '3D:\ \Phi' + '_{' + line_time + line_phi + '}'
    if opt_var == 'potsc-antenna-init':
        t_point = ovar[1]
        file_name = ovar[2]
        phi_point = 0.0
        flag_orb_shape = False
        if len(ovar) > 3:
            phi_point = ovar[3]
        if len(ovar) > 4:
            flag_orb_shape = ovar[4]

        structure_name = 'antenna-init'

        init_antenna(dd, name_antenna_file=file_name, structure_name=structure_name,
                     flag_orb_shape=flag_orb_shape)
        ff = dd['3d'][structure_name]
        s   = ff['s']
        chi = ff['chi']

        vvar, t1 = signal_t1_phi1(ff, t_point, phi_point)
        line_time = 't = {:0.3e},\ '.format(t1)
        line_phi = '\\varphi = {:0.1f}'.format(phi_point)
        tit_var = '3D:\ Init.\ Antenna:\ \Phi' + '_{' + line_time + line_phi + '}'

    # results
    R    = dd['3d']['r']
    Z    = dd['3d']['z']

    res.update({
        'data': vvar.T,  # (s, chi)
        'r': R,
        'z': Z,
        's': s,
        'chi': chi,
        'tit': tit_var
    })

    return res

