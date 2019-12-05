import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import zf_gam as zf
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
    mix.reload_module(wr)
    mix.reload_module(transport)
    mix.reload_module(write_data)


# # Get toroidal modes, which are actually simulated in a plasma system
# def get_n_modes_old(n_modes, n_modes_flux):
#     ids_n_global, n_global = [], []
#     for i_n_mode, ntor in enumerate(n_modes):
#         if ntor not in n_modes_flux:
#             continue
#         ids_n_global.append(i_n_mode)
#         n_global.append(ntor)
#     return np.array(n_global), np.array(ids_n_global)
#
#
# # Get antenna toroidal modes, which are actually simulated in a plasma system
# def get_antenna_n_modes(n_modes, n_modes_flux, n_saved_ant):
#     # antenna toroidal modes, which have been actually simulated in a plasma system
#     n_ant = []
#
#     # indices of the antenna toroidal modes wrt n_saved_ant
#     ids_n_ant = []
#
#     # indices of the antenna toroidal modes wrt n_modes
#     ids_n_ant_global = []
#
#     count_id_n_ant = -1
#     n_global, ids_n_global = get_n_modes(n_modes, n_modes_flux)
#     for i_n_mode, ntor in enumerate(n_global):
#         if ntor in n_saved_ant:
#             count_id_n_ant += 1
#             n_ant.append(ntor)
#             ids_n_ant.append(count_id_n_ant)
#             ids_n_ant_global.append(ids_n_global[i_n_mode])
#     return np.array(n_ant), np.array(ids_n_ant), np.array(ids_n_ant_global)


# Get toroidal modes:
def get_n_modes(n_modes, n_modes_flux, n_mask=None):
    if n_mask is None:
        n_mask = n_modes

    n, ids_n, ids_n_mask = [], [], []
    count_n_mask = -1
    for i_n_mode, ntor in enumerate(n_modes):
        if ntor not in n_modes_flux:
            continue
        if ntor in n_mask:
            count_n_mask += 1
            n.append(ntor)
            ids_n.append(i_n_mode)
            ids_n_mask.append(count_n_mask)
    return np.array(n), np.array(ids_n), np.array(ids_n_mask)


# Sum radial Bspline in radial direction
def sum_rad_bspline(sel_opt_bspl_s, d3, ff, s_bspl_ft):
    s = d3['s']
    ds = (d3['sfmax'] - d3['sfmin']) / (np.size(s) - 1)
    nidbas = d3['nidbas']

    data_ft = None

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

    return data_ft  # [t, n, s, m]


# Save antenna structure
def save_rhsf(dd, n_save, freqs_ant_wc, A_scaling=None, name_file_antenna='user_antenna.h5'):
    # n_save - array of toroidal modes to save.
    # freqs_ant_wc - array of frequencies of a toroidal modes n_save.
    # A_scaling - array of scalings of the result structures.
    # name_file_antenna - name of a file where one is going to save
    #                       the antenna's radial structure

    def save_coef_bspl_dft(ids_t):
        coef_bspl_dft = np.zeros([len(n_save_res), ns_bspl, nm], dtype=np.complex128)
        for id_n_mode in range(len(n_save_res)):
            coef_bspl_dft[id_n_mode, :, :] = \
                ff['coef_bspl_dft'][ids_t[id_n_mode], id_n_mode, :, :] * A_scaling_res[id_n_mode]
        return coef_bspl_dft

    # --- read 3d-field data ---
    init(dd)
    d3 = dd['3d']
    ff = d3['fields']

    mmins        = ff['m_mins']
    ns_bspl = np.shape(mmins)[1]
    nm = 2 * d3['deltam'] + 1

    # antenna toroidal modes to save
    n_save_res, ids_n_ant_global, ids_n_ant = get_n_modes(d3['n_sim'], d3['n_sim_flux'], n_save)

    # minimum poloidal modes
    mmins_res = np.zeros( [len(n_save_res), np.shape(mmins)[1]] )
    for id_n1_ant, n1_ant in enumerate(n_save_res):
        mmins_res[id_n1_ant, :] = mmins_res[ids_n_ant_global[id_n1_ant], :]

    # --- result antenna frequencies and periods ---
    freqs_ant_res = np.zeros(len(n_save_res))
    for counter_loc, id_ant in enumerate(ids_n_ant):
        freqs_ant_res[counter_loc] = freqs_ant_wc[id_ant]
    T_periods_res = 2 * np.pi / freqs_ant_res

    # time points of the quarter of antenna toroidal mode's periods:
    ids_t_quarters, ts_quarters = np.zeros(len(n_save_res), dtype=np.int), np.zeros(len(n_save_res))
    for counter_loc, T_period in enumerate(T_periods_res):
        ids_t_quarters[counter_loc], ts_quarters[counter_loc], _ = \
            mix.get_ids(ff['t'], ff['t'][-1] - 0.75 * T_period)

    # antenna scaling
    A_scaling_res = np.zeros(len(n_save_res))
    if A_scaling is None:
        A_scaling_res = np.ones(len(n_save_res))
    else:
        for counter_loc, id_ant in enumerate(ids_n_ant):
            A_scaling_res[counter_loc] = A_scaling[id_ant]

    # radial structure at last time point
    coef_bspl_dft_BEGIN = save_coef_bspl_dft([-1] * len(n_save_res))

    # radial structure in a quarter of the n1 mode's period
    coef_bspl_dft_QUARTER = save_coef_bspl_dft(ids_t_quarters)

    # --- SAVE ANTENNA'S STRUCTURE ---
    path_to_file = dd['path'] + '/' + name_file_antenna
    ff = h5.File(path_to_file, 'w')
    try:
        # save toroidal modes, their frequencies and scalings
        grp = ff.create_group("parameters")

        ddata = grp.create_dataset('n_saved', data=n_save_res)
        ddata.attrs[u'descr'] = 'Saved toroidal modes.'

        ddata = grp.create_dataset('w', data=np.array([freqs_ant_res]))
        ddata.attrs[u'descr'] = 'Frequencies of chosen toroidal modes.'

        ddata = grp.create_dataset('A', data=np.array([A_scaling]))
        ddata.attrs[u'descr'] = 'Scales of chosen toroidal modes.'

        # save radial structure
        grp = ff.create_group("data/var3d/generic/phi_antenna/run.1")

        ddata = grp.create_dataset('mmin', data=mmins_res)
        ddata.attrs[u'descr'] = 'Minimal poloidal mode numbers.'

        ddata = grp.create_dataset('bspl_dft_t0', data=coef_bspl_dft_BEGIN)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients: t = 0.'

        ddata = grp.create_dataset('bspl_dft_tquarter', data=coef_bspl_dft_QUARTER)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients: t = 0.25 T.'

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

        # array of shape [nn, neq] (when ORB5 reads the .h5 file, it will automatically transpose the matrix)
        temp_t0 = np.zeros(np.shape(new_bspl_dft_BEGIN),   dtype=comp_datatype)
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
    except:
        ff.close()
        sys.exit(-1)

    ff.close()

    return


def init(dd):
    if 'fields' in dd['3d']:
        return

    # --- read necessary additional parameters ---
    read_add(dd)

    # --- create 3d-fields field ---
    d3 = dd['3d']
    d3['fields'] = {}
    ff = dd['3d']['fields']

    ff['structure_name'] = 'fields'

    # --- read 3d-field data ---
    path_to_file = dd['path'] + '/orb5_res_parallel.h5'
    f = h5.File(path_to_file, 'r')

    ff['t']       = np.array(f['/data/var3d/generic/pot3d/time'])
    ff['m_mins']  = np.array(f['/data/var3d/generic/pot3d/mmin'])  # [n, s]
    coef_bspl_dft = np.array(f['/data/var3d/generic/pot3d/data'])  # [t, n, s, m]

    ff['coef_bspl_dft'] = coef_bspl_dft[:]['real'] + 1j * coef_bspl_dft[:]['imaginary']

    f.close()

    # --- create poloidal and toroidal Fourier transform of the field ---
    form_fields_nm(d3, ff)


def init_antenna(dd, name_antenna_file='', structure_name='antenna', flag_orb_shape=False):
    def read_coef_bspl(id_t):
        name_bspl = 'bspl_dft_t0'
        if id_t == 1:
            name_bspl = 'bspl_dft_tquarter'
        temp_bspl = np.array(f[path_antenna + name_bspl])[:, :, :]  # [dummy_t, n, s, m]

        if temp_bspl.dtype == np.complex128:
            ffa['coef_bspl_dft'][id_t, :, :, :] = temp_bspl
        else:
            ffa['coef_bspl_dft'][id_t, :, :, :] = \
                temp_bspl[:]['real'] + 1j * temp_bspl[:]['imaginary']

    def read_coef_bspl_orb_shape(id_t):
        name_bspl = 'orb_bspl_dft_t0'
        if id_t == 1:
            name_bspl = 'orb_bspl_dft_tquarter'
        orb_shape_bspl = np.array(f[path_antenna + name_bspl])[:, :]  # [ns * nm, nn]
        orb_shape_bspl = orb_shape_bspl[:]['real'] + 1j * orb_shape_bspl[:]['imaginary']
        ffa['coef_bspl_dft'][id_t, :, :, :] = orb_shape_bspl.reshape((nn, ns, nm))

    # if structure_name in dd['3d']:
    #     return

    # --- read necessary additional parameters ---
    read_add(dd)

    # --- create antenna field ---
    d3 = dd['3d']
    d3[structure_name] = {}
    ffa = dd['3d'][structure_name]

    ffa['structure_name'] = structure_name

    # --- read input antenna parameters ---
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')
    try:
        # antenna's type
        ffa['nsel_type'] = 'none'
        path_parameter = '/parameters/antenna/nsel_type'
        if path_parameter in f:
            antenna_type = f[path_parameter].attrs
            ids_attr = list(antenna_type)
            ffa['nsel_type'] = [antenna_type[name].decode("utf-8")
                                   for name in ids_attr[0:len(ids_attr) - 1]][0]

        # type of antenna's radial profile
        ffa['nsel_profile'] = ''
        path_parameter = '/parameters/antenna/nsel_profile'
        if path_parameter in f:
            antenna_type_profile = f[path_parameter].attrs
            ids_attr = list(antenna_type_profile)
            ffa['nsel_profile'] = [antenna_type_profile[name].decode("utf-8")
                                for name in ids_attr[0:1]][0]

        # name of .h5 file with antenna initial structure
        ffa['file_name_rhsf'] = ''
        path_parameter = '/parameters/antenna/file_name_rhsf'
        if path_parameter in f:
            antenna_type = f[path_parameter].attrs
            ids_attr = list(antenna_type)
            ffa['file_name_rhsf'] = [antenna_type[name].decode("utf-8")
                                for name in ids_attr[0:1]][0]
        if ffa['nsel_type'] == 'none':
            ffa['file_name_rhsf'] = ''

        if structure_name == 'antenna-init':
            ffa['nsel_type'] = 'phi'
            ffa['nsel_profile'] = 'rhsf'

    except:
        f.close()
    f.close()

    # --- file name of antenna structure ---
    if name_antenna_file == '':
        if structure_name == 'antenna':
            name_antenna_file = 'orb5_res_parallel.h5'
        if structure_name == 'antenna-init':
            if ffa['file_name_rhsf'] == '':
                name_antenna_file = 'user_antenna.h5'
                ffa['file_name_rhsf'] = name_antenna_file
            else:
                name_antenna_file = ffa['file_name_rhsf']
    if ffa['file_name_rhsf'] == '' and structure_name == 'antenna-init':
        ffa['file_name_rhsf'] = name_antenna_file

    # --- read initial antenna structure ---
    path_to_file = dd['path'] + '/' + name_antenna_file
    f = h5.File(path_to_file, 'r')
    try:
        # read antenna parameters (from user_antenna.h5)
        if ffa['nsel_profile'] == 'rhsf':
            path_to_file = dd['path'] + '/' + ffa['file_name_rhsf']
            finit = h5.File(path_to_file, 'r')
            try:
                path_parameter = '/parameters/w'
                if path_parameter in finit:
                    ffa['w'] = np.array(finit[path_parameter])

                path_parameter = '/parameters/n_saved'
                if path_parameter in finit:
                    ffa['n_saved'] = np.array(finit[path_parameter])
            except:
                finit.close()
                sys.exit(-1)
            finit.close()

        # read antenna structure (from user_antenna.h5 or orb5_res_parallel.h5)
        path_antenna  = '/data/var3d/generic/' + ffa['nsel_type'] + '_antenna/run.1/'
        if path_antenna in f:
            ffa['m_mins'] = np.array(f[path_antenna + 'mmin'])  # [n, s]

            nn = np.shape(ffa['m_mins'])[0]
            ns = np.shape(ffa['m_mins'])[1]
            nm = 2 * d3['deltam'] + 1

            ffa['coef_bspl_dft'] = np.zeros([2, nn, ns, nm], dtype=np.complex128)
            if not flag_orb_shape:
                read_coef_bspl(0)
                read_coef_bspl(1)
            else:
                read_coef_bspl_orb_shape(0)
                read_coef_bspl_orb_shape(1)
        else:
            mix.error_mes('<' + ffa['nsel_type'] + '> is a wrong antenna type')
    except:
        f.close()
        sys.exit(-1)
    f.close()

    # T        = 2 * np.pi / frequency
    # ffa['t'] = np.array([0.0, 0.25 * T])

    # since antenna time evolution is defined only by its initial frequencies,
    # we can use any time grid for the antenna
    rd.phibar(dd)
    ffa['t'] = dd['phibar']['t']

    # --- create poloidal and toroidal Fourier transform of the field ---
    form_fields_nm(d3, ffa)


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

    # toroidal modes
    n_modes = [i for i in range(nfilt1, nfilt2 + 1)]
    n_modes_flux = [
        i for i in range(nfilt1, nfilt2 + 1) if i % n_flux_tube == 0
    ]

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
        'n_sim': n_modes,
        'n_sim_flux': n_modes_flux,
    })


def form_fields_nm(d3, ff, sel_opt_bspl_s='opt2'):
    coef_bspl_dft = ff['coef_bspl_dft']

    n_nmodes = np.shape(coef_bspl_dft)[1]  # number of toroidal modes
    ns       = np.shape(coef_bspl_dft)[2]  # number of radial points
    m_width  = np.shape(coef_bspl_dft)[3]  # number of poloidal modes for every n-mode

    # --- find absolute values of poloidal modes ---
    m_offsets = np.meshgrid(np.arange(m_width), np.arange(ns))[0]
    m_modes   = np.zeros([n_nmodes, ns, m_width])
    for id_n_mode in np.arange(n_nmodes):
        loc_m_mins = np.array([ff['m_mins'][id_n_mode]]).transpose()
        m_modes[id_n_mode, :, :] = loc_m_mins + m_offsets

    # --- get toroidal mode numbers ---
    # here, len(n) must be equal to n_nmodes
    n_mask = None
    if ff['structure_name'] == 'antenna-init' or ff['structure_name'] == 'antenna':
        n_mask = ff['n_saved']
    n, ids_n, ids_n_mask = get_n_modes(d3['n_sim'], d3['n_sim_flux'], n_mask)

    # --- Take into account difference in field3d structures ---
    if ff['structure_name'] == 'antenna-init':
        # for antenna-init, only antenna toroidal modes are saved
        ids_n_res = range(len(n))
    else:
        # for 3d diagnostic and antenna (in orb5_res_parallel.h5), all toroidal modes are saved
        ids_n_res = np.array(ids_n)

    # --- get poloidal and toroidal Fourier transformation of signal, but still Bspline in radial direction ---
    s_bspl_ft = np.array(coef_bspl_dft)
    for count_n, ntor in enumerate(n):
        m_modes_n1 = m_modes[ids_n_res[count_n], :, :]  # [s, m_width]
        s_bspl_ft[:, ids_n_res[count_n], :, :] *= \
            ( np.sinc(ntor/d3['nphi']) * np.sinc(m_modes_n1/d3['nchi']) )**(d3['nidbas']+1) * \
                np.exp(1j * np.pi * (d3['nidbas'] - 1) * (m_modes_n1 / d3['nchi'] + ntor / d3['nphi']))

    # Include negative n-modes:
    for count_n, ntor in enumerate(n):
        if ntor > 0:
            s_bspl_ft[:, ids_n_res[count_n], :, :] = 2 * s_bspl_ft[:, ids_n_res[count_n], :, :]

    # --- sum radial Bspline in radial direction ---
    data_ft = sum_rad_bspline(sel_opt_bspl_s, d3, ff, s_bspl_ft)

    # results:
    ff.update({
        'm': m_modes,  # [n, s, m_width]
        's': d3['s'], 'chi': d3['chi'], 'phi': d3['phi'],
        'data_ft': data_ft,    # [t, n, s, m]
    })


# def form_fields_nm_old(d3, ff, sel_opt_bspl_s='opt2'):
#     coef_bspl_dft = ff['coef_bspl_dft']
#
#     n_nmodes = np.shape(coef_bspl_dft)[1]  # number of toroidal modes
#     ns       = np.shape(coef_bspl_dft)[2]  # number of radial points
#     m_width  = np.shape(coef_bspl_dft)[3]  # number of poloidal modes for every n-mode
#
#     # --- find absolute values of poloidal modes ---
#     m_offsets = np.meshgrid(np.arange(m_width), np.arange(ns))[0]
#     m_modes   = np.zeros([n_nmodes, ns, m_width])
#     for id_n_mode in np.arange(n_nmodes):
#         loc_m_mins = np.array([ff['m_mins'][id_n_mode]]).transpose()
#         m_modes[id_n_mode, :, :] = loc_m_mins + m_offsets
#
#     # --- get toroidal mode numbers ---
#     # here, len(n_global) must be equal to n_nmodes
#     n_global, ids_n_global = get_n_modes(d3['n_sim'], d3['n_sim_flux'])
#
#     # --- get poloidal and toroidal Fourier transformation of signal, but still Bspline in radial direction ---
#     s_bspl_ft = np.array(coef_bspl_dft)
#     for count_n, ntor in enumerate(n_global):
#         m_modes_n1 = m_modes[ids_n_global[count_n], :, :]  # [s, m_width]
#         s_bspl_ft[:, ids_n_global[count_n], :, :] *= \
#             ( np.sinc(ntor/d3['nphi']) * np.sinc(m_modes_n1/d3['nchi']) )**(d3['nidbas']+1) * \
#                 np.exp(1j * np.pi * (d3['nidbas'] - 1) * (m_modes_n1 / d3['nchi'] + ntor / d3['nphi']))
#
#     # Include negative n-modes:
#     for count_n, ntor in enumerate(n_global):
#         if ntor > 0:
#             s_bspl_ft[:, ids_n_global[count_n], :, :] = 2 * s_bspl_ft[:, ids_n_global[count_n], :, :]
#
#     # --- sum radial Bspline in radial direction ---
#     data_ft = sum_rad_bspline(sel_opt_bspl_s, d3, ff, s_bspl_ft)
#
#     # results:
#     ff.update({
#         'm': m_modes,  # [n, s, m_width]
#         's': d3['s'], 'chi': d3['chi'], 'phi': d3['phi'],
#         'data_ft': data_ft,    # [t, n, s, m]
#     })
#
#
# def form_antenna_nm(d3, ff, sel_opt_bspl_s='opt2'):
#     coef_bspl_dft = ff['coef_bspl_dft']
#
#     n_nmodes = np.shape(coef_bspl_dft)[1]  # number of toroidal modes
#     ns       = np.shape(coef_bspl_dft)[2]  # number of radial points
#     m_width  = np.shape(coef_bspl_dft)[3]  # number of poloidal modes for every n-mode
#
#     # --- find absolute values of poloidal modes ---
#     m_offsets = np.meshgrid(np.arange(m_width), np.arange(ns))[0]
#     m_modes   = np.zeros([n_nmodes, ns, m_width])
#     for id_n_mode in np.arange(n_nmodes):
#         loc_m_mins = np.array([ff['m_mins'][id_n_mode]]).transpose()
#         m_modes[id_n_mode, :, :] = loc_m_mins + m_offsets
#
#     # --- get antenna toroidal modes, which are actually simulated: ---
#     # len(n_ant) must be equal to n_nmodes:
#     n_ant, ids_n_ant, ids_n_ant_global = \
#         get_antenna_n_modes(d3['n_sim'], d3['n_sim_flux'], ff['n_save'])
#
#     # --- get poloidal and toroidal Fourier transformation of signal, but still Bspline in radial direction ---
#     s_bspl_ft = np.array(coef_bspl_dft)
#     for count_n, ntor in enumerate(n_ant):
#         m_modes_n1 = m_modes[ids_n_ant_global[count_n], :, :]  # [s, m_width]
#         s_bspl_ft[:, ids_n_ant_global[count_n], :, :] *= \
#             ( np.sinc(ntor/d3['nphi']) * np.sinc(m_modes_n1/d3['nchi']) )**(d3['nidbas']+1) * \
#                 np.exp(1j * np.pi * (d3['nidbas'] - 1) * (m_modes_n1 / d3['nchi'] + ntor / d3['nphi']))
#
#     # Include negative n-modes:
#     for count_n, ntor in enumerate(n_ant):
#         if ntor > 0:
#             s_bspl_ft[:, ids_n_ant_global[count_n], :, :] = 2 * s_bspl_ft[:, ids_n_ant_global[count_n], :, :]
#
#     # --- sum radial Bspline in radial direction ---
#     data_ft = sum_rad_bspline(sel_opt_bspl_s, d3, ff, s_bspl_ft)
#
#     # results:
#     ff.update({
#         'm': m_modes,  # [n, s, m_width]
#         's': d3['s'], 'chi': d3['chi'], 'phi': d3['phi'],
#         'data_ft': data_ft,    # [t, n, s, m]
#     })
#
#
# def form_antenna_init_nm(d3, ff, sel_opt_bspl_s='opt2'):
#     coef_bspl_dft = ff['coef_bspl_dft']
#
#     n_nmodes = np.shape(coef_bspl_dft)[1]  # number of toroidal modes
#     ns       = np.shape(coef_bspl_dft)[2]  # number of radial points
#     m_width  = np.shape(coef_bspl_dft)[3]  # number of poloidal modes for every n-mode
#
#     s = d3['s']
#
#     # antenna toroidal modes, which are actually simulated:
#     # len(n_ant) must be equal to n_nmodes:
#     n_ant, ids_n_ant, ids_n_ant_global = \
#         get_antenna_n_modes(d3['n_sim'], d3['n_sim_flux'], ff['n_save'])
#
#     # --- find absolute values of poloidal modes ---
#     m_offsets = np.meshgrid(np.arange(m_width), np.arange(ns))[0]
#     m_modes   = np.zeros([n_nmodes, ns, m_width])
#     for id_n_mode in np.arange(n_nmodes):
#         loc_m_mins = np.array([ff['m_mins'][id_n_mode]]).transpose()
#         m_modes[id_n_mode, :, :] = loc_m_mins + m_offsets
#
#     # --- get poloidal and toroidal Fourier transformation of signal, but still Bspline in radial direction ---
#     s_bspl_ft = np.array(coef_bspl_dft)
#     for count_n, ntor in enumerate(n_ant):
#         m_modes_n1 = m_modes[count_n, :, :]  # [s, m_width]
#         s_bspl_ft[:, count_n, :, :] *= \
#             ( np.sinc(ntor/d3['nphi']) * np.sinc(m_modes_n1/d3['nchi']) )**(d3['nidbas']+1) * \
#                 np.exp(1j * np.pi * (d3['nidbas'] - 1) * (m_modes_n1 / d3['nchi'] + ntor / d3['nphi']))
#
#     # Include negative n-modes:
#     for count_n, ntor in enumerate(n_ant):
#         if ntor > 0:
#             s_bspl_ft[:, count_n, :, :] = 2 * s_bspl_ft[:, count_n, :, :]
#
#     # --- sum radial Bspline in radial direction ---
#     data_ft = sum_rad_bspline(sel_opt_bspl_s, d3, ff, s_bspl_ft)
#
#     # results:
#     ff.update({
#         'm_modes': m_modes,  # [n, s, m_width]
#         's': s, 'chi': d3['chi'], 'phi': d3['phi'],
#         'data_ft': data_ft,    # [t, n, s, m]
#     })


def signal_chi1_phi1(d3, ff, chi_point, phi_point=0.0):
    data_ft = ff['data_ft']  # [t, n, s, m]

    s = ff['s']
    t = ff['t']
    m_modes  = ff['m']  # [n, s, m_width]

    # field toroidal modes
    n, ids_n, _ = get_n_modes(d3['n_sim'], d3['n_sim_flux'])

    # Get a signal in real poloidal and toroidal coordinates:
    signal_ts = np.zeros([np.size(t), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        for count_n, ntor in enumerate(n):
            id_n = ids_n[count_n]
            for id_m, m_mode in enumerate(m_modes[id_n, id_s1, :]):
                signal_ts[:, id_s1] += data_ft[:, id_n, id_s1, id_m] \
                            * np.exp(1j * m_mode * chi_point + 1j * ntor * phi_point)
    return np.real(signal_ts)


def signal_chi1_phi1_antenna(d3, ff, chi_point, phi_point=0.0):
    data_ft = ff['data_ft']  # [t, n, s, m]

    s = ff['s']
    t = ff['t']
    m_modes = ff['m']  # [n, s, m_width]
    w = ff['w']  # antenna frequencies

    # field toroidal modes
    n, ids_n, ids_n_ant = get_n_modes(d3['n_sim'], d3['n_sim_flux'], ff['n_saved'])

    # --- Take into account difference in field3d structures ---
    if ff['structure_name'] == 'antenna-init':
        # for antenna-init, only antenna toroidal modes are saved
        ids_n_res = range(len(n))
    else:
        # for 3d diagnostic and antenna (in orb5_res_parallel.h5), all toroidal modes are saved
        ids_n_res = np.array(ids_n)

    # # Get a signal in real poloidal and toroidal coordinates:
    # signal_ts = np.zeros([np.size(t), np.size(s)], dtype=np.complex64)
    # for id_s1, s1 in enumerate(s):
    #     for count_n, ntor in enumerate(n):
    #         id_n = ids_n_res[count_n]
    #         id_n_ant = ids_n_ant[count_n]
    #         for id_m, m_mode in enumerate(m_modes[id_n, id_s1, :]):
    #             signal_ts[:, id_s1] += \
    #                 ( data_ft[0, id_n, id_s1, id_m] * np.cos(w[id_n_ant] * t) +
    #                   data_ft[1, id_n, id_s1, id_m] * np.sin(w[id_n_ant] * t) ) * \
    #                             np.cos(m_mode * chi_point + ntor * phi_point)
    #                           # np.exp(1j * m_mode * chi_point + 1j * ntor * phi_point)

    # Get a signal in real poloidal and toroidal coordinates:
    signal_ts = np.zeros([np.size(t), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        for count_n, ntor in enumerate(n):
            id_n = ids_n_res[count_n]
            id_n_ant = ids_n_ant[count_n]
            for id_m, m_mode in enumerate(m_modes[id_n, id_s1, :]):
                signal_ts[:, id_s1] += data_ft[0, id_n, id_s1, id_m] * np.cos(w[id_n_ant] * t) *\
                        np.exp(1j * m_mode * chi_point + 1j * ntor * phi_point)
    return np.real(signal_ts)


def signal_t1_phi1(d3, ff, t_point, phi_point=0.0):
    data_ft = ff['data_ft']  # [t, n, s, m]

    t = ff['t']
    s = ff['s']
    chi = ff['chi']
    m_modes = ff['m']  # [n, s, m_width]

    id_t1, t1, _ = mix.get_ids(t, t_point)

    # field toroidal modes
    n, ids_n, _ = get_n_modes(d3['n_sim'], d3['n_sim_flux'])

    # Get a signal in real poloidal and toroidal coordinates:
    signal_chis = np.zeros([np.size(chi), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        for count_n, ntor in enumerate(n):
            id_n = ids_n[count_n]
            for id_m, m_mode in enumerate(m_modes[id_n, id_s1, :]):
                signal_chis[:, id_s1] += data_ft[id_t1, id_n, id_s1, id_m] \
                        * np.exp(1j * m_mode * chi) * np.exp(1j * ntor * phi_point)
    return np.real(signal_chis), t1


def signal_n_allm(d3, ff, n_mode_chosen, chi_point):
    data_ft = ff['data_ft']  # [t, n, s, m]

    s = ff['s']
    t = ff['t']
    m_modes = ff['m']       # [n, s, m_width]

    # field toroidal modes
    n, ids_n, _ = get_n_modes(d3['n_sim'], d3['n_sim_flux'])

    id_n_chosen = np.argwhere(n == n_mode_chosen)[0][0]
    id_n_global = ids_n[id_n_chosen]

    # Get a signal in real poloidal and toroidal coordinates:
    signal_ts = np.zeros([np.size(t), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        for id_m, m_mode in enumerate(m_modes[id_n_global, id_s1, :]):
            signal_ts[:, id_s1] += data_ft[:, id_n_global, id_s1, id_m] \
                        * np.exp(1j * m_mode * chi_point)
    return np.real(signal_ts)


def choose_one_var_ts(one_signal):
    dd = one_signal['dd']
    opt_var = one_signal['variable']
    vvar, tit_var = None, ''
    s, t = None, None
    res = {}

    if opt_var == 'n1':
        n_mode_chosen = one_signal['n1']
        chi_point = one_signal.get('chi-point', 0.0)

        init(dd)
        ff = dd['3d']['fields']
        s = ff['s']
        t = ff['t']

        vvar = signal_n_allm(dd['3d'], ff, n_mode_chosen, chi_point)
        tit_var = '3D:\ <\Phi' + '(n = {:d})'.format(n_mode_chosen) + '>_m'
    if opt_var == 'potsc':
        chi_point = one_signal.get('chi-point', 0.0)
        phi_point = one_signal.get('phi-point', 0.0)

        init(dd)
        ff = dd['3d']['fields']
        s = ff['s']
        t = ff['t']

        vvar = signal_chi1_phi1(dd['3d'], ff, chi_point, phi_point)
        line_chi = '\chi = {:0.1f},\ '.format(chi_point)
        line_phi = '\\varphi = {:0.1f}'.format(phi_point)
        tit_var = '3D:\ \Phi' + '_{' + line_chi + line_phi + '}'
    if opt_var == 'potsc-antenna':
        chi_point = one_signal.get('chi-point', 0.0)
        phi_point = one_signal.get('phi-point', 0.0)

        structure_name = 'antenna'

        init_antenna(dd, structure_name=structure_name)
        ff = dd['3d'][structure_name]
        s = ff['s']
        t = ff['t']

        vvar = signal_chi1_phi1_antenna(dd['3d'], ff, chi_point, phi_point)
        line_chi = '\chi = {:0.1f},\ '.format(chi_point)
        line_phi = '\\varphi = {:0.1f}'.format(phi_point)
        tit_var = '3D:\ Antenna:\ \Phi' + '_{' + line_chi + line_phi + '}'
    if opt_var == 'potsc-antenna-init':
        chi_point = one_signal.get('chi-point', 0.0)
        phi_point = one_signal.get('phi-point', 0.0)
        file_name = one_signal('file_name', 'user_antenna.h5')
        flag_orb_shape = one_signal('flag_orb_shape', False)

        structure_name = 'antenna-init'
        init_antenna(dd, name_antenna_file=file_name, structure_name=structure_name,
                     flag_orb_shape=flag_orb_shape)
        ff = dd['3d'][structure_name]
        s = ff['s']
        t = ff['t']

        vvar = signal_chi1_phi1_antenna(dd['3d'], ff, chi_point, phi_point)
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


def choose_one_var_rz(one_signal):
    dd = one_signal['dd']
    opt_var = one_signal['variable']

    vvar, s, chi, tit_var = None, None, None, ''
    res = {}

    if opt_var == 'potsc':
        t_point = one_signal['t-point']
        phi_point = one_signal.get('phi-point', 0.0)

        init(dd)
        ff = dd['3d']['fields']
        s   = ff['s']
        chi = ff['chi']

        vvar, t1 = signal_t1_phi1(dd['3d'], ff, t_point, phi_point)
        line_time = 't = {:0.3e},\ '.format(t1)
        line_phi = '\\varphi = {:0.1f}'.format(phi_point)
        tit_var = '3D:\ \Phi' + '_{' + line_time + line_phi + '}'
    if opt_var == 'potsc-antenna':
        t_point = one_signal['t-point']
        phi_point = one_signal.get('phi-point', 0.0)

        structure_name = 'antenna'

        init_antenna(dd['3d'], dd, structure_name=structure_name)
        ff = dd['3d'][structure_name]
        s = ff['s']
        chi = ff['chi']

        vvar, t1 = signal_t1_phi1(ff, t_point, phi_point)
        line_time = 't = {:0.3e},\ '.format(t1)
        line_phi = '\\varphi = {:0.1f}'.format(phi_point)
        tit_var = '3D:\ Antenna:\ \Phi' + '_{' + line_time + line_phi + '}'
    if opt_var == 'potsc-antenna-init':
        t_point = one_signal['t-point']
        phi_point = one_signal.get('phi-point', 0.0)
        file_name = one_signal('file_name', 'user_antenna.h5')
        flag_orb_shape = one_signal('flag_orb_shape', False)

        structure_name = 'antenna-init'

        init_antenna(dd, name_antenna_file=file_name, structure_name=structure_name,
                     flag_orb_shape=flag_orb_shape)
        ff = dd['3d'][structure_name]
        s   = ff['s']
        chi = ff['chi']

        vvar, t1 = signal_t1_phi1(dd['3d'], ff, t_point, phi_point)
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

