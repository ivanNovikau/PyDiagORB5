import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import zf_gam as zf
import write_data as wr
import transport
import write_data
import Global_variables as GLO
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
    mix.reload_module(GLO)


# read 3D data from orb5_res_parallel
def init(dd, nn=None, flag_antenna=False):
    # nn - [...], n-modes to read, if None, read all n-modes

    # --- read necessary additional parameters ---
    read_add(dd)

    # --- create 3d.fields field ---
    d3 = dd['3d']
    if 'fields' not in d3:
        d3['fields'] = {}
    ff = dd['3d']['fields']
    ff['structure_name'] = ff.get('structure_name', 'fields')

    # --- n-modes to read ---
    if nn is None:
        nn_res = np.array(d3['n_all_possible'])
    else:
        nn_res = np.array(nn)
    nn_res = np.array([
        one_n_mode
        for one_n_mode in nn_res
        if d3['nfilt1'] <= one_n_mode <= d3['nfilt2']
    ])

    if len(nn_res) == 0:
        mix.error_mes('--- Wrong n mode. ---')

    ids_nn_res = list(nn_res - d3['nfilt1'])

    ff['n'] = nn_res

    # --- read 3d-field data ---
    path_to_file = dd['path'] + '/orb5_res_parallel.h5'
    f = h5.File(path_to_file, 'r')

    ff['t'] = ff.get(
        't',
        np.array(f['/data/var3d/generic/pot3d/time'])
    )
    ff['m_mins']  = np.array(
        f['/data/var3d/generic/pot3d/mmin'][ids_nn_res, :]
    )  # [n, s]

    if not flag_antenna:
        coef_bspl_dft = np.array(
            f['/data/var3d/generic/pot3d/data'][:, ids_nn_res, :, :]
        )  # [t, n, s, m]
    else:
        coef_bspl_dft = np.array(
            f['/data/var3d/generic/pot3d_antenna/data'][:, ids_nn_res, :, :]
        )  # [t, n, s, m]

    ff['coef_bspl_dft'] = coef_bspl_dft[:]['real'] + 1j * coef_bspl_dft[:]['imaginary']

    f.close()

    # --- create poloidal and toroidal Fourier transform of the field ---
    form_fields_nm(d3, ff)


# read 3D density data from orb5_res_parallel
def init_den(dd, name_species, nn=None):
    # nn - [...], n-modes to read, if None, read all n-modes

    # --- read necessary additional parameters ---
    read_add(dd)

    # --- create 3d.fields field ---
    d3 = dd['3d']
    name_structure = 'density-' + name_species
    if name_structure not in d3:
        d3[name_structure] = {}
    ff = dd['3d'][name_structure]
    ff['structure_name'] = ff.get('structure_name', name_structure)

    # --- n-modes to read ---
    if nn is None:
        nn_res = np.array(d3['n_all_possible'])
    else:
        nn_res = np.array(nn)
    nn_res = np.array([
        one_n_mode
        for one_n_mode in nn_res
        if d3['nfilt1'] <= one_n_mode <= d3['nfilt2']
    ])

    if len(nn_res) == 0:
        mix.error_mes('--- Wrong n mode. ---')

    ids_nn_res = list(nn_res - d3['nfilt1'])

    ff['n'] = nn_res

    # --- read 3d-field data ---
    path_to_file = dd['path'] + '/orb5_res_parallel.h5'
    f = h5.File(path_to_file, 'r')

    root_path = '/data/var3d/' + name_species + '/rho/'

    ff['t'] = ff.get('t', np.array(f[root_path + 'time']))
    ff['m_mins']  = np.array(f[root_path + 'mmin'][ids_nn_res, :])  # [n, s]
    coef_bspl_dft = np.array(f[root_path + 'data'][:, ids_nn_res, :, :])  # [t, n, s, m]
    ff['coef_bspl_dft'] = coef_bspl_dft[:]['real'] + 1j * coef_bspl_dft[:]['imaginary']

    f.close()

    # --- create poloidal and toroidal Fourier transform of the field ---
    form_fields_nm(d3, ff)


# read coordinate grids, read filter parameters
def read_add(dd):
    if 'n_all_possible_flux' in dd['3d']:
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

    # toroidal grid
    phi_grid = np.array(f['/data/var2d/generic/potsp/coord2'])

    # step of the flux-tube filter
    n_flux_tube = np.array(f['/parameters/fields/n_flux_tube'])[0].astype(np.int)

    # radial boundaries
    sfmin = np.array(f['/parameters/fields/sfmin'])[0]
    sfmax = np.array(f['/parameters/fields/sfmax'])[0]

    # close file
    f.close()

    # radial, poloidal grids:
    rd.potsc_grids(dd)
    s   = dd['potsc_grids']['s']
    chi = dd['potsc_grids']['chi']
    r = dd['potsc_grids']['r']
    z = dd['potsc_grids']['z']

    # all toroidal modes allowed by the imposed filter
    n_modes = [i for i in range(nfilt1, nfilt2 + 1)]

    # take into account step on the toroidal modes
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
        'n_all_possible': n_modes,
        'n_all_possible_flux': n_modes_flux,
    })


def form_fields_nm(d3, ff, sel_opt_bspl_s='opt2'):
    coef_bspl_dft = ff['coef_bspl_dft']

    n_nmodes = np.shape(coef_bspl_dft)[1]  # number of toroidal modes
    ns       = np.shape(coef_bspl_dft)[2]  # number of radial points
    m_width  = np.shape(coef_bspl_dft)[3]  # number of poloidal modes for every n-mode

    # --- create shifts of the poloidal number ---
    # shape(ns, m_width):
    # [
    #   [0, 1, ..., m_width-1],
    #   [0, 1, ..., m_width-1]
    # ]
    m_offsets = np.meshgrid(np.arange(m_width), np.arange(ns))[0]

    # --- find absolute values of poloidal modes ---
    m_modes   = np.zeros([n_nmodes, ns, m_width])
    for id_n_mode in np.arange(n_nmodes):
        loc_m_mins = np.array([ff['m_mins'][id_n_mode]]).transpose()
        m_modes[id_n_mode, :, :] = loc_m_mins + m_offsets

    # --- get toroidal mode numbers ---
    # here, len(n) must be equal to n_nmodes
    n, ids_n, _ = get_n_modes(ff['n'], d3['n_all_possible_flux'])

    # --- Take into account difference in field3d structures ---
    if ff['structure_name'] == 'antenna-init':
        # for antenna-init, only antenna toroidal modes are saved
        ids_n_res = range(len(n))
    else:
        ids_n_res = np.array(ids_n)

    # --- get poloidal and toroidal Fourier transformation of signal
    s_bspl_ft = np.array(coef_bspl_dft)  # [t, n, s, m]
    for count_n, ntor in enumerate(n):
        m_modes_n1 = m_modes[ids_n_res[count_n], :, :]  # [s, m_width]
        s_bspl_ft[:, ids_n_res[count_n], :, :] *= \
            ( np.sinc(ntor/d3['nphi']) *
              np.sinc(GLO.DEF_SIGN_M * m_modes_n1/d3['nchi']) )**(d3['nidbas']+1) * \
                np.exp(1j * np.pi * (d3['nidbas'] - 1) *
                       (GLO.DEF_SIGN_M * m_modes_n1 / d3['nchi'] + ntor / d3['nphi']))

    # Include negative n-modes:
    for count_n, ntor in enumerate(n):
        if ntor > 0:
            s_bspl_ft[:, ids_n_res[count_n], :, :] = 2 * s_bspl_ft[:, ids_n_res[count_n], :, :]

    data_ft = s_bspl_ft
    # # --- sum radial Bspline in radial direction ---
    # data_ft = sum_rad_bspline(sel_opt_bspl_s, d3, ff, s_bspl_ft)

    # results:
    ff.update({
        'm': m_modes,  # [n, s, m_width]
        's': d3['s'], 'chi': d3['chi'], 'phi': d3['phi'],
        'data_ft': data_ft,    # [t, n, s, m]
    })


# get n-modes taking into account a step and a mask;
# return also ids of the result modes within initial array of modes
#   and within an imposed mask
def get_n_modes(n_modes, n_modes_flux, n_mask=None):
    # n_modes - all available modes
    # n_modes_flux - all n-modes, where a finite step between the modes
    #                   is taken into account
    # n_maks - modes, which are allowed to be taken into account
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


# Save antenna structure into a .h5 file
def save_rhsf(oo):
    # oo = {} with the following fields:
    # n_save - array of toroidal modes to save.
    # flag_one_freq - True: every antenna toroidal mode can have only one frequency;
    # freqs_ant_wc:
    #   if flag_one_freq = True: 1d array of frequencies for every toroidal mode n_save.
    #   if flag_one_freq = False: 2d array of frequencies for every m and n antenna modes.
    # gamma_ant_wc - array of damping (growth) rates of toroidal modes n_save.
    # A_scaling - array of scalings of the result structures for every toroidal mode.
    # name_file_antenna - name of a file where one is going to save
    #                       the antenna's radial structure
    # t_moment - when to save 3d structure, if None, when save it at the end of the simulation
    # m_mask - to indicate which m-modes are taken into account in antenna;
    #   the same mask is applied for every n mode:  [-deltam, ..., -1, 0, 1,... deltam] or None;
    #   if None, all m-modes are applied [-deltam,...deltam]

    def save_coef_bspl_dft(ids_t, ids_n_ant_global_local, gammas=None, T_periods=None):
        if gammas is None:
            gammas = np.zeros(len(n_save_res))
            T_periods = np.zeros(len(n_save_res))

        coef_bspl_dft = np.zeros([len(n_save_res), ns_bspl, nm], dtype=np.complex128)  # [n, s, m]
        for id_n_mode in range(len(n_save_res)):
            coef_bspl_dft[id_n_mode, :, :] = \
                ff['coef_bspl_dft'][ids_t[id_n_mode], ids_n_ant_global_local[id_n_mode], :, :] * \
                A_scaling_res[id_n_mode] \
                * np.exp(gammas[id_n_mode] * coef_quarter * T_periods[id_n_mode])
        return coef_bspl_dft

    # get parameters
    dd = oo['dd']
    n_save = oo['n_save']
    freqs_ant_wc = oo['freqs_ant_wc']
    gamma_ant_wc = oo['gamma_ant_wc']

    A_scaling = oo.get('A_scaling', None)
    name_file_antenna = oo.get('name_file_antenna', 'user_antenna.h5')
    t_moments = oo.get('t_moments', None)

    # --- read 3d-field data ---
    init(dd, nn=n_save)
    d3 = dd['3d']
    ff = d3['fields']
    t = ff['t']

    mmins   = ff['m_mins']  # [n, s]
    ns_bspl = np.shape(mmins)[1]
    nm = 2 * d3['deltam'] + 1

    # antenna toroidal modes to save: take into account that not every n-mode in
    # n_save array might be present in the considered simulation (in the array ff['n'])
    n_save_res, ids_n_ant_global, ids_n_ant = get_n_modes(ff['n'], d3['n_all_possible_flux'], n_save)

    # time moments when to save 3d structure:
    t_moments_res, ids_t_antenna = \
        np.zeros(len(n_save_res)), np.zeros(len(n_save_res), dtype=np.int)
    for id_n, t_moment in enumerate(t_moments):
        if t_moments is not None:
            ids_t_antenna[id_n], t_moments_res[id_n], _ = mix.get_ids(t, t_moment)
        else:
            ids_t_antenna[id_n], t_moments_res[id_n] = -1, t[-1]
        # print('Save n={:d} antenna structure at t = {:0.1e}'.format(
        #     n_save_res[id_n], t_moments_res[id_n])
        # )

    # --- minimum poloidal modes ---
    mmins_res = np.zeros( [len(n_save_res), np.shape(mmins)[1]] )
    for id_n1_ant, _ in enumerate(n_save_res):
        mmins_res[id_n1_ant, :] = mmins[ids_n_ant_global[id_n1_ant], :]

    # --- result antenna frequencies and periods, and damping/growth rates ---
    freqs_ant_res = np.zeros(len(n_save_res))
    gamma_ant_res = np.zeros(len(n_save_res))
    for counter_loc, id_ant in enumerate(ids_n_ant):
        freqs_ant_res[counter_loc] = freqs_ant_wc[id_ant]
        gamma_ant_res[counter_loc] = gamma_ant_wc[id_ant]
    T_periods_res = 2 * np.pi / abs(freqs_ant_res)

    # where we need an additional signal:
    coef_quarter = 0.25  # to save a signal at (1 - coef_quarter) * T

    # time points of the quarter of antenna toroidal mode's periods:
    ids_t_quarters, ts_quarters = \
        np.zeros(len(n_save_res), dtype=np.int), np.zeros(len(n_save_res))
    for id_n, T_period in enumerate(T_periods_res):
        ids_t_quarters[id_n], ts_quarters[id_n], _ = \
            mix.get_ids(t, t_moments_res[id_n] - coef_quarter * T_period)

    # antenna scaling
    A_scaling_res = np.zeros(len(n_save_res))
    if A_scaling is None:
        A_scaling_res = np.ones(len(n_save_res))
    else:
        for counter_loc, id_ant in enumerate(ids_n_ant):
            A_scaling_res[counter_loc] = A_scaling[id_ant]

    # print some antenna parameters:
    for id_n, n1 in enumerate(n_save_res):
        print('ANT n={:d}: t = {:0.1e}, A = {:0.1f}, w[wc] = {:0.2e}, g[wc] = {:0.2e}'.format(
                n_save_res[id_n], t_moments_res[id_n], A_scaling_res[id_n],
                freqs_ant_res[id_n], gamma_ant_res[id_n]
            )
        )

    # radial structure at last time point
    coef_bspl_dft_BEGIN = save_coef_bspl_dft(ids_t_antenna, ids_n_ant_global)

    # radial structure at a quarter of the n1 mode's period
    coef_bspl_dft_QUARTER = save_coef_bspl_dft(ids_t_quarters, ids_n_ant_global,
                                               gamma_ant_res, T_periods_res)

    # --- SAVE ANTENNA'S STRUCTURE ---
    path_to_file = dd['path'] + '/' + name_file_antenna
    ffile = h5.File(path_to_file, 'w')
    try:
        # save toroidal modes, their frequencies and scalings
        grp = ffile.create_group("parameters")

        ddata = grp.create_dataset('n_saved', data=n_save_res)
        ddata.attrs[u'descr'] = 'Antenna toroidal modes.'

        ddata = grp.create_dataset('w', data=np.array(freqs_ant_res), dtype=np.float)
        ddata.attrs[u'descr'] = 'Frequencies [norm. to wci] of antenna toroidal modes.'

        ddata = grp.create_dataset('A', data=np.array(A_scaling), dtype=np.float)
        ddata.attrs[u'descr'] = 'Scales of antenna toroidal modes.'

        # save radial structure and time moments
        grp = ffile.create_group("data/var3d/generic/phi_antenna/run.1")

        ddata = grp.create_dataset('t0', data=t_moments_res)
        ddata.attrs[u'descr'] = \
            'Time moments, where initial antenna structures for different n-modes are saved.'

        ddata = grp.create_dataset('tq', data=ts_quarters)
        ddata.attrs[u'descr'] = 'Actual time moments at three quarters of n-mode periods, ' \
                                'where antenna structures for different n-modes are saved.'

        tq_theory = t_moments_res - coef_quarter * T_periods_res
        ddata = grp.create_dataset('tq-theory', data=tq_theory)
        ddata.attrs[u'descr'] = 'Theoretical time moments at three quarters of n-mode periods, ' \
                                'where antenna structures for different n-modes are saved.'

        ddata = grp.create_dataset('mmin', data=mmins_res)
        ddata.attrs[u'descr'] = 'Minimal poloidal mode numbers.'

        ddata = grp.create_dataset('bspl_dft_t0', data=coef_bspl_dft_BEGIN)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients at t0'

        ddata = grp.create_dataset('bspl_dft_tquarter', data=coef_bspl_dft_QUARTER)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients at tq'

        # change the shape:
        nn_bspl = np.shape(coef_bspl_dft_BEGIN)[0]
        ns_bspl = np.shape(coef_bspl_dft_BEGIN)[1]
        wm_bspl = np.shape(coef_bspl_dft_BEGIN)[2]
        neq_f = ns_bspl*wm_bspl

        # adjust shape of the matrix to the shape used in ORB5
        new_bspl_dft_BEGIN      = coef_bspl_dft_BEGIN.reshape((nn_bspl, neq_f))
        new_bspl_dft_QUARTER    = coef_bspl_dft_QUARTER.reshape((nn_bspl, neq_f))

        # array of shape [nn, neq] (when ORB5 reads the .h5 file, it will automatically transpose the matrix)
        temp_t0 = np.zeros(np.shape(new_bspl_dft_BEGIN),   dtype=GLO.comp_datatype)
        temp_tq = np.zeros(np.shape(new_bspl_dft_QUARTER), dtype=GLO.comp_datatype)
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
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients at t0' \
                                'Shape is [(ns+nidbas)*width_m, number_of_n_modes]'

        ddata = grp.create_dataset('orb_bspl_dft_tquarter', data=temp_tq)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients at tq' \
                                'Shape is [(ns+nidbas)*width_m, number_of_n_modes]'
    except:
        ffile.close()
        sys.exit(-1)

    ffile.close()

    return


# Combine two antenna .h5 files
def merge_several_rhsf(full_file_names, res_path, res_name_file):
    # antennae have to have different n-modes
    # antennas have to be defined on the same radial space grids, and have the same deltam

    res_file = res_path + '/' + res_name_file

    nres = []
    na, wa, t0a, tqa, tq_theory_a = [], [], [], [], []
    mmin_init = []
    bspl_t0, bspl_tq, orb_bspl_t0, orb_bspl_tq = [], [], [], []
    for id_file, fname in enumerate(full_file_names):
        nloc = rd.read_array(fname, '/parameters/n_saved')

        nres += list(nloc)
        na.append(np.array(nloc, dtype=np.int))
        wa.append(rd.read_array(fname, '/parameters/w'))

        bline = 'data/var3d/generic/phi_antenna/run.1/'
        t0a.append(rd.read_array(fname, bline + 't0'))
        tqa.append(rd.read_array(fname, bline + 'tq'))
        tq_theory_a.append(rd.read_array(fname, bline + 'tq-theory'))

        mmin_init.append(rd.read_array(fname, bline + 'mmin'))

        bspl_t0.append(rd.read_array(fname, bline + 'bspl_dft_t0'))
        bspl_tq.append(rd.read_array(fname, bline + 'bspl_dft_tquarter'))
        orb_bspl_t0.append(rd.read_array(fname, bline + 'orb_bspl_dft_t0'))
        orb_bspl_tq.append(rd.read_array(fname, bline + 'orb_bspl_dft_tquarter'))

    # form a resulting array with antenna' toroidal mode numbers
    nres.sort()
    nres = np.array(nres, dtype=np.int)
    number_n_modes = len(nres)
    if number_n_modes != len(np.unique(nres)):
        print('Error: antenna files should contain different modes.')
        return

    ws_res = np.zeros(number_n_modes)
    t0_res = np.zeros(number_n_modes)
    tq_res = np.zeros(number_n_modes)
    tq_theory_res = np.zeros(number_n_modes)
    mmin_res = np.zeros([ number_n_modes, np.shape(mmin_init[0])[1] ])

    ns = np.shape(bspl_t0[0])[1]
    nm = np.shape(bspl_t0[0])[2]
    orb_shape = np.shape(orb_bspl_t0[0])
    bspl_t0_res = np.zeros([number_n_modes, ns,  nm], dtype=np.complex128)
    bspl_tq_res = np.zeros([number_n_modes, ns,  nm], dtype=np.complex128)
    orb_bspl_t0_res = np.zeros([number_n_modes, orb_shape[1]], dtype=GLO.comp_datatype)
    orb_bspl_tq_res = np.zeros([number_n_modes, orb_shape[1]], dtype=GLO.comp_datatype)
    for id_n, n1 in enumerate(nres):
        for id_file, n_file in enumerate(na):
            if n1 in n_file:
                id_n_pos = np.where(n_file == n1)[0][0]

                # form a resulting array with antenna modes' frequencies
                ws_res[id_n] = wa[id_file][id_n_pos]

                # form resulting arrays with time moments
                t0_res[id_n] = t0a[id_file][id_n_pos]
                tq_res[id_n] = tqa[id_file][id_n_pos]
                tq_theory_res[id_n] = tq_theory_a[id_file][id_n_pos]

                # form resulting matrix with minimal m-numbers
                mmin_res[id_n, :] = mmin_init[id_file][id_n_pos]

                # form resulting DFT of finite elements in shape (n, s, m)
                bspl_t0_res[id_n, :, :] = bspl_t0[id_file][id_n_pos]
                bspl_tq_res[id_n, :, :] = bspl_tq[id_file][id_n_pos]

                # form resulting DFT of finite elements in shape (n, s*m)
                orb_bspl_t0_res[id_n, :] = orb_bspl_t0[id_file][id_n_pos]
                orb_bspl_tq_res[id_n, :] = orb_bspl_tq[id_file][id_n_pos]

    # save the antenna parameters and structure in the resulting file:
    ffile = h5.File(res_file, 'w')
    try:
        # save toroidal modes, their frequencies and scalings
        grp = ffile.create_group("parameters")

        ddata = grp.create_dataset('n_saved', data=nres)
        ddata.attrs[u'descr'] = 'Antenna toroidal modes.'

        ddata = grp.create_dataset('w', data=ws_res, dtype=np.float)
        ddata.attrs[u'descr'] = 'Frequencies [norm. to wci] of antenna toroidal modes.'

        # save radial structure and time moments
        grp = ffile.create_group("data/var3d/generic/phi_antenna/run.1")

        ddata = grp.create_dataset('t0', data=t0_res)
        ddata.attrs[u'descr'] = \
            'Time moments, where initial antenna structures for different n-modes are saved.'

        ddata = grp.create_dataset('tq', data=tq_res)
        ddata.attrs[u'descr'] = 'Actual time moments at three quarters of n-mode periods, ' \
                                'where antenna structures for different n-modes are saved.'

        ddata = grp.create_dataset('tq-theory', data=tq_theory_res)
        ddata.attrs[u'descr'] = 'Theoretical time moments at three quarters of n-mode periods, ' \
                                'where antenna structures for different n-modes are saved.'

        ddata = grp.create_dataset('mmin', data=mmin_res)
        ddata.attrs[u'descr'] = 'Minimal poloidal mode numbers.'

        ddata = grp.create_dataset('bspl_dft_t0', data=bspl_t0_res)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients at t0'

        ddata = grp.create_dataset('bspl_dft_tquarter', data=bspl_tq_res)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients at tq'

        ddata = grp.create_dataset('orb_bspl_dft_t0', data=orb_bspl_t0_res)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients at t0' \
                                'Shape is [(ns+nidbas)*width_m, number_of_n_modes]'

        ddata = grp.create_dataset('orb_bspl_dft_tquarter', data=orb_bspl_tq_res)
        ddata.attrs[u'descr'] = 'Discrete Fourier Transform of Bspline coefficients at tq' \
                                'Shape is [(ns+nidbas)*width_m, number_of_n_modes]'
    except:
        ffile.close()
        sys.exit(-1)

    ffile.close()


def init_antenna(dd, name_antenna_file='', structure_name='antenna', flag_orb_shape=False):

    def read_coef_bspl(id_t):
        # id_t = 0 or 1
        name_bspl = 'bspl_dft_tquarter' if id_t == 1 else 'bspl_dft_t0'
        temp_bspl = np.array(f[path_antenna + name_bspl])[:, :, :]  # [dummy_t, n, s, m]

        if temp_bspl.dtype == np.complex128:
            ffa['coef_bspl_dft'][id_t, :, :, :] = temp_bspl
        else:
            ffa['coef_bspl_dft'][id_t, :, :, :] = \
                temp_bspl[:]['real'] + 1j * temp_bspl[:]['imaginary']

    def read_coef_bspl_orb_shape(id_t):
        # id_t = 0 or 1
        name_bspl = 'orb_bspl_dft_tquarter' if id_t == 1 else 'orb_bspl_dft_t0'
        orb_shape_bspl = np.array(f[path_antenna + name_bspl])[:, :]  # [ns * nm, nn]
        orb_shape_bspl = orb_shape_bspl[:]['real'] + 1j * orb_shape_bspl[:]['imaginary']
        ffa['coef_bspl_dft'][id_t, :, :, :] = orb_shape_bspl.reshape((nn, ns, nm))

    # --- read necessary additional parameters ---
    read_add(dd)

    # --- create antenna field ---
    d3 = dd['3d']
    d3[structure_name] = {}
    ffa = dd['3d'][structure_name]

    ffa['structure_name'] = structure_name

    # path to files:
    path_to_files = dd['path']

    # --- read input antenna parameters ---
    path_to_file = path_to_files + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')
    try:
        # antenna's type
        ffa['nsel_type'] = 'none'
        path_parameter = '/parameters/antenna/nsel_type'
        if path_parameter in f:
            antenna_type = f[path_parameter].attrs
            ids_attr = list(antenna_type)
            ffa['nsel_type'] = [antenna_type[name].decode("utf-8")
                for name in ids_attr[0:len(ids_attr) - 1]][2*(dd['n_starts'] - 1)]

        # type of antenna's radial profile
        ffa['nsel_profile'] = ''
        path_parameter = '/parameters/antenna/nsel_profile'
        if path_parameter in f:
            antenna_type_profile = f[path_parameter].attrs
            ids_attr = list(antenna_type_profile)
            ffa['nsel_profile'] = [antenna_type_profile[name].decode("utf-8")
                for name in ids_attr[0:len(ids_attr) - 1]][2*(dd['n_starts'] - 1)]

        # name of .h5 file with antenna initial structure
        ffa['file_name_rhsf'] = ''
        path_parameter = '/parameters/antenna/file_name_rhsf'
        if path_parameter in f:
            antenna_type = f[path_parameter].attrs
            ids_attr = list(antenna_type)
            ffa['file_name_rhsf'] = [antenna_type[name].decode("utf-8")
                for name in ids_attr[0:len(ids_attr) - 1]][dd['n_starts'] - 1]

        if ffa['nsel_type'] == 'none':
            ffa['file_name_rhsf'] = ''

        if structure_name == 'antenna-init':
            ffa['nsel_type'] = 'phi'
            ffa['nsel_profile'] = 'rhsf'
    except:
        f.close()
        sys.exit(-1)
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

    # read antenna parameters (from a file dedicated to initial antenna)
    if ffa['nsel_profile'] == 'rhsf':
        path_to_file = path_to_files + '/' + ffa['file_name_rhsf']
        finit = h5.File(path_to_file, 'r')
        try:
            path_parameter = '/parameters/w'
            if path_parameter in finit:
                ffa['w'] = np.array(finit[path_parameter])

            path_parameter = '/parameters/n_saved'
            if path_parameter in finit:
                ffa['n'] = np.array(finit[path_parameter])
        except:
            finit.close()
            sys.exit(-1)
        finit.close()

    # read antenna structure (from user_antenna.h5 or orb5_res_parallel.h5)
    path_to_file = path_to_files + '/' + name_antenna_file
    f = h5.File(path_to_file, 'r')
    try:
        line_id_start = dd['n_starts_string']
        if structure_name == 'antenna-init':
            line_id_start = '1'

        path_antenna  = '/data/var3d/generic/' + ffa['nsel_type'] + \
                        '_antenna/run.' + line_id_start + '/'
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

    # artificial time grid:
    # where
    # 0 - initial antenna structure,
    # 1 - antenna structure at a quarter of a n mode's period
    ffa['t'] = [0, 1]

    # --- create poloidal and toroidal Fourier transform of the field ---
    form_fields_nm(d3, ffa)


def signal_chi1_phi1(d3, ff, chi_point, phi_point=0.0):
    data_ft = ff['data_ft']  # [t, n, s, m]

    s = ff['s']
    t = ff['t']
    m_modes  = ff['m']  # [n, s, m_width]

    # field toroidal modes
    n, ids_n, _ = get_n_modes(ff['n'], d3['n_all_possible_flux'])

    # --- Take into account difference in fields3d structures ---
    ids_n_res = range(len(n)) \
        if ff['structure_name'] == 'antenna-init' \
        else np.array(ids_n)

    # Get a signal in real poloidal and toroidal coordinates:
    signal_ts = np.zeros([np.size(t), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        for count_n, ntor in enumerate(n):
            id_n = ids_n_res[count_n]
            for id_m, m_mode in enumerate(m_modes[id_n, id_s1, :]):
                signal_ts[:, id_s1] += data_ft[:, id_n, id_s1, id_m] \
                            * np.exp(GLO.DEF_SIGN_M * 1j * m_mode * chi_point + 1j * ntor * phi_point)
    return np.real(signal_ts)


def signal_t1_phi1(d3, ff, t_point, phi_point=0.0):
    data_ft = ff['data_ft']  # [t, n, s, m]

    t = ff['t']
    s = ff['s']
    chi = ff['chi']
    m_modes = ff['m']  # [n, s, m_width]

    id_t1, t1, _ = mix.get_ids(t, t_point)

    # field toroidal modes
    n, ids_n, _ = get_n_modes(ff['n'], d3['n_all_possible_flux'])

    # --- Take into account difference in field3d structures ---
    ids_n_res = range(len(n)) \
        if ff['structure_name'] == 'antenna-init' \
        else np.array(ids_n)

    # Get a signal in real poloidal and toroidal coordinates:
    signal_chis = np.zeros([np.size(chi), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        for count_n, ntor in enumerate(n):
            id_n = ids_n_res[count_n]
            for id_m, m_mode in enumerate(m_modes[id_n, id_s1, :]):
                signal_chis[:, id_s1] += data_ft[id_t1, id_n, id_s1, id_m] \
                        * np.exp(GLO.DEF_SIGN_M * 1j * m_mode * chi) * np.exp(1j * ntor * phi_point)
    return np.real(signal_chis), t1


def signal_n_allm_ts(d3, ff, n_mode_chosen, chi_point, sel_e='none'):
    data_ft = ff['data_ft']  # [t, n, s, m]

    s = ff['s']
    t = ff['t']
    chi = ff['chi']
    m_modes = ff['m']       # [n, s, m_width]

    # field toroidal modes
    n, ids_n, _ = get_n_modes(ff['n'], d3['n_all_possible_flux'], [n_mode_chosen])
    if len(n) == 0:
        mix.error_mes('--- Wrong n mode in a signal constrution. ---')
    id_n_res = ids_n[0]

    # Get a signal in real poloidal and toroidal coordinates:
    signal_ts = np.zeros([np.size(t), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        for id_m, m_mode in enumerate(m_modes[id_n_res, id_s1, :]):
            # signal_ts[:, id_s1] += data_ft[:, id_n_res, id_s1, id_m] \
            #                        * np.exp(GLO.DEF_SIGN_M * 1j * m_mode * chi_point)

            temp = data_ft[:, id_n_res, id_s1, id_m]

            if sel_e == 'echi':  # poloidal component of electric field: E_chi
                id_chi, chi_current, _ = mix.get_ids(chi, chi_point)

                temp = temp[:, None] * np.exp(GLO.DEF_SIGN_M * 1j * m_mode * chi[None, :])
                temp = - np.gradient(temp, chi, axis=1) / s1

                signal_ts[:, id_s1] += temp[:, id_chi]
            else:
                signal_ts[:, id_s1] += temp * np.exp(GLO.DEF_SIGN_M * 1j * m_mode * chi_point)
    return np.real(signal_ts)


def signal_n1_m1_ts(d3, ff, n_mode_chosen, m_mode_chosen, sel_e='none', chi_point=0):
    data_ft = ff['data_ft']  # [t, n, s, m]

    s = ff['s']
    t = ff['t']
    m_modes = ff['m']       # [n, s, m_width]
    chi = ff['chi']

    # field toroidal modes
    n, ids_n, _ = get_n_modes(ff['n'], d3['n_all_possible_flux'], [n_mode_chosen])
    if len(n) == 0:
        mix.error_mes('--- Wrong n mode in a signal construction. ---')
    id_n_res = ids_n[0]

    # Get a signal in real poloidal and toroidal coordinates:
    signal_ts = np.zeros([np.size(t), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        current_m_modes = m_modes[id_n_res, id_s1, :]
        id_m = mix.get_id_int_strict(current_m_modes, m_mode_chosen)

        temp = 0 if np.isnan(id_m) else data_ft[:, id_n_res, id_s1, id_m]

        if sel_e == 'echi':  # poloidal component of electric field: E_chi
            if np.isnan(id_m):
                signal_ts[:, id_s1] += 0
                continue

            id_chi, chi_current, _ = mix.get_ids(chi, chi_point)

            temp = temp[:, None] * np.exp(GLO.DEF_SIGN_M * 1j * m_mode_chosen * chi[None, :])
            temp = - np.gradient(temp, chi, axis=1) / s1

            signal_ts[:, id_s1] += temp[:, id_chi]
        else:
            signal_ts[:, id_s1] += temp
    return np.real(signal_ts)


def signal_n1_m1_t1(d3, ff, n_mode_chosen, m_mode_chosen, t_point):
    data_ft = ff['data_ft']  # [t, n, s, m]

    s = ff['s']
    t = ff['t']
    chi = ff['chi']
    m_modes = ff['m']  # [n, s, m_width]

    id_t1, t1, _ = mix.get_ids(t, t_point)

    # field toroidal modes
    n, ids_n, _ = get_n_modes(ff['n'], d3['n_all_possible_flux'], [n_mode_chosen])
    if len(n) == 0:
        mix.error_mes('--- Wrong n mode in a signal constrution. ---')
    id_n_res = ids_n[0]

    # Get a signal in real poloidal and toroidal coordinates:
    signal_chis = np.zeros([np.size(chi), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        current_m_modes = m_modes[id_n_res, id_s1, :]
        id_m, m_current, _ = mix.get_ids(current_m_modes, m_mode_chosen)
        signal_chis[:, id_s1] += data_ft[id_t1, id_n_res, id_s1, id_m] \
                                     * np.exp(GLO.DEF_SIGN_M * 1j * m_mode_chosen * chi)
    return np.real(signal_chis), t1


def signal_n_allm_t1(d3, ff, n_mode_chosen, t_point):
    data_ft = ff['data_ft']  # [t, n, s, m]

    s = ff['s']
    t = ff['t']
    chi = ff['chi']
    m_modes = ff['m']  # [n, s, m_width]

    id_t1, t1, _ = mix.get_ids(t, t_point)

    # field toroidal modes
    n, ids_n, _ = get_n_modes(ff['n'], d3['n_all_possible_flux'], [n_mode_chosen])
    if len(n) == 0:
        mix.error_mes('--- Wrong n mode in a signal constrution. ---')
    id_n_res = ids_n[0]

    # Get a signal in real poloidal and toroidal coordinates:
    signal_chis = np.zeros([np.size(chi), np.size(s)], dtype=np.complex64)
    for id_s1, s1 in enumerate(s):
        for id_m, m_mode in enumerate(m_modes[id_n_res, id_s1, :]):
            signal_chis[:, id_s1] += data_ft[id_t1, id_n_res, id_s1, id_m] \
                                     * np.exp(GLO.DEF_SIGN_M * 1j * m_mode * chi)
    return np.real(signal_chis), t1


def choose_one_var_ts(one_signal):
    dd = one_signal['dd']
    opt_var = one_signal['variable']
    vvar, tit_var = None, ''
    s, t = None, None
    res = {}

    name_field = one_signal.get('name_field', 'fields')
    name_species = one_signal.get('name_species', None)
    sel_e = one_signal.get('sel_e', 'none')

    if opt_var == 'n1':
        n_mode_chosen = one_signal['n1']
        chi_point = one_signal.get('chi-point', 0.0)

        ff, name_var = None, None
        if name_field == 'fields':
            init(dd, nn=[n_mode_chosen])
            ff = dd['3d']['fields']
            name_var = '\Phi'
        elif name_field == 'density':
            init_den(dd, name_species, nn=[n_mode_chosen])
            ff = dd['3d']['density-' + name_species]
            name_var = 'n({:s})'.format(name_species)
        else:
            mix.error_mes('Wrong structure name.')
        s = ff['s']
        t = ff['t']

        vvar = signal_n_allm_ts(dd['3d'], ff, n_mode_chosen, chi_point, sel_e=sel_e)

        if name_field == 'fields':
            if sel_e == 'none':
                name_var = '\Phi'
            if sel_e == 'er':
                vvar = - np.gradient(vvar, s, axis=1)
                name_var = 'E_r'
            if sel_e == 'echi':
                name_var = 'E_\chi' + '(\chi = {:0.3f})'.format(chi_point)

        tit_var = '3D:\ <' + name_var + '(n = {:d})'.format(n_mode_chosen) + '>_m'
    if opt_var == 'n1-antenna':
        n_mode_chosen = one_signal['n1']
        chi_point = one_signal.get('chi-point', 0.0)

        init(dd, nn=[n_mode_chosen], flag_antenna=True)
        ff = dd['3d']['fields']
        s = ff['s']
        t = ff['t']

        vvar = signal_n_allm_ts(dd['3d'], ff, n_mode_chosen, chi_point)
        tit_var = '3D:\ <\Phi_{ant}' + '(n = {:d})'.format(n_mode_chosen) + '>_m'
    if opt_var == 'n1-m1':
        n_mode_chosen = one_signal['n1']
        m_mode_chosen = one_signal['m1']
        chi_point = one_signal.get('chi-point', 0.0)

        ff, name_var = None, None
        if name_field == 'fields':
            init(dd, nn=[n_mode_chosen])
            ff = dd['3d']['fields']
            name_var = '\Phi'
        elif name_field == 'density':
            init_den(dd, name_species, nn=[n_mode_chosen])
            ff = dd['3d']['density-' + name_species]
            name_var = 'n({:s})'.format(name_species)
        else:
            mix.error_mes('Wrong structure name.')

        s = ff['s']
        t = ff['t']

        vvar = signal_n1_m1_ts(dd['3d'], ff, n_mode_chosen, m_mode_chosen,
                               sel_e=sel_e, chi_point=chi_point)

        if name_field == 'fields':
            if sel_e == 'none':
                name_var = '\Phi'
            if sel_e == 'er':
                vvar = - np.gradient(vvar, s, axis=1)
                name_var = 'E_r'
            if sel_e == 'echi':
                name_var = 'E_\chi' + '(\chi = {:0.3f})'.format(chi_point)

        tit_var = '3D:\ ' + name_var + \
            '(n = {:d},'.format(n_mode_chosen) + \
            'm = {:d})'.format(m_mode_chosen)
    if opt_var == 'n1-antenna-init':
        n_mode_chosen = one_signal['n1']
        chi_point = one_signal.get('chi-point', 0.0)
        file_name = one_signal.get('file_name', 'user_antenna.h5')
        flag_orb_shape = one_signal.get('flag_orb_shape', False)

        structure_name = 'antenna-init'
        init_antenna(dd, name_antenna_file=file_name, structure_name=structure_name,
                     flag_orb_shape=flag_orb_shape)
        ff = dd['3d'][structure_name]
        s = ff['s']
        t = ff['t']

        vvar = signal_n_allm_ts(dd['3d'], ff, n_mode_chosen, chi_point)
        tit_var = '3D:\ Init.\ Antenna\ :\ <\Phi' + '(n = {:d})'.format(n_mode_chosen) + '>_m'
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

        vvar = signal_chi1_phi1(dd['3d'], ff, chi_point, phi_point)
        line_chi = '\\chi = {:0.1f},\ '.format(chi_point)
        line_phi = '\\varphi = {:0.1f}'.format(phi_point)
        tit_var = '3D:\ Antenna:\ \\Phi' + '_{' + line_chi + line_phi + '}'
    if opt_var == 'potsc-antenna-init':
        chi_point = one_signal.get('chi-point', 0.0)
        phi_point = one_signal.get('phi-point', 0.0)
        file_name = one_signal.get('file_name', 'user_antenna.h5')
        flag_orb_shape = one_signal.get('flag_orb_shape', False)

        structure_name = 'antenna-init'
        init_antenna(dd, name_antenna_file=file_name, structure_name=structure_name,
                     flag_orb_shape=flag_orb_shape)
        ff = dd['3d'][structure_name]
        s = ff['s']
        t = ff['t']

        vvar = signal_chi1_phi1(dd['3d'], ff, chi_point, phi_point)
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

    name_field = one_signal.get('name_field', 'fields')
    name_species = one_signal.get('name_species', None)
    flag_er = one_signal.get('flag_er', False)

    if opt_var == 'n1':
        n_mode_chosen = one_signal['n1']
        t_point = one_signal['t-point']

        init(dd, nn=[n_mode_chosen])
        ff = dd['3d']['fields']
        s = ff['s']
        chi = ff['chi']

        vvar, t1 = signal_n_allm_t1(dd['3d'], ff, n_mode_chosen, t_point)
        tit_var = '3D:\ <\Phi' + \
                  '(n = {:d}, t = {:0.3e})'.format(n_mode_chosen, t1) \
                   + '>_m'
    if opt_var == 'n1-antenna':
        n_mode_chosen = one_signal['n1']
        t_point = one_signal['t-point']

        init(dd, nn=[n_mode_chosen], flag_antenna=True)
        ff = dd['3d']['fields']
        s = ff['s']
        chi = ff['chi']

        vvar, t1 = signal_n_allm_t1(dd['3d'], ff, n_mode_chosen, t_point)
        tit_var = '3D:\ <\Phi_{ant}' + \
                  '(n = {:d}, t = {:0.3e})'.format(n_mode_chosen, t1) \
                   + '>_m'
    if opt_var == 'n1-m1':
        n_mode_chosen = one_signal['n1']
        m_mode_chosen = one_signal['m1']
        t_point = one_signal['t-point']

        ff, name_var = None, None
        if name_field == 'fields':
            init(dd, nn=[n_mode_chosen])
            ff = dd['3d']['fields']
            name_var = '\Phi'
        elif name_field == 'density':
            init_den(dd, name_species, nn=[n_mode_chosen])
            ff = dd['3d']['density-' + name_species]
            name_var = 'n({:s})'.format(name_species)
        else:
            mix.error_mes('Wrong structure name.')

        s = ff['s']
        chi = ff['chi']

        vvar, t1 = signal_n1_m1_t1(dd['3d'], ff, n_mode_chosen, m_mode_chosen, t_point)

        if flag_er and name_field == 'fields':
            vvar = - np.gradient(vvar, s, axis=1)
            name_var = 'E_r'

        tit_var = '3D:\ ' + name_var + \
            '(n = {:d}, m = {:d}, t = {:0.3e})'.format(
                n_mode_chosen, m_mode_chosen, t_point
            )
    if opt_var == 'n1-antenna-init':
        n_mode_chosen = one_signal['n1']
        t_point = one_signal['t-point']
        file_name = one_signal.get('file_name', 'user_antenna.h5')
        flag_orb_shape = one_signal.get('flag_orb_shape', False)

        structure_name = 'antenna-init'

        init_antenna(dd, name_antenna_file=file_name, structure_name=structure_name,
                     flag_orb_shape=flag_orb_shape)
        ff = dd['3d'][structure_name]
        s = ff['s']
        chi = ff['chi']

        vvar, t1 = signal_n_allm_t1(dd['3d'], ff, n_mode_chosen, t_point)
        tit_var = '3D:\ Init.\ Antenna:\ <\Phi' + \
                  '(n = {:d}, t = {:0.3e})'.format(n_mode_chosen, t1) \
                  + '>_m'
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
        file_name = one_signal.get('file_name', 'user_antenna.h5')
        flag_orb_shape = one_signal.get('flag_orb_shape', False)

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



