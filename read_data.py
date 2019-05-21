import Mix as mix
import Constants as cn
import Species as Species
import ymath
import write_data as wr
import numpy as np
import h5py as h5
import math
import platform


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cn)
    mix.reload_module(Species)
    mix.reload_module(ymath)
    mix.reload_module(wr)


def read_signal(path_to_read, var_path):
    # so far, works only for dataset and one-stage group

    f = h5.File(path_to_read, 'r')
    if var_path in f:
        vvar = f[var_path]
        if isinstance(vvar, h5.Dataset):
            vvar = np.array(vvar)
        if isinstance(vvar, h5.Group):
            vvar = dict(vvar)
            for one_key in vvar.keys():
                vvar[one_key] = np.array(vvar[one_key])
    else:
        return None
    f.close()
    return vvar


def potsc(dd):
    if 'potsc' in dd:
        return
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')

    t = np.array(f['/data/var2d/generic/potsc/time'])
    s = np.array(f['/data/var2d/generic/potsc/coord1'])
    chi = np.array(f['/data/var2d/generic/potsc/coord2'])
    r = np.array(f['/data/var2d/generic/potsc/rsc'])
    z = np.array(f['/data/var2d/generic/potsc/zsc'])
    potsc_data = np.array(f['/data/var2d/generic/potsc/data'])
    dd['potsc'] = {
        't': t,
        's': s,
        'chi': chi,
        'r': r,
        'z': z,
        'data': potsc_data}


def potsc_grids(dd):
    var_name = 'potsc_grids'
    if var_name in dd:
        return

    data = read_signal(dd['path_ext'], var_name)
    if data is None:
        data = {}

        # read data from orb5 output file
        f = h5.File(dd['path_orb'], 'r')
        data['t']   = np.array(f['/data/var2d/generic/potsc/time'])
        data['chi'] = np.array(f['/data/var2d/generic/potsc/coord2'])
        data['s']   = np.array(f['/data/var2d/generic/potsc/coord1'])
        data['r']   = np.array(f['/data/var2d/generic/potsc/rsc'])
        data['z']   = np.array(f['/data/var2d/generic/potsc/zsc'])
        f.close()

        # save data to an external file
        desc = 'potsc grids'
        wr.save_data_adv(dd['path_ext'], data, {'name': var_name, 'desc': desc})

    # save data to the structure
    dd[var_name] = data


# full potential at some angle chi
def potsc_chi(dd, oo):
    potsc_grids(dd)
    chi_s = oo.get('chi_s', [0.0])
    nchi, names = np.size(chi_s), []
    for count_chi in range(nchi):
        one_chi = chi_s[count_chi]
        var_name = 'potsc-chi-' + '{:0.3f}'.format(one_chi)
        names.append(var_name)
        if var_name in dd:
            continue

        data = read_signal(dd['path_ext'], var_name)
        if data is None:
            data = {}

            # read data from orb5 output file
            f = h5.File(dd['path_orb'], 'r')
            data['id_chi_1'], data['chi_1'] = \
                mix.find(dd['potsc_grids']['chi'], one_chi)
            data['data'] = np.array(
                f['/data/var2d/generic/potsc/data'][:, data['id_chi_1'], :])
            f.close()

            # save data to an external file
            desc = 'full potential at phi=0 and t = {:0.3f}'.format(data['chi_1'])
            wr.save_data_adv(dd['path_ext'], data, {'name': var_name, 'desc': desc})

        # save data to the structure
        dd[var_name] = data
    return names


# potential at some time moment
def potsc_t(dd, oo):
    potsc_grids(dd)
    t_points = oo.get('t_points', [0.0])
    nt_points = np.size(t_points)
    names = []
    for count_t in range(nt_points):
        one_t = t_points[count_t]
        var_name = 'potsc-t-' + '{:0.3e}'.format(one_t)
        names.append(var_name)
        if var_name in dd:
            continue

        data = read_signal(dd['path_ext'], var_name)
        if data is None:
            data = {}

            # read data from orb5 output file
            f = h5.File(dd['path_orb'], 'r')
            data['id_t_point'], data['t_point'] = \
                mix.find(dd['potsc_grids']['t'], one_t)
            data['data'] = \
                np.array(f['/data/var2d/generic/potsc/data'][data['id_t_point'], :, :])
            f.close()

            # save data to an external file
            desc = 'full potential at phi=0 and t = {:0.3e}'.format(data['t_point'])
            wr.save_data_adv(dd['path_ext'], data, {'name': var_name, 'desc': desc})

        # save data to the structure
        dd[var_name] = data
    return names


# potential at some radial points (data are saved):
def potsc_s(dd, oo):
    potsc_grids(dd)
    ss = oo.get('ss', [0.0])
    ns = np.size(ss)
    names = []
    for count_s in range(ns):
        one_s = ss[count_s]
        var_name = 'potsc-s-' + '{:0.3f}'.format(one_s)
        names.append(var_name)
        if var_name in dd:
            continue

        data = read_signal(dd['path_ext'], var_name)
        if data is None:
            data = {}

            # read data from orb5 output file
            f = h5.File(dd['path_orb'], 'r')
            data['id_s_1'], dd['s_1'] = \
                mix.find(dd['potsc_grids']['s'], one_s)
            data['data'] = np.array(
                f['/data/var2d/generic/potsc/data'][:, :, data['id_s_1']])
            f.close()

            # save data to an external file
            desc = 'full potential at phi=0 and s = {:0.3f}'.format(dd['s_1'])
            wr.save_data_adv(dd['path_ext'], data, {'name': var_name, 'desc': desc})

        # save data to the structure
        dd[var_name] = data
    return names


def phibar(dd):
    if 'phibar' in dd:
        return
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')

    t = np.array(f['/data/var1d/generic/phibar/time'])
    s = np.array(f['/data/var1d/generic/phibar/coord1'])
    phibar_data = np.array(f['/data/var1d/generic/phibar/data'])
    dd['phibar'] = {
        't': t,
        's': s,
        'data': phibar_data}


def radial_heat_flux(dd):
    if 'efluxw_rad' in dd:
        return
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')
    for sp_name in dd['kin_species_names']:
        dd[sp_name].radial_heat_flux(f)


def q(dd):
    if 'q' in dd:
        return
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')
    s = np.array(f['/equil/profiles/generic/sgrid_eq'])
    data = np.array(f['/equil/profiles/generic/q'])

    dd['q'] = {
        's': s,
        'data': data
    }


def vti(dd):
    if 'vti' in dd:
        return
    mass_pf = dd['pf'].mass
    Ti_peak = dd['pf'].T_speak(dd)
    dd['vti'] = ymath.find_vt(Ti_peak, mass_pf)


def elongation(dd):
    if 'elong' in dd:
        return
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')
    dd['elong'] = f['/equil/scalars/generic/e_mid'][0]


def nT_evol(dd, species_name):
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')
    if species_name is 'pf':
        dd['pf'].find_nT_evol(dd, f)
    else:
        dd['kin_species'][species_name].find_nT_evol(dd, f)


def krpert(dd, species_name):
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')
    if species_name is 'pf':
        dd['pf'].find_krpert(dd, f)
    else:
        dd['kin_species'][species_name].find_krpert(dd, f)


def init(dd):
    path_to_file = dd['path'] + '/orb5_res.h5'
    dd['path_orb'] = path_to_file
    f = h5.File(path_to_file, 'r')

    # max amount of memory to occupy for one array (in gigabytes):
    dd['max_size_Gb'] = 1.5

    # define an operational system
    sys_current = platform.system()
    dd['oper_system'] = {
        'Windows'.lower(): cn.SYS_WINDOWS,
        'Linux'.lower(): cn.SYS_WINDOWS
    }.get(sys_current.lower(), cn.SYS_UNDEFINED)

    # define path to a file to save data from this project:
    dd['path_ext'] = dd['path'] + '/saved_data.h5'
    wr.open_file(dd['path_ext'], False)

    # number of starts:
    input_files = f['/run_info/input']
    dd['n_starts'] = len(list(input_files))

    # initialization of the species:
    species(dd, f)

    # read basic parameters
    dd['Lx'] = f['/parameters/equil/lx'][0]
    dd['sfmin'] = f['/parameters/fields/sfmin'][0]
    dd['sfmax'] = f['/parameters/fields/sfmax'][0]

    # calculate basic variables:
    mass_pf = dd['pf'].mass
    Z_pf = dd['pf'].Z
    B0 = dd['B0']
    Te_peak = dd['electrons'].T_speak(dd)

    dd['wc'] = ymath.find_wc(B0, mass_pf, Z_pf)
    dd['cs'] = ymath.find_cs(Te_peak, mass_pf)

    dd['Lwork'] = (dd['sfmax'] - dd['sfmin']) * dd['a0']

    # read profiles:
    for sp_name in dd['species_names']:
        dd[sp_name].nT(dd, f)

    # T,n dependent variables:
    dd['T_speak'] = dd['pf'].T_speak(dd)
    dd['rhoL_speak'] = ymath.find_rhoL(dd['T_speak'], dd['B0'], dd['pf'].mass, dd['pf'].Z)


# init ->
def species(dd, f):
    if 'kin_species' in dd:
        return

    # read species' names
    names_all_species = f['/parameters/species_names'].attrs
    ids_attr = list(names_all_species)

    dd['species_names'] = [names_all_species[name].decode("utf-8")
                           for name in ids_attr[1:len(ids_attr)]]

    # correct species' names
    if dd['oper_system'] == cn.SYS_WINDOWS:
        correct_species_names_win(dd, f)

    # create species' objects
    count_species = 0
    dd['kin_species_names'] = []
    dd['kin_species'] = {}
    for namesp in dd['species_names']:
        count_species += 1
        one_species = Species.Species(namesp, dd, f)
        if one_species.is_kinetic:
            dd['kin_species_names'].append(namesp)
            dd['kin_species'][namesp] = one_species
        if count_species == 1:
            dd['pf'] = one_species
        dd[namesp] = one_species


# init -> species ->
def correct_species_names_win(dd, f):
    keys = list(f['/parameters'].keys())
    count_sp = -1
    for spname in dd['species_names']:
        count_sp += 1
        for key1 in keys:
            if spname.lower() in key1.lower():
                dd['species_names'][count_sp] = key1
                break













