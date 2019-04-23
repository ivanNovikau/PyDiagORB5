import Mix as mix
import Constants as cn
import Species as Species
import ymath
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
    if 'potsc_grids' in dd:
        return
    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')
    t = np.array(f['/data/var2d/generic/potsc/time'])
    chi = np.array(f['/data/var2d/generic/potsc/coord2'])
    s = np.array(f['/data/var2d/generic/potsc/coord1'])
    r = np.array(f['/data/var2d/generic/potsc/rsc'])
    z = np.array(f['/data/var2d/generic/potsc/zsc'])

    dd['potsc_grids'] = {
        't': t, 'chi': chi, 's': s, 'r': r, 'z': z
    }
    f.close()


def potsc_chi(dd, oo):
    potsc_grids(dd)
    chi_s = oo.get('chi_s', [0.0])

    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')

    nchi = np.size(chi_s)
    names_potsc_chi = []
    for count_chi in range(nchi):
        one_chi = chi_s[count_chi]
        name_potsc_chi = 'potsc-chi-' + '{:0.3f}'.format(one_chi)
        names_potsc_chi.append(name_potsc_chi)

        if name_potsc_chi in dd:
            continue

        id_chi, chi1 = mix.find(dd['potsc_grids']['chi'], one_chi)

        potsc_data = np.array(f['/data/var2d/generic/potsc/data'][:, id_chi, :])
        dd[name_potsc_chi] = {
            'chi_1': chi1, 'id_chi_1': id_chi,
            'data': potsc_data}

    f.close()

    return names_potsc_chi


def potsc_t(dd, oo):
    potsc_grids(dd)
    t_points = oo.get('t_points', [0.0])

    path_to_file = dd['path'] + '/orb5_res.h5'
    f = h5.File(path_to_file, 'r')

    nt_points = np.size(t_points)
    names_potsc_t = []
    for count_t in range(nt_points):
        one_t = t_points[count_t]
        name_potsc_t = 'potsc-t-' + '{:0.3f}'.format(one_t)
        names_potsc_t.append(name_potsc_t)

        if name_potsc_t in dd:
            continue

        id_t, t_point = mix.find(dd['potsc_grids']['t'], one_t)

        potsc_data = np.array(f['/data/var2d/generic/potsc/data'][id_t, :, :])
        dd[name_potsc_t] = {
            't_point': t_point, 'id_t_point': id_t,
            'data': potsc_data}

    f.close()

    return names_potsc_t


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
    f = h5.File(path_to_file, 'r')

    # define an operational system
    sys_current = platform.system()
    dd['oper_system'] = {
        'Windows'.lower(): cn.SYS_WINDOWS,
        'Linux'.lower(): cn.SYS_WINDOWS
    }.get(sys_current.lower(), cn.SYS_UNDEFINED)

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













