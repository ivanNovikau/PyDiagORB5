import Mix as mix
import Constants as cn
import Species
import numpy as np
import h5py as h5
import math
import platform


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cn)
    mix.reload_module(Species)


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


# init ->
def species(dd, f):
    if 'species' in dd:
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
    for namesp in dd['oper_system']:
        count_species += 1
        one_species = Species.Species(namesp, dd, f)
        if one_species.is_kinetic:
            dd['kin_species_names'].append(namesp)
            dd['kin_species'][namesp] = one_species
        if count_species == 1:
            dd['pf_species'] = one_species
        dd['species'][namesp] = one_species


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









