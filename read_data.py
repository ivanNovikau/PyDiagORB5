import Mix as mix
import numpy as np
import h5py as h5
import math


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)


def potsc(dd, path_to_folder):
    if 'potsc' in dd:
        return dd

    path_to_file = path_to_folder + '/orb5_res.h5'

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
    return dd