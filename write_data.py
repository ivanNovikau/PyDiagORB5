import Mix as mix
import Constants as cn
import Species as Species
import ymath
import numpy as np
import h5py as h5
import math


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cn)
    mix.reload_module(Species)
    mix.reload_module(ymath)


def create_file(path_file):
    ff = h5.File(path_file, 'w')
    ff.close()


def save_data(path_file, name, data, desc=''):
    ff = h5.File(path_file, 'a')
    if name not in ff:
        ddata = ff.create_dataset(name, data=data)
        ddata.attrs[u'descr'] = desc
    ff.close()


