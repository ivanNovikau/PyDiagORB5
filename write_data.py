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


def create_open_file(dd):
    path_to_write = dd['path_to_write']
    ff = h5.File(path_to_write, 'w')
    return ff


def save_data(ff, name, data):
    ff.create_dataset(name, data=data)


def close_file(ff):
    ff.close()
