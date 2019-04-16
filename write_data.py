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


def save_data(path_file, name, data):
    ff = h5.File(path_file, 'w')
    ff.create_dataset(name, data=data)
    ff.close()


