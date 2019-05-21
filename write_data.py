import Mix as mix
import Constants as cn
import Species as Species
import ymath
import numpy as np
import h5py as h5
import math
import os.path


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


def open_file(path_file, flag_new_file):
    if flag_new_file:
        create_file(path_file)
    else:
        if os.path.isfile(path_file):
            ff = h5.File(path_file, 'a')
            ff.close()
        else:
            create_file(path_file)


def save_data_adv(path_file, data, oo):
    TYPE_ARRAY = 'array'
    TYPE_DICT = 'dict'

    name = oo.get('name', [])
    desc = oo.get('desc', None)

    sel_var = TYPE_ARRAY
    if isinstance(data, dict):
        sel_var = 'dict'

    ff = h5.File(path_file, 'a')
    if sel_var == TYPE_ARRAY:
        if name in ff:
            del ff[name]
        ddata = ff.create_dataset(name, data=data)
        ddata.attrs[u'descr'] = desc
    if sel_var == TYPE_DICT:
        if name in ff:
            del ff[name]
        grp = ff.create_group(name)
        grp.attrs[u'descr'] = desc
        for one_key in data.keys():
            ddata = grp.create_dataset(one_key, data=data[one_key])
            ddata.attrs[u'descr'] = desc

    ff.close()


