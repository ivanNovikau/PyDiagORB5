import numpy as np
from scipy import constants


def reload():
    # Important: put here all modules that you want to reload
    dummy = 0


class Species:
    name = ''
    is_kinetic = None
    is_electrons = None
    mass_rel = np.nan
    mass = np.nan
    Z = np.nan

    # in general, rd.init -> rd.species ->
    def __init__(self, name, dd, f):
        self.name = name
        self.is_kinetic = f['/parameters/' + name + '/is_kinetic'][0]
        self.is_electrons = f['/parameters/' + name + '/is_electrons'][0]
        self.mass_rel = f['/parameters/' + name + '/mass'][0]
        self.Z = dd['/parameters/' + name + '/charge'][0]
        self.mass = dd['mass_pf'] * self.mass_rel
