import Mix as mix
import ymath
import numpy as np
from scipy import constants


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(ymath)


class Species:
    name = ''
    is_kinetic = None
    is_electrons = None
    mass_rel = np.nan
    mass = np.nan
    Z = np.nan

    tau = np.nan

    nsel_profile = None

    nT_equil = None
    nT_evol = None

    efluxw_rad = None

    krpert = None

    # in general, rd.init -> rd.species ->
    def __init__(self, name, dd, f):
        self.name = name
        self.is_kinetic = f['/parameters/' + name + '/is_kinetic'][0]
        self.is_electrons = f['/parameters/' + name + '/is_electrons'][0]
        self.mass_rel = f['/parameters/' + name + '/mass'][0]
        self.Z = f['/parameters/' + name + '/charge'][0]
        self.mass = dd['mass_pf'] * self.mass_rel

        if self.is_electrons:
            self.tau = 1
        else:
            self.tau = f['/parameters/' + name + '/tau'][0]

    def nT(self, dd, f):
        if self.nT_equil is not None:
            return

        self.nsel_profile = mix.get_attribute(f, '/parameters/' + self.name + '/nsel_profile')
        self.nsel_profile = self.nsel_profile[1]

        if self.nsel_profile == '2':
            s = np.array(f['/equil/profiles/generic/sgrid_eq'])
            psi = s**2
        else:
            s = np.array(f['/equil/profiles/' + self.name + '/s_prof'])
            psi = np.array(f['/equil/profiles/' + self.name + '/psi_prof'])
        n = np.array(f['/equil/profiles/' + self.name + '/n_pic'])
        T = np.array(f['/equil/profiles/' + self.name + '/t_pic'])
        vp = np.array(f['/equil/profiles/' + self.name + '/v_pic'])

        T_J   = T * dd['electrons'].T_speak(dd) * self.tau
        T_keV = T_J / (1e3 * constants.elementary_charge)

        gradT    = np.gradient(T, s)
        grad_logT = np.gradient(np.log(T), s)
        gradn    = np.gradient(n, s)
        grad_logn = np.gradient(np.log(n), s)

        self.nT_equil = {
            's': s,
            'psi': psi,
            'T': T, 'T_J': T_J, 'T_keV': T_keV,
            'gradT': gradT, 'grad_logT': grad_logT,
            'n': n,
            'gradn': gradn, 'grad_logn': grad_logn,
            'vp': vp
        }

    def find_nT_evol(self, dd, f):
        if self.nT_evol is not None:
            return
        t = np.array(f['/data/var1d/' + self.name + '/f_av/time'])
        s  = np.array(f['/data/var1d/' + self.name + '/f_av/coord1'])
        n    = np.array(f['/data/var1d/' + self.name + '/f_av/data'])
        vperp2f = np.array(f['/data/var1d/' + self.name + '/v_perp2_av/data'])
        u2f = np.array(f['/data/var1d/' + self.name + '/v_par2_av/data'])
        uf = np.array(f['/data/var1d/' + self.name + '/v_par_av/data'])

        # Remark: at the time of creating this script,
        # n, vperp2f, u2f, uf have the same s-grid

        p = 1./3 * self.mass_rel * (vperp2f + u2f - uf**2/n)
        T = self.tau * p / n

        grad_T    = np.gradient(T, s, axis=1)
        grad_logT = np.gradient(np.log(T), s, axis=1)
        grad_n    = np.gradient(n, s, axis=1)
        grad_logn = np.gradient(np.log(n), s, axis=1)

        # results:
        self.nT_evol = {
            's': s, 't': t,
            'n': n, 'gradn': grad_n, 'grad_logn': grad_logn,
            'T': T, 'gradT': grad_T, 'grad_logT': grad_logT
        }

    def find_krpert(self, dd, f):
        if self.krpert is not None:
            return
        self.krpert = f['/parameters/' + self.name + '/kr_pert'][0]

    def T_speak(self, dd):
        # get T at speak (in J) for the species
        mass_pf = dd['pf'].mass
        Z_pf = dd['pf'].Z
        B0 = dd['B0']
        a0 = dd['a0']
        Lx = dd['Lx']
        T_peak = ymath.find_Ti(self.tau, Lx, a0, B0, mass_pf, Z_pf)
        return T_peak

    def radial_heat_flux(self, f):
        if self.efluxw_rad is not None:
            return
        s = np.array(f['/data/var1d/' + self.name + '/efluxw_rad/coord1'])
        t = np.array(f['/data/var1d/' + self.name + '/efluxw_rad/time'])
        data = np.array(f['/data/var1d/' + self.name + '/efluxw_rad/data'])
        self.efluxw_rad = {
            's': s,
            't': t,
            'data': data
        }



