import Mix as mix
import ControlPlot as cpr
import ymath
import curve as crv
import numpy as np


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)


def read_profiles(path_to_read, file_names):
    # open initial files
    ff = []
    full_names_read = []
    for one_file_name in file_names:
        one_full_name = path_to_read + '/' + one_file_name
        full_names_read.append(one_full_name)
        ff.append(open(one_full_name, "r"))

    # read data from initial files
    res = {}
    count_file = -1
    for f in ff:
        count_file += 1
        ns = int(f.readline())
        psi, T, n, vp = \
            np.zeros(ns), np.zeros(ns), np.zeros(ns), np.zeros(ns)
        id_el = -1
        for one_line in f:
            id_el += 1
            psi[id_el], T[id_el], n[id_el], vp[id_el] = \
                [float(x) for x in one_line.split()]
        res[file_names[count_file]] = {
            'ns': ns,
            'psi': psi,
            'T': T,
            'n': n,
            'vp': vp
        }

    # close files
    for one_ff in ff:
        one_ff.close()

    # results
    return res


def write_profiles(path_to_write, file_names, profs):
    # open/create files to write data to
    ff = []
    full_names_read = []
    for one_file_name in file_names:
        one_full_name = path_to_write + '/' + one_file_name
        full_names_read.append(one_full_name)
        ff.append(open(one_full_name, "w"))

    # write data to the files
    count_file = -1
    for f in ff:
        count_file += 1
        prof = profs[file_names[count_file]]
        f.write('{:d}\n'.format(prof['ns']))
        for id_s in range(prof['ns']):
            f.write('{:.6e}   {:.6e}   {:.6e}   {:.6e}\n'.format(
                prof['psi'][id_s], prof['T'][id_s],
                prof['n'][id_s], prof['vp'][id_s]
            ))

    # close files
    for one_ff in ff:
        one_ff.close()


def plot_profiles(path_to_read, file_names):
    # read initial profiles
    profs = read_profiles(path_to_read, file_names)

    # reorganise data
    psi, T, n, vp = {}, {}, {}, {}
    for one_name in file_names:
        prof = profs[one_name]
        psi[one_name], T[one_name], n[one_name], vp[one_name] = \
            prof['psi'], prof['T'], prof['n'], prof['vp']

    # plot data:
    curves_T = crv.Curves().xlab('\psi/\psi_{edge}').ylab('T')
    curves_n = crv.Curves().xlab('\psi/\psi_{edge}').ylab('n')
    curves_vp = crv.Curves().xlab('\psi/\psi_{edge}').ylab('v_\parallel')
    id_curve = -1
    for one_name in file_names:
        id_curve += 1
        prof = profs[one_name]
        leg_name = '\_'.join(one_name.split('_'))
        curves_T.new(one_name).XS(prof['psi']).YS(prof['T'])\
            .leg(leg_name).new_sty(id_curve)
        curves_n.new(one_name).XS(prof['psi']).YS(prof['n'])\
            .leg(leg_name).new_sty(id_curve)
        curves_vp.new(one_name).XS(prof['psi']).YS(prof['vp'])\
            .leg(leg_name).new_sty(id_curve)
    cpr.plot_curves(curves_T)
    cpr.plot_curves(curves_n)
    cpr.plot_curves(curves_vp)


def build_vpp_profile(path_to_read, path_to_write, file_names):
    # read initial profiles
    init_profs = read_profiles(path_to_read, file_names)

    # modify profiles
    res_profs = dict(init_profs)
    for one_name in file_names:
        prof = init_profs[one_name]

        # *** Flat profile ***
        vp_value = -5
        prof['vp'] = flat_profile(prof['ns'], vp_value)

        # # *** Gaussian profile ***
        # A0  = -1.0
        # mu  =  0.8
        # sig =  0.05
        # prof['vp'] = gaussian_profile(prof['psi'], A0, mu, sig)

        # # *** Double Gaussian profile ***
        # A1, A2     = 1, -1
        # mu1, mu2   = 0.8, 0.6
        # sig1, sig2 = 0.05, 0.05
        # prof['vp'] = double_gaussian_profile(prof['psi'],
        #                               A1, mu1, sig1,
        #                               A2, mu2, sig2)

    # write result profiles
    write_profiles(path_to_write, file_names, res_profs)


def flat_profile(ns, value):
    res = np.array([value for i in range(ns)])
    return res


def gaussian_profile(psi, A0, mu, sig):
    def gaussian(x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    res = np.array([A0 * gaussian(one_psi, mu, sig) for one_psi in psi])
    return res


def double_gaussian_profile(psi, A1, mu1, sig1, A2, mu2, sig2):
    def gaussian(x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    res = np.array([
        A1*gaussian(one_psi, mu1, sig1) + A2*gaussian(one_psi, mu2, sig2)
        for one_psi in psi
    ])
    return res
