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


def save_profiles_from_project(dd, path_to_write):
    file_names = []
    profs = {}
    for sp_name in dd['species_names']:
        file_names.append(sp_name + '_profiles.dat')
        psi = dd[sp_name].nT_equil['s']**2
        ns = np.size(psi)
        T = dd[sp_name].nT_equil['T']
        n = dd[sp_name].nT_equil['n']
        vp = dd[sp_name].nT_equil['vp']
        profs[file_names[-1]] = {
            'ns': ns,
            'psi': psi,
            'T': T,
            'n': n,
            'vp': vp
        }
    write_profiles(path_to_write, file_names, profs)


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

    # curves_T = crv.Curves().xlab('s').ylab('T')
    curves_n = crv.Curves().xlab('s').ylab('n')
    # curves_vp = crv.Curves().xlab('s').ylab('v_\parallel')
    id_curve = -1
    for one_name in file_names:
        id_curve += 1
        prof = profs[one_name]
        leg_name = '\_'.join(one_name.split('_'))
        # curves_T.new(one_name).XS(prof['psi']).YS(prof['T']) \
        #     .leg(leg_name).new_sty(id_curve)
        curves_n.new(one_name).XS(np.sqrt(prof['psi'])).YS(prof['n']) \
            .leg(leg_name).new_sty(id_curve)
        # curves_vp.new(one_name).XS(prof['psi']).YS(prof['vp']) \
        #     .leg(leg_name).new_sty(id_curve)
    # cpr.plot_curves(curves_T)
    cpr.plot_curves(curves_n)
    # cpr.plot_curves(curves_vp)


def compare_profiles_project_files(dd, path_to_files, oo):
    # function searches for species_profiles.dat files
    # (such as deuterium_profiles.dat etc) and compare profiles
    # from these files to the profiles from the project dd:

    species_to_compare = oo.get('species_to_compare', [])

    # read profiles from files
    file_names = []
    for sp_name in dd['species_names']:
        if sp_name in species_to_compare:
            file_names.append(sp_name + '_profiles.dat')
    profs_files = read_profiles(path_to_files, file_names)

    project_name = dd.get('project_name', 'project1')
    project_file_name = oo.get('project_file_name', 'project2')

    # plot profiles:
    count_file = -1
    for sp_name in species_to_compare:
        count_file = count_file + 1
        file_name = file_names[count_file]
        proff = profs_files[file_name]
        leg_name_proj = project_name + ':\ ' + sp_name
        leg_name_file = project_file_name + ':\ ' + sp_name

        curves_T  = crv.Curves().xlab('s').ylab('T')
        curves_n  = crv.Curves().xlab('s').ylab('n')
        curves_vp = crv.Curves().xlab('s').ylab('v_\parallel')

        curves_T.new() \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['T']) \
            .leg(leg_name_proj).new_sty(0)
        curves_T.new()\
            .XS(np.sqrt(proff['psi']))\
            .YS(proff['T']) \
            .leg(leg_name_file).new_sty(1)

        curves_n.new() \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['n']) \
            .leg(leg_name_proj).new_sty(0)
        curves_n.new() \
            .XS(np.sqrt(proff['psi'])) \
            .YS(proff['n']) \
            .leg(leg_name_file).new_sty(1)

        curves_vp.new() \
            .XS(dd[sp_name].nT_equil['s']) \
            .YS(dd[sp_name].nT_equil['vp']) \
            .leg(leg_name_proj).new_sty(0)
        curves_vp.new() \
            .XS(np.sqrt(proff['psi'])) \
            .YS(proff['vp']) \
            .leg(leg_name_file).new_sty(1)

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


def build_new_profiles(path_to_read, path_to_write, file_names):
    # read initial profiles
    sp_profs = read_profiles(path_to_read, file_names)

    # modify profiles
    for one_name in file_names:
        prof = sp_profs[one_name]

        # region *** Flat profiles ***
        n_value, T_value, vp_value = 1, 1, 1
        # prof['T']  = flat_profile(prof['ns'], T_value)
        # prof['n']  = flat_profile(prof['ns'], n_value)
        prof['vp'] = flat_profile(prof['ns'], vp_value)

        # endregion

        # region *** Gaussian n-profiles (edge) for fast species ***

        # # * edge *
        # A0, mu, sig = 1.0, 0.8, 0.05
        # prof['n'] = gaussian_profile(prof['psi'], A0, mu, sig)

        # # * center *
        # A0, mu, sig = 1.0, 0.3, 0.09
        # prof['n'] = gaussian_profile(prof['psi'], A0, mu, sig)

        # # * s = 0.8 *
        # A0, mu, sig, Ab = 1.0, 0.64, 0.05, 1e-3
        # prof['n'] = gaussian_profile_WITH_background(
        #     prof['psi'], A0, mu, sig, Ab)

        # # * s = 0.8, ds = 2 *
        # A0, mu, sig, Ab = 1.0, 0.64, 0.1, 1e-3
        # prof['n'] = gaussian_profile_WITH_background(
        #     prof['psi'], A0, mu, sig, Ab)

        # # * s = 0.75, ds = 0.25 *
        # A0, mu, sig, Ab = 1.0, 0.5625, 0.15, 1e-3
        # prof['n'] = gaussian_profile_WITH_background(
        #     prof['psi'], A0, mu, sig, Ab)

        # # * s = 0.7, ds = 0.3 *
        # A0, mu, sig, Ab = 1.0, 0.49, 0.19, 1e-3
        # prof['n'] = gaussian_profile_WITH_background(
        #     prof['psi'], A0, mu, sig, Ab)

        # # * s = 0.4, ds = 0.6 *
        # # A0, mu, sigL, sigR, Ab = 1.0, 0.36, 0.12, 0.2, 1e-3
        # A0, mu, sigL, sigR, Ab = 1.0, 0.16, 0.07, 0.3, 1e-3
        # prof['n'] = gaussian_profile_asym(prof['psi'], A0, mu, sigL, sigR, Ab)

        # # * s = 0.35 *
        # A0, mu, sig, Ab = 1.0, 0.1225, 0.04, 1e-3
        # prof['ns'] = 201
        # prof['psi'] = rescale_psi(prof['ns'])
        # prof['n'] = gaussian_profile_WITH_background(
        #     prof['psi'], A0, mu, sig, Ab)

        # # * s = 0.45 *
        # A0, mu, sig, Ab = 1.0, 0.2025, 0.05, 1e-3
        # prof['ns'] = 121
        # prof['psi'] = rescale_psi(prof['ns'])
        # prof['n'] = gaussian_profile_WITH_background(
        #     prof['psi'], A0, mu, sig, Ab)
        #
        # # and flat T, vp profiles
        # T_value, vp_value = 1, 0
        # prof['T'] = flat_profile(prof['ns'], T_value)
        # prof['vp'] = flat_profile(prof['ns'], vp_value)

        # endregion

        # region *** Modify T at the edge ***
        # s_edge = 0.98
        # coef_mod = 50
        #
        # prof['ns'] = 256
        # psi_init = prof['psi']
        # prof['psi'] = rescale_psi(prof['ns'])
        # prof['T']  = np.interp(prof['psi'], psi_init, prof['T'])
        # prof['n']  = np.interp(prof['psi'], psi_init, prof['n'])
        # prof['vp'] = np.interp(prof['psi'], psi_init, prof['vp'])
        #
        # psi_edge = s_edge ** 2
        # for id_psi1 in range(prof['ns']):
        #     if prof['psi'][id_psi1] >= psi_edge:
        #         # prof['T'][id_psi1] = prof['T'][id_psi1]\
        #         #     * np.exp((psi_edge - prof['psi'][id_psi1]) * coef_mod)
        #         prof['T'][id_psi1] = prof['T'][id_psi1] + \
        #             (psi_edge - prof['psi'][id_psi1])**2 * coef_mod
        # endregion

    # write result profiles
    write_profiles(path_to_write, file_names, sp_profs)


def copy_T_from_project_to_file(dd, path_to_file, file_name,
                                path_to_res, species_name):
    # species_name - name of a species in project dd, which T will be copied

    # read initial profile
    sp_profs = read_profiles(path_to_file, [file_name])

    # copy profiles
    profs = sp_profs[file_name]
    proj_profs = dd[species_name]

    psi = proj_profs.nT_equil['s']**2
    T   = proj_profs.nT_equil['T']
    profs['T'] = np.interp(profs['psi'], psi, T)

    # write result profile
    write_profiles(path_to_res, [file_name], sp_profs)


def rescale_psi(npsi_new):
    psi_new = np.linspace(0, 1, npsi_new)
    return psi_new


def flat_profile(ns, value):
    res = np.array([value for i in range(ns)])
    return res


def gaussian_profile(psi, A0, mu, sig):
    def gaussian(x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    res = np.array([A0 * gaussian(one_psi, mu, sig) for one_psi in psi])
    return res


def gaussian_profile_WITH_background(psi, A0, mu, sig, Ab):
    def gaussian(x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    res = np.array([Ab + A0 * gaussian(one_psi, mu, sig)
                    for one_psi in psi])
    return res


def gaussian_profile_asym(psi, A0, mu, sig1, sig2, Ab):
    def gaussian(x, mu, sig1, sig2):
        if x <= mu:
            return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig1, 2.)))
        if x > mu:
            return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig2, 2.)))
    res = np.array([Ab + A0 * gaussian(one_psi, mu, sig1, sig2)
                    for one_psi in psi])
    return res


def double_gaussian_profile(psi, A1, mu1, sig1, A2, mu2, sig2):
    def gaussian(x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    res = np.array([
        A1*gaussian(one_psi, mu1, sig1) + A2*gaussian(one_psi, mu2, sig2)
        for one_psi in psi
    ])
    return res

