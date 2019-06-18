import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import zf_gam as zf
import general as gn
import write_data as wr
import transport
import Geom
import numpy as np
import scipy.io
from scipy.stats import norm as stat_norm
import matplotlib.mlab as mlab
import scipy.signal


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(zf)
    mix.reload_module(gn)
    mix.reload_module(wr)
    mix.reload_module(transport)
    mix.reload_module(Geom)


# NEW: take signal (vpar,mu)
def choose_one_var_vparmu(ovar, dd):
    opt_var      = ovar[0]
    species_name = ovar[1]

    vvar, vpar, mu, tit_var = None, None, None, ''
    if opt_var == 'jdote_es-mean-t':
        t_int = ovar[2]

        jdote_es_dict = dd[species_name].jdote_es(dd)
        t, vpar, mu = jdote_es_dict['t'], jdote_es_dict['vpar'], jdote_es_dict['mu']
        ids_t, t_int, line_t = mix.get_ids(t, t_int)

        vvar = np.mean(jdote_es_dict['data'][ids_t, :, :], axis=0)

        tit_var  = species_name + ':\ <J*E>_{t}'
        tit_var = 't = ' + line_t + ':\ ' + tit_var

    res = {
        'data': vvar,  # (t, mu, vpar)
        'vpar': vpar,
        'mu': mu,
        'tit': tit_var
    }

    return res


# NEW: variable (t)
def choose_one_var_t(ovar, dd):
    # --- read energy transfer signal ---
    def jdote_es_species(dd, species_name, mu_int, vpar_int):
        jdote_es_dict = dd[species_name].jdote_es(dd)
        t, vpar, mu = jdote_es_dict['t'], jdote_es_dict['vpar'], jdote_es_dict['mu']

        mu_int_work, vpar_int_work = mu_int, vpar_int
        if len(mu_int_work) == 0:
            mu_int_work = mu
        if len(vpar_int_work) == 0:
            vpar_int_work = vpar

        ids_mu,   mu_int_work,   line_mu   = mix.get_ids(mu,   mu_int_work, '{:0.2f}')
        ids_vpar, vpar_int_work, line_vpar = mix.get_ids(vpar, vpar_int_work, '{:0.1f}')

        vvar_loc = \
            jdote_es_dict['data'][:, ids_mu[0]:ids_mu[-1]+1, ids_vpar[0]:ids_vpar[-1]+1]
        vvar_loc = np.sum(vvar_loc, axis=2)
        vvar_loc = np.sum(vvar_loc[:, :], axis=1)

        return vvar_loc, t, line_mu, line_vpar

    # --- read field energy ---
    def efield_species(dd, species_name):
        efield_dict = dd[species_name].efield(dd)
        t        = efield_dict['t']
        vvar_loc = efield_dict['data']

        return vvar_loc, t

    # --- name of a signal ---
    opt_var = ovar[0]
    species_name = ovar[1]

    # --- velocity domains ---
    mu_int, vpar_int = [], []
    if len(ovar) > 2:
        mu_int = ovar[2]
        vpar_int = ovar[3]

    # --- Create empty result dictionary ---
    res = {}

    # --- choose signal ---
    vvar, t, tit_var = [], [], []
    if opt_var == 'jdote_es':
        if species_name is not 'total':
            vvar, t, line_mu, line_vpar = \
                jdote_es_species(dd, species_name, mu_int, vpar_int)
        else:
            sp_names = dd['kin_species_names']
            vvar, t, line_mu, line_vpar = \
                jdote_es_species(dd, sp_names[0], mu_int, vpar_int)
            for sp_name in sp_names[1:]:
                vvar_sp, _, _, _ = jdote_es_species(dd, sp_name, mu_int, vpar_int)
                vvar += vvar_sp
        tit_var = species_name + ':\ <J*E>_{\mu, v_{\parallel}}'
        res['line_sum'] = '\mu = ' + line_mu + '$\n$ v_{\parallel} = ' + line_vpar

    if opt_var == 'efield':
        if species_name is not 'total':
            vvar, t = efield_species(dd, species_name)
        else:
            sp_names = dd['kin_species_names']
            vvar, t = efield_species(dd, sp_names[0])
            for sp_name in sp_names[1:]:
                vvar_sp, _ = efield_species(dd, sp_name)
                vvar += vvar_sp
        tit_var = species_name + ':\ \mathcal{E}}'

    # --- save signal ---
    res.update({
        'data': vvar,
        'x': t,
        'tit': tit_var
    })

    return res


def data_GENE_CPS2019():
    path_data = 'd:/Work-Projects/MyProgs/ORB_data/MPR/nled/kin/linear/GENE-CPS2019-DF/data.mat'
    dGENE = scipy.io.loadmat(path_data)

    Te_div_B0 = dGENE['Te_B0'][0, 0]
    cs = dGENE['cs'][0, 0]
    Lref = dGENE['lref'][0, 0]
    vpar = dGENE['out_Ivan_v'][0]
    mu   = dGENE['out_Ivan_mu'][0]
    data_t_jet_jed_jef = dGENE['out_Ivan']

    t    = data_t_jet_jed_jef[:, 0]
    Ptot = data_t_jet_jed_jef[:, 1]  # total J*E
    Pth  = data_t_jet_jed_jef[:, 2]  # J*E of thermal species
    Pf   = data_t_jet_jed_jef[:, 3]  # J*E of fast species

    gamma_vmu = dGENE['out_Ivan_gamma_vm']  # (vpar, mu)

    res = {
        'Te_div_B0': Te_div_B0, 'cs': cs, 'Lref': Lref,
        'vpar': vpar, 'mu': mu, 't': t,
        'Ptot': Ptot, 'Pth': Pth, 'Pf': Pf,
        'gamma_vmu': gamma_vmu
    }

    return res


# Comparison with GENE for the paper CPS-2019:
def comparison_with_GENE_CPS2019(dd, oo):
    # Parameters
    dd_f0   = oo.get('dd_f0', None)

    labx_ms = 'v_{\parallel}(m/s)'
    labx_norm = 'v_{\parallel}'

    leg_gene = 'GENE:\ ' + '-\langle \mathcal{P} \\rangle_{\mu}'
    leg_orb  = 'ORB5:\ ' + '\langle \mathcal{P} \\rangle_{\mu}'

    # Find iniitial distribution functions of fast species
    rd.distribution_1d(dd_f0, 'fast')

    # velocity normalization:
    Te_max = np.max(dd['electrons'].nT_equil['T_J'])
    cs_max = np.sqrt(Te_max / dd['pf'].mass)
    norm_v = 1. / cs_max

    # --- GENE data ---
    gd = data_GENE_CPS2019()
    gvpar_norm = gd['vpar'] * norm_v
    gamma_vmu_aver_mu = np.mean(gd['gamma_vmu'], axis=1)
    gamma_vmu_aver_mu = - gamma_vmu_aver_mu

    # --- ORB5 data ---
    # Initial distribution function:
    f0       = dd_f0['fast'].f_1d['f_vel_1d'][0]
    ovpar_f0 = dd_f0['fast'].f_1d['vpar']

    # <J*E>_mu
    t_point_csR0 = 73.5  # time moment where GENE gamma is taken
    renorm_t_coef = dd['wc'] / (gd['cs']/gd['Lref'])
    t_point = t_point_csR0 * renorm_t_coef

    jdote_es_dict = dd['fast'].jdote_es(dd)
    id_t, _, line_t = mix.get_ids(jdote_es_dict['t'], t_point)

    oje_t1 = np.squeeze(jdote_es_dict['data'][id_t, :, :])
    oje_aver_mu = np.mean(oje_t1, axis=0)
    ovpar = jdote_es_dict['vpar']

    # --- Plot <J*E>_mu for fast species ---
    curves = crv.Curves().xlab(labx_norm).ylab('norm.\ values')\
        .tit('fast\ deuterium')
    curves.flag_norm = True
    curves.new()\
        .XS(gvpar_norm)\
        .YS(gamma_vmu_aver_mu)\
        .leg(leg_gene)
    curves.new() \
        .XS(ovpar) \
        .YS(oje_aver_mu) \
        .leg(leg_orb)
    curves.new() \
        .XS(ovpar_f0) \
        .YS(f0)\
        .sty('--').col('grey').leg('f_0')
    cpr.plot_curves(curves)


# Test histogram:
def test_histogram():
    # # --- TEST plotting of a histogram ---
    # data = np.array([1, 2, 3, 12, 14, 16, 18, 28, 29])
    # data2 = np.array([1, 2, 3, 5, 12, 12, 12, 14, 16, 18, 20])
    # bbins = np.array([0, 10, 20, 30])
    # hist = np.histogram(data, bins=bbins)
    #
    # curves = crv.Curves().xlab('x').ylab('n')
    # curves.new().XS(bbins).YS(data).set_hist().alpha(1).leg('data1')
    # curves.new().XS(10).YS(data2).set_hist().alpha(0.8).leg('data2')
    # cpr.plot_curves(curves)

    # --- TEST fitting to gaussian function ---
    mean_x = 2
    sigma_x = 1
    data_x = np.random.normal(mean_x, sigma_x, 1000)

    n_bins = 20
    hist = np.histogram(data_x, bins=n_bins)
    bins = hist[1]
    hist_res = hist[0]

    fit_mean_x, fit_sigma_x = stat_norm.fit(data_x)
    f_data_x = mlab.normpdf(bins, fit_mean_x, fit_sigma_x)

    test_mean_x = data_x.mean()
    test_sigma_x = data_x.std()

    curves = crv.Curves().xlab('x').ylab('n')
    curves.new().XS(n_bins).YS(data_x).set_hist().alpha(1).leg('data-x')
    curves.new().XS(bins).YS(f_data_x).col('red')
    cpr.plot_curves(curves)


# Calculate the growth/damping rate, using the MPR diagnostic
def calc_gamma(oo_wg, oo_plot, oo_desc):
    # --- DESCRIPTION ---
    # See common.MPR_gamma for description

    def find_gamma(je_loc, ef_loc):
        g_loc = np.mean(- 0.5 * je_loc / ef_loc)
        return g_loc

    # --- Flags ---
    flag_naive_t_peaks = oo_wg.get('flag_naive_t_peaks', True)
    naive_n_periods = oo_wg.get('naive_n_periods', 3)
    w_init = oo_wg.get('gam-w', None)

    flag_t_right = oo_wg.get('flag_t_right', False)
    flag_stat = oo_wg.get('flag_stat', False)

    # --- DATA ---
    je_dict = oo_wg.get('je_dict', None)
    ef_dict = oo_wg.get('ef_dict', None)
    je, ef, t = je_dict['data'], ef_dict['data'], je_dict['x']

    # - domain of plotting -
    t_plot = oo_plot.get('t_plot', [])
    if len(t_plot) is 0:
        t_plot = t
    ids_t_plot, t_plot, _ = mix.get_ids(t, t_plot)

    # - work domain: where intervals will be chosen -
    t_work = oo_wg.get('t_work', [])
    if len(t_work) is 0:
        t_work = t
    ids_t_work, t_work, _ = mix.get_ids(t, t_work)
    je_work = je[ids_t_work]
    ef_work = ef[ids_t_work]

    # - initial domain of calculation -
    ids_je_peaks, _     = scipy.signal.find_peaks(je_work)  # here, DON'T include np.abs
    ids_all_je_peaks, _ = scipy.signal.find_peaks(np.abs(je_work))

    if flag_naive_t_peaks:
        id_last_peak = 2 * np.int(naive_n_periods)

        id_je_peak_start = ids_je_peaks[0]
        id_je_peak_end = ids_je_peaks[id_last_peak]
    else:
        T_init = 2*np.pi/w_init
        t_end  = t_work[0] + naive_n_periods * T_init
        id_je_peak_start = 0
        id_je_peak_end, _, _ = mix.get_ids(t_work, t_end)

    t_calc  = t_work[id_je_peak_start:id_je_peak_end+1]
    je_calc = je_work[id_je_peak_start:id_je_peak_end+1]
    ef_calc = ef_work[id_je_peak_start:id_je_peak_end+1]

    # --- PLOTTING: field energy and energy transfer signal (t) ---
    area_work = Geom.Fill()
    area_work.xs = [t_work[0], t_work[-1], t_work[-1], t_work[0]]
    area_work.ys = ['limb', 'limb', 'limu', 'limu']
    area_work.color = 'grey'
    area_work.alpha = 0.3

    area_calc = Geom.Fill()
    area_calc.xs = [t_calc[0], t_calc[-1], t_calc[-1], t_calc[0]]
    area_calc.ys = ['limb', 'limb', 'limu', 'limu']
    area_calc.color = 'grey'
    area_calc.alpha = 0.6

    # - Ef -
    curves = crv.Curves() \
        .xlab('t[\omega_{ci}^{-1}]').tit(ef_dict['tit'])
    curves.newg(area_work)
    curves.newg(area_calc)
    curves.new() \
        .XS(t_plot).YS(ef[ids_t_plot])
    cpr.plot_curves(curves)

    # - JE -
    curves = crv.Curves() \
        .xlab('t[\omega_{ci}^{-1}]').tit(je_dict['tit'])
    curves.newg(area_work)
    curves.newg(area_calc)
    curves.new() \
        .XS(t_plot).YS(je[ids_t_plot]) \
        .leg(je_dict['line_sum'])
    curves.new() \
        .XS(t_work[ids_je_peaks]).YS(je_work[ids_je_peaks]) \
        .leg('peaks').sty('o')
    cpr.plot_curves(curves)

    # --- Trivial calculation of damping/growth rate ---
    g_init = find_gamma(je_calc, ef_calc)

    # --- VARIATION of the RIGHT BOUNDARY ---
    if flag_t_right:
        n_points = oo_wg.get('right-n', None)
        gs_right = np.zeros(2*n_points + 1)
        ts_right = np.zeros(2*n_points + 1)

        id_je_peak_middle = id_je_peak_end

        ids_peaks_right = range(
            id_je_peak_middle - n_points,
            id_je_peak_middle + n_points + 1
        )

        g_middle = None
        count_peak = -1
        for id_peak_right in ids_peaks_right:
            count_peak += 1

            je_current = je_work[id_je_peak_start:id_peak_right]
            ef_current = ef_work[id_je_peak_start:id_peak_right]

            gs_right[count_peak] = find_gamma(je_current, ef_current)
            ts_right[count_peak] = t_work[id_peak_right]

            if id_peak_right == id_je_peak_middle:
                g_middle = gs_right[count_peak]

        ids_all_right_peaks = ids_all_je_peaks[6:]
        gs_all = np.zeros(len(ids_all_right_peaks))
        count_peak = -1
        for id_peak_right in ids_all_right_peaks:
            count_peak += 1

            je_current = je_work[id_je_peak_start:id_peak_right]
            ef_current = ef_work[id_je_peak_start:id_peak_right]

            gs_all[count_peak] = find_gamma(je_current, ef_current)

        curves = crv.Curves().xlab('right\ t').ylab('g[\omega_{ci}]')
        curves.new().XS(ts_right).YS(gs_right)
        curves.new().XS(t_work[ids_all_right_peaks]).YS(gs_all).sty('o')
        curves.new().XS(t_work[id_je_peak_middle]).YS(g_middle).sty('s')
        cpr.plot_curves(curves)

        curves = crv.Curves().xlab('right\ t').ylab('JE')
        curves.new().XS(t_work[ids_peaks_right]).YS(je_work[ids_peaks_right])
        curves.new().XS(t_work[id_je_peak_middle]).YS(je_work[id_je_peak_middle]).sty('s')
        cpr.plot_curves(curves)

    # --- STATISTICS ---
    stat_res = {}
    if flag_stat:
        n_samples         = oo_wg.get('n_samples', None)
        min_gam_n_periods = oo_wg.get('min_gam_n_periods', None)

        gam_T_init = 2 * np.pi / w_init

        flag_print_stat_full = False
        if n_samples <= 10:
            flag_print_stat_full = True

        dict_intervals = mix.get_t_intervals(
            n_samples, t_work, min_gam_n_periods, gam_T_init
        )

        t_intervals     = dict_intervals['t_intervals']
        ids_t_intervals = dict_intervals['ids_intervals']

        if flag_print_stat_full:
            for one_t_interval in t_intervals:
                area_calc_chosen = Geom.Fill()
                area_calc_chosen.xs = [
                    one_t_interval[0], one_t_interval[-1],
                    one_t_interval[-1], one_t_interval[0]
                ]
                area_calc_chosen.ys = ['limb', 'limb', 'limu', 'limu']
                area_calc_chosen.color = 'grey'
                area_calc_chosen.alpha = 0.6

                curves = crv.Curves() \
                    .xlab('t[\omega_{ci}^{-1}]').tit('STAT.:\ ' + je_dict['tit'])
                curves.newg(area_work)
                curves.newg(area_calc_chosen)
                curves.new() \
                    .XS(t_plot).YS(je[ids_t_plot]) \
                    .leg(je_dict['line_sum'])
                curves.new() \
                    .XS(t_work[ids_je_peaks]).YS(je_work[ids_je_peaks]) \
                    .leg('peaks').sty('o')

        # calculate gamma in all time intervals:
        nsamples = len(t_intervals)
        gs = np.zeros(nsamples)
        for count_sample in range(nsamples):
            ids_one_t_interval = ids_t_intervals[count_sample]
            # noinspection PyTypeChecker
            je_calc = je_work[
                ids_one_t_interval[0]:ids_one_t_interval[-1] + 1
            ]
            # noinspection PyTypeChecker
            ef_calc = ef_work[
                ids_one_t_interval[0]:ids_one_t_interval[-1] + 1
            ]
            gs[count_sample] = find_gamma(je_calc, ef_calc)

            if flag_print_stat_full:
                tx_left  = t_work[ids_one_t_interval[0]]
                tx_right = t_work[ids_one_t_interval[-1]]
                n_pers = (tx_right - tx_left) / gam_T_init
                print('TEST: g = {:0.3e}, ngam = {:f}, t1 = {:0.3e}, t2 = {:0.3e}'.format(
                    gs[count_sample], n_pers, tx_left, tx_right)
                )

        # build a histogram
        n_bins = 40
        hist = np.histogram(gs, bins=n_bins)
        bins = hist[1]

        fit_mean_x, fit_sigma_x = stat_norm.fit(gs)
        f_data_x = mlab.normpdf(bins, fit_mean_x, fit_sigma_x)

        stat_res = {
            'gamma': fit_mean_x,
            'conf_int_95': 1.96 * fit_sigma_x,
            'old_errorbar': 3 * fit_sigma_x
        }

        curves = crv.Curves().xlab('\gamma[\omega_{ci}]').ylab('n')
        curves.new().XS(n_bins).YS(gs).set_hist().alpha(1).leg('\gamma[\omega_{ci}]')
        curves.new().XS(bins).YS(f_data_x).col('red')
        cpr.plot_curves(curves)

    # --- PRINTING RESULT ---
    # - trivial result -
    sp_name = oo_desc['species'] + ': '
    line_init = sp_name + 'Initial: '
    line_init += 'g[wci] = ' + '{:0.3e}'.format(g_init)
    print(line_init)

    # - with statistics -
    if flag_stat:
        line_stat = sp_name + 'Statistics: '
        line_stat += 'g[wci] = ' + '{:0.3e}'.format(stat_res['gamma'])\
                     + '+-' + '{:0.3e}'.format(stat_res['conf_int_95']) + '\n'
        line_stat += 'OLD confidence interval (3*sigma) = ' \
                     + '{:0.3e}'.format(stat_res['old_errorbar'])
        print(line_stat)

    return


