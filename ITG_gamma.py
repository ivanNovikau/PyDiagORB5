import Mix as mix
import read_data as rd
import ControlPlot as cpr
import ymath
import curve as crv
import zf_gam as zf
import write_data as wr
import transport
import numpy as np
from scipy import interpolate
import h5py as h5


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)
    mix.reload_module(rd)
    mix.reload_module(cpr)
    mix.reload_module(ymath)
    mix.reload_module(crv)
    mix.reload_module(zf)
    mix.reload_module(wr)
    mix.reload_module(transport)


# radial derivative of the full potential at a particular chi:
def ersc(dd, oo):
    chi_s = oo.get('chi_s', [0.0])

    names_potsc = rd.potsc_chi(dd, oo)

    nchi = np.size(chi_s)
    names_ersc_chi = []
    for count_chi in range(nchi):
        one_chi       = chi_s[count_chi]
        name_ersc_chi = 'ersc-chi-' + '{:0.3f}'.format(one_chi)
        names_ersc_chi.append(name_ersc_chi)
        if name_ersc_chi in dd:
            continue

        one_name_potsc = names_potsc[count_chi]
        data = - np.gradient(
            dd[one_name_potsc]['data'], dd['potsc_grids']['s'], axis=1)
        dd[name_ersc_chi] = {
            'chi_1':    dd[one_name_potsc]['chi_1'],
            'id_chi_1': dd[one_name_potsc]['id_chi_1'],
            'data': data}
    return names_ersc_chi


# phibar along potsc grids:
def phibar_interp(dd):
    if 'phibar_interp' in dd:
        return

    zf.phibar(dd)
    rd.potsc_grids(dd)

    t = dd['phibar']['t']
    s = dd['phibar']['s']

    f_interp = interpolate.interp2d(t, s, dd['phibar']['data'].T)
    dd['phibar_interp'] = {
        'data': f_interp(dd['potsc_grids']['t'], dd['potsc_grids']['s']).T
    }


# non-zonal potential at some poloidal angles
def phinz(dd, oo):
    chi_s = oo.get('chi_s', [0.0])

    phibar_interp(dd)
    names_potsc = rd.potsc_chi(dd, oo)

    nchi = np.size(chi_s)
    names_phinz_chi = []
    for count_chi in range(nchi):
        one_chi        = chi_s[count_chi]
        name_phinz_chi = 'phinz-chi-' + '{:0.3f}'.format(one_chi)
        names_phinz_chi.append(name_phinz_chi)
        if name_phinz_chi in dd:
            continue

        one_name_potsc = names_potsc[count_chi]
        data = dd[one_name_potsc]['data'] - dd['phibar_interp']['data']
        dd[name_phinz_chi] = {
            'chi_1':    dd[one_name_potsc]['chi_1'],
            'id_chi_1': dd[one_name_potsc]['id_chi_1'],
            'data': data}
    return names_phinz_chi


# non-zonal potential at a time point
def phinz_t(dd, oo):
    names_potsc = rd.potsc_t(dd, oo)
    phibar_interp(dd)

    t_points = oo.get('t_points', [0.0])
    nt_points = np.size(t_points)

    names = []
    for count_t in range(nt_points):
        one_t = t_points[count_t]
        var_name = 'phinz-t-' + '{:0.3e}'.format(one_t)
        names.append(var_name)
        if var_name in dd:
            continue

        data = {}
        one_name_potsc = names_potsc[count_t]
        data['id_t_point'], data['t_point'] = \
            mix.find(dd['potsc_grids']['t'], one_t)
        phibar_t1   = dd['phibar_interp']['data'][data['id_t_point'], :]
        data['data'] = dd[one_name_potsc]['data'] - phibar_t1[None, :]
        dd[var_name] = data
    return names


# non-zonal potential at some radial point:
def phinz_s(dd, oo):
    names_potsc = rd.potsc_s(dd, oo)
    phibar_interp(dd)

    s_points = oo.get('s_points', [0.0])
    ns_points = np.size(s_points)

    names = []
    for count_s in range(ns_points):
        one_s = s_points[count_s]
        var_name = 'phinz-s-' + '{:0.3e}'.format(one_s)
        names.append(var_name)
        if var_name in dd:
            continue

        data = {}
        one_name_potsc = names_potsc[count_s]
        data['id_s_point'] = dd[one_name_potsc]['id_s_1']
        data['s_point']    = dd[one_name_potsc]['s_point']
        phibar_s1    = dd['phibar_interp']['data'][:, data['id_s_point']]
        data['data'] = dd[one_name_potsc]['data'] - phibar_s1[:, None]
        dd[var_name] = data
    return names


#  radial derivative of the non-zonal potential at a particular chi:
def ernz_r(dd, oo):
    chi_s = oo.get('chi_s', [0.0])

    names_phinz = phinz(dd, oo)

    nchi = np.size(chi_s)
    names_ernz_chi = []

    for count_chi in range(nchi):
        one_chi       = chi_s[count_chi]
        name_ernz_chi = 'ernz_r-chi-{:0.3f}'.format(one_chi)
        names_ernz_chi.append(name_ernz_chi)
        if name_ernz_chi in dd:
            continue

        one_name_phinz = names_phinz[count_chi]
        data = - np.gradient(
            dd[one_name_phinz]['data'], dd['potsc_grids']['s'], axis=1)
        dd[name_ernz_chi] = {
            'chi_1':    dd[one_name_phinz]['chi_1'],
            'id_chi_1': dd[one_name_phinz]['id_chi_1'],
            'data': data}
    return names_ernz_chi


# poloidal derivative of the non-zonal potential at a particular chi:
def ernz_chi(dd, oo):
    phibar_interp(dd)

    chi_s = oo.get('chi_s', [0.0])
    nchi_points = np.size(chi_s)
    t, s, chi = dd['potsc_grids']['t'], dd['potsc_grids']['s'], dd['potsc_grids']['chi']
    nt, ns, nchi = np.size(t), np.size(s), np.size(chi)

    # number of radial points to read at once
    max_allowed_size = dd['max_size_Gb'] * 1024 ** 3  # in bytes
    float_size = 8  # in bytes
    id_s_step = int( round( max_allowed_size / (float_size * nt * nchi) ) )

    names = []
    for count_chi in range(nchi_points):
        one_chi       = chi_s[count_chi]
        var_name = 'ernz_chi-chi-{:0.3f}'.format(one_chi)
        names.append(var_name)
        if var_name in dd:
            continue

        data = rd.read_signal(dd['path_ext'], var_name)
        if data is None:
            data = {}

            # --- calculate data ---
            f = h5.File(dd['path_orb'], 'r')
            data['data'] = np.zeros([ns, nt])  # transposed w.r.t final matrix
            data['id_chi_1'], data['chi_1'] = mix.find(chi, one_chi)
            line_read = '/data/var2d/generic/potsc/data'
            id_s_end = 0
            while id_s_end < ns:
                # read Phi at several points along s: Phi(t,chi,ids_current)
                id_s_begin = id_s_end
                id_s_end  += id_s_step
                if id_s_end >= ns:
                    id_s_end = ns
                ids_current = [i for i in range(id_s_begin, id_s_end)]
                potsc_s1 = np.array(f[line_read][:, :, ids_current])

                # find non-zonal Phi(t,chi,ids_current)
                phibar_s1 = dd['phibar_interp']['data'][:, ids_current]
                phinz_s1 = potsc_s1 - phibar_s1[:, None, :]
                del potsc_s1

                # chi-derivative of the non-zonal Phi(t,chi,ids_current)
                loc = np.gradient(phinz_s1, chi, axis=1)
                del phinz_s1

                # save non-zonal poloidal electric field E_chi(t,chi1,ids_current)
                count_id_s = -1
                for id_s_loc in ids_current:
                    count_id_s += 1
                    data['data'][id_s_loc] = \
                        - loc[:, data['id_chi_1'], count_id_s] / s[id_s_loc]
                del loc
            f.close()

            # for the sake of generality, transpose data (ns,nt) -> (nt,ns)
            data['data'] = data['data'].T

            # save the data to an external file
            desc = 'non-zonal poloidal electric field at chi = {:0.3f}'.format(data['chi_1'])
            wr.save_data_adv(dd['path_ext'], data, {'name': var_name, 'desc': desc})

        # --- save data to the structure ---
        dd[var_name] = data
    return names


# NEW: max of the non-zonal potential along chi at every s-point:
def phinz_abs_max_along_chi(dd):
    phibar_interp(dd)
    t, s, chi = dd['potsc_grids']['t'], dd['potsc_grids']['s'], dd['potsc_grids']['chi']
    nt, ns, nchi = np.size(t), np.size(s), np.size(chi)

    # --- check the data in the project structure ---
    var_name = 'phinz-max-along-chi'
    if var_name in dd:
        return var_name

    # --- check the data in the external file ---
    data = rd.read_signal(dd['path_ext'], var_name)
    if data is not None:
        dd[var_name] = data
        return var_name

    # number of radial points to read at once
    max_allowed_size = dd['max_size_Gb'] * 1024 ** 3  # in bytes
    float_size = 8  # in bytes
    id_s_step = int(round(max_allowed_size / (float_size * nt * nchi)))

    # --- calculate the data ---
    f = h5.File(dd['path_orb'], 'r')
    data = {
        's': s, 't': t, 'data': np.zeros([ns, nt]),
        'chi_max': np.zeros([ns, nt])
    }
    line_read = '/data/var2d/generic/potsc/data'
    id_s_end = 0
    while id_s_end < ns:
        # read Phi at several points along s: Phi(t,chi,ids_current)
        id_s_begin = id_s_end
        id_s_end += id_s_step
        if id_s_end >= ns:
            id_s_end = ns
        ids_current = [i for i in range(id_s_begin, id_s_end)]

        # read Phi at several points along s: Phi(t,chi,ids_current)
        potsc_s1 = np.array(f[line_read][:, :, ids_current])

        # find non-zonal Phi(t,chi,ids_current)
        phibar_s1 = dd['phibar_interp']['data'][:, ids_current]
        phinz_s1 = potsc_s1 - phibar_s1[:, None, :]
        del potsc_s1, phibar_s1

        # find absolute value
        abs_phinz_s1 = np.abs(phinz_s1)
        del phinz_s1

        # find maximum
        count_id_s = -1
        for id_s_loc in ids_current:
            count_id_s += 1
            phinz_max = np.amax(abs_phinz_s1[:, :, count_id_s], axis=1)
            data['data'][id_s_loc] = phinz_max

            ids_max = abs_phinz_s1[:, :, count_id_s].argmax(axis=1)
            data['chi_max'][id_s_loc] = chi[ids_max]

    # for the sake of generality, transpose data (ns,nt) -> (nt,ns)
    data['data']    = data['data'].T
    data['chi_max'] = data['chi_max'].T

    f.close()

    # save data to an external file
    desc = 'max along chi of non-zonal potential at phi = 0'
    wr.save_data_adv(dd['path_ext'], data, {'name': var_name, 'desc': desc})

    # save to the project structure
    dd[var_name] = data

    return var_name


#  NEW: radial derivative of the nonzonal potential at s1:
def ernz_r_abs_max_along_chi(dd):
    phibar_interp(dd)
    t,   s,  chi = dd['potsc_grids']['t'], dd['potsc_grids']['s'], dd['potsc_grids']['chi']
    nt, ns, nchi = np.size(t), np.size(s), np.size(chi)

    # --- check the data in the project structure ---
    var_name = 'ernz_r-max-along-chi'
    if var_name in dd:
        return var_name

    # --- check the data in the external file ---
    data = rd.read_signal(dd['path_ext'], var_name)
    if data is not None:
        dd[var_name] = data
        return var_name

    # number of poloidal points to read at once
    max_allowed_size = dd['max_size_Gb'] * 1024 ** 3  # in bytes
    float_size = 8  # in bytes
    id_chi_step = int(round(max_allowed_size / (float_size * nt * ns)))

    # --- calculate data ---
    data = {
        's': s, 't': t, 'data': [],
        'chi_max': np.zeros([nt, ns])
    }

    f = h5.File(dd['path_orb'], 'r')
    line_read = '/data/var2d/generic/potsc/data'
    id_chi_end = 0
    global_max     = np.zeros([nt, ns])
    global_ids_max = np.zeros([nt, ns])
    while id_chi_end < nchi:
        # read Phi at several points along s: Phi(t,chi,ids_current)
        id_chi_begin = id_chi_end
        id_chi_end  += id_chi_step
        if id_chi_end >= nchi:
            id_chi_end = nchi
        ids_current = [i for i in range(id_chi_begin, id_chi_end)]
        potsc_chi1 = np.array(f[line_read][:, ids_current, :])

        # find non-zonal Phi(t,ids_current,s)
        phibar_data = dd['phibar_interp']['data'][:, :]
        phinz_chi1  = potsc_chi1 - phibar_data[:, None, :]
        del potsc_chi1, phibar_data

        # chi-derivative of the non-zonal Phi(t,ids_current,s)
        loc = - np.gradient(phinz_chi1, s, axis=2)
        del phinz_chi1

        # find absolute value
        abs_loc = np.abs(loc)
        del loc

        # find maximum in a given chi-interval
        current_max     = np.amax(abs_loc[:, :, :], axis=1)
        current_ids_max = abs_loc[:, :, :].argmax(axis=1)
        del abs_loc

        # correct global maximum
        mask_comp = current_max > global_max
        global_max     = np.where(mask_comp, current_max, global_max)
        global_ids_max = np.where(mask_comp, current_ids_max, global_ids_max)

    f.close()
    data['data']    = global_max

    for id_s in range(ns):
        for id_t in range(nt):
            data['chi_max'][id_t, id_s] = chi[int(global_ids_max[id_t, id_s])]

    # save the data to an external file
    desc = 'max along chi of non-zonal radial electric field at phi = 0'
    wr.save_data_adv(dd['path_ext'], data, {'name': var_name, 'desc': desc})

    # --- save data to the structure ---
    dd[var_name] = data
    return var_name


# NEW: nonzonal phi at chi1
def phinz_chi1(dd, chi_point):
    phibar_interp(dd)
    name_potsc = rd.potsc_chi1(dd, chi_point)

    name_phinz_chi = 'phinz-chi-' + '{:0.3f}'.format(chi_point)
    if name_phinz_chi in dd:
        return name_phinz_chi

    data = dd[name_potsc]['data'] - dd['phibar_interp']['data']
    dd[name_phinz_chi] = {
        'chi_1':    dd[name_potsc]['chi_1'],
        'id_chi_1': dd[name_potsc]['id_chi_1'],
        'data': data}
    return name_phinz_chi


# NEW: nonzonal phi at s1
def phinz_s1(dd, s_point):
    phibar_interp(dd)
    name_potsc = rd.potsc_s1(dd, s_point)

    name_phinz = 'phinz-s-' + '{:0.3f}'.format(s_point)
    if name_phinz in dd:
        return name_phinz

    data = dd[name_potsc]['data'] - dd['phibar_interp']['data'][:, dd[name_potsc]['id_s1'], None]
    dd[name_phinz] = {
        's1':    dd[name_potsc]['s1'],
        'id_s1': dd[name_potsc]['id_s1'],
        'data':  data}
    return name_phinz


# NEW: nonzonal phi at t1
def phinz_t1(dd, t_point):
    phibar_interp(dd)
    name_potsc = rd.potsc_t1(dd, t_point)

    name_phinz = 'phinz-t-' + '{:0.3e}'.format(t_point)
    if name_phinz in dd:
        return name_phinz

    phibar_t1 = dd['phibar_interp']['data'][dd[name_potsc]['id_t1'], :]
    data      = dd[name_potsc]['data'] - phibar_t1[None, :]
    dd[name_phinz] = {
        't1':    dd[name_potsc]['t1'],
        'id_t1': dd[name_potsc]['id_t1'],
        'data':  data  # (chi, s)
    }
    return name_phinz


#  NEW: radial derivative of the nonzonal potential at chi1:
def ernz_r_chi1(dd, chi_point):
    name_phinz = phinz_chi1(dd, chi_point)
    name_ernz_chi = 'ernz_r-chi-{:0.3f}'.format(chi_point)
    if name_ernz_chi in dd:
        return name_ernz_chi

    data = - np.gradient(
        dd[name_phinz]['data'], dd['potsc_grids']['s'], axis=1)
    dd[name_ernz_chi] = {
        'chi_1':    dd[name_phinz]['chi_1'],
        'id_chi_1': dd[name_phinz]['id_chi_1'],
        'data': data}
    return name_ernz_chi


#  NEW: radial derivative of the full potential at chi1:
def er_r_chi1(dd, chi_point):
    name_phi = rd.potsc_chi1(dd, chi_point)
    name_er_chi = 'er_r-chi-{:0.3f}'.format(chi_point)
    if name_er_chi in dd:
        return name_er_chi

    data = - np.gradient(
        dd[name_phi]['data'], dd['potsc_grids']['s'], axis=1)
    dd[name_er_chi] = {
        'chi_1':    dd[name_phi]['chi_1'],
        'id_chi_1': dd[name_phi]['id_chi_1'],
        'data': data}
    return name_er_chi


#  NEW: radial derivative of the nonzonal potential at s1:
def ernz_r_s1(dd, s_point):
    phibar_interp(dd)
    t, s, chi    = dd['potsc_grids']['t'], dd['potsc_grids']['s'], dd['potsc_grids']['chi']
    nt, ns, nchi = np.size(t), np.size(s), np.size(chi)

    # number of poloidal points to read at once
    max_allowed_size = dd['max_size_Gb'] * 1024 ** 3  # in bytes
    float_size = 8  # in bytes
    id_chi_step = int(round(max_allowed_size / (float_size * nt * ns)))

    var_name = 'ernz_r-s-{:0.3f}'.format(s_point)
    if var_name in dd:
        return var_name

    data = rd.read_signal(dd['path_ext'], var_name)
    if data is None:
        data = {}

        # --- calculate data ---
        f = h5.File(dd['path_orb'], 'r')
        data['data'] = np.zeros([nchi, nt])  # transposed w.r.t final matrix
        data['id_s1'], data['s1'] = mix.find(s, s_point)
        line_read = '/data/var2d/generic/potsc/data'
        id_chi_end = 0
        while id_chi_end < nchi:
            # read Phi at several points along s: Phi(t,chi,ids_current)
            id_chi_begin = id_chi_end
            id_chi_end  += id_chi_step
            if id_chi_end >= nchi:
                id_chi_end = nchi
            ids_current = [i for i in range(id_chi_begin, id_chi_end)]
            potsc_chi1 = np.array(f[line_read][:, ids_current, :])

            # find non-zonal Phi(t,ids_current,s)
            phibar_data = dd['phibar_interp']['data'][:, :]
            phinz_chi1 = potsc_chi1 - phibar_data[:, None, :]
            del potsc_chi1, phibar_data

            # chi-derivative of the non-zonal Phi(t,ids_current,s)
            loc = np.gradient(phinz_chi1, s, axis=2)
            del phinz_chi1

            # save non-zonal poloidal electric field E_chi(t,ids_current,s)
            count_id_chi = -1
            for id_chi_loc in ids_current:
                count_id_chi += 1
                data['data'][id_chi_loc] = - loc[:, count_id_chi, data['id_s1']]
            del loc
        f.close()

        # for the sake of generality, transpose data (nchi,nt) -> (nt,nchi)
        data['data'] = data['data'].T

        # save the data to an external file
        desc = 'non-zonal poloidal electric field at s = {:0.3f}'.format(data['s1'])
        wr.save_data_adv(dd['path_ext'], data, {'name': var_name, 'desc': desc})

    # --- save data to the structure ---
    dd[var_name] = data
    return var_name


# NEW: poloidal derivative of the non-zonal potential at chi1:
def ernz_chi_chi1(dd, chi_point):
    phibar_interp(dd)
    t, s, chi = dd['potsc_grids']['t'], dd['potsc_grids']['s'], dd['potsc_grids']['chi']
    nt, ns, nchi = np.size(t), np.size(s), np.size(chi)

    # number of radial points to read at once
    max_allowed_size = dd['max_size_Gb'] * 1024 ** 3  # in bytes
    float_size = 8  # in bytes
    id_s_step = int( round( max_allowed_size / (float_size * nt * nchi) ) )

    var_name = 'ernz_chi-chi-{:0.3f}'.format(chi_point)
    if var_name in dd:
        return var_name

    data = rd.read_signal(dd['path_ext'], var_name)
    if data is None:
        data = {}

        # --- calculate data ---
        f = h5.File(dd['path_orb'], 'r')
        data['data'] = np.zeros([ns, nt])  # transposed w.r.t final matrix
        data['id_chi_1'], data['chi_1'] = mix.find(chi, chi_point)
        line_read = '/data/var2d/generic/potsc/data'
        id_s_end = 0
        while id_s_end < ns:
            # read Phi at several points along s: Phi(t,chi,ids_current)
            id_s_begin = id_s_end
            id_s_end  += id_s_step
            if id_s_end >= ns:
                id_s_end = ns
            ids_current = [i for i in range(id_s_begin, id_s_end)]
            potsc_s1 = np.array(f[line_read][:, :, ids_current])

            # find non-zonal Phi(t,chi,ids_current)
            phibar_s1 = dd['phibar_interp']['data'][:, ids_current]
            phinz_s1 = potsc_s1 - phibar_s1[:, None, :]
            del potsc_s1

            # chi-derivative of the non-zonal Phi(t,chi,ids_current)
            loc = np.gradient(phinz_s1, chi, axis=1)
            del phinz_s1

            # save non-zonal poloidal electric field E_chi(t,chi1,ids_current)
            count_id_s = -1
            for id_s_loc in ids_current:
                count_id_s += 1
                data['data'][id_s_loc] = \
                    - loc[:, data['id_chi_1'], count_id_s] / s[id_s_loc]
            del loc
        f.close()

        # for the sake of generality, transpose data (ns,nt) -> (nt,ns)
        data['data'] = data['data'].T

        # save the data to an external file
        desc = 'non-zonal poloidal electric field at chi = {:0.3f}'.format(data['chi_1'])
        wr.save_data_adv(dd['path_ext'], data, {'name': var_name, 'desc': desc})

    # --- save data to the structure ---
    dd[var_name] = data
    return var_name


def choose_one_var_ts(one_signal):
    dd = one_signal['dd']
    opt_var = one_signal['variable']

    var_name, tit_var, line_chi = '', '', ''
    res = {}
    type_potsc = 'potsc'
    if opt_var == 'phinz':
        chi_point = one_signal['chi-point']
        var_name = phinz_chi1(dd, chi_point)
        tit_var  = '\widetilde{\Phi}'
        line_chi = '\chi = {:0.1f}'.format(dd[var_name]['chi_1'])
        line_chi = '_{' + line_chi + '}'
    elif opt_var == 'potsc':
        chi_point = one_signal['chi-point']
        var_name = rd.potsc_chi1(dd, chi_point)
        tit_var  = '\Phi'
        line_chi = '\chi = {:0.1f}'.format(dd[var_name]['chi_1'])
        line_chi = '_{' + line_chi + '}'
    elif opt_var == 'potsc-antenna':
        type_potsc = 'potsc_ant'
        chi_point = one_signal['chi-point']
        var_name = rd.potsc_chi1(dd, chi_point, type=type_potsc)
        tit_var = '\Phi'
        line_chi = '\chi = {:0.1f}'.format(dd[var_name]['chi_1'])
        line_chi = '_{ant, ' + line_chi + '}'
    elif opt_var == 'ernz_r':
        chi_point = one_signal['chi-point']
        var_name = ernz_r_chi1(dd, chi_point)
        tit_var  = '\widetilde{E}'
        line_chi = '\chi = {:0.1f}'.format(dd[var_name]['chi_1'])
        line_chi = '_{r,' + line_chi + '}'
    elif opt_var == 'er_r':  # full radial electric field along radial coordinate
        chi_point = one_signal['chi-point']
        var_name = er_r_chi1(dd, chi_point)
        tit_var  = 'E'
        line_chi = '\chi = {:0.1f}'.format(dd[var_name]['chi_1'])
        line_chi = '_{r, ' + line_chi + '}'
    elif opt_var == 'ernz_chi':
        chi_point = one_signal['chi-point']
        var_name = ernz_chi_chi1(dd, chi_point)
        tit_var  = '\widetilde{E}'
        line_chi = '\chi = {:0.1f}'.format(dd[var_name]['chi_1'])
        line_chi = '_{\chi,' + line_chi + '}'
    elif opt_var == 'phinz-max-chi':
        var_name = phinz_abs_max_along_chi(dd)
        tit_var  = 'max_{\chi}:\ \widetilde{\Phi}'
        line_chi = ''
        res.update({
            'chi_max': dd[var_name]['chi_max']
        })
    elif opt_var == 'ernz_r-max-chi':
        var_name = ernz_r_abs_max_along_chi(dd)
        tit_var  = 'max_{\chi}:\ \widetilde{E}_r'
        line_chi = ''
        res.update({
            'chi_max': dd[var_name]['chi_max']
        })
    else:
        mix.error_mes('Wrong name of signal')

    vvar = dd[var_name]['data']
    s    = dd[type_potsc + '_grids']['s']
    t    = dd[type_potsc + '_grids']['t']

    tit_var += line_chi

    res.update({
        'data': vvar,
        's': s,
        't': t,
        'tit': tit_var
    })

    return res


def choose_one_var_rz(one_signal):
    dd = one_signal['dd']
    opt_var = one_signal['variable']

    var_name, tit_var, line_1 = '', '', ''
    res = {}
    type_potsc = 'potsc'
    if opt_var == 'phinz':
        t_point = one_signal['t-point']
        var_name = phinz_t1(dd, t_point)
        tit_var  = '\widetilde{\Phi}'
        line_1 = 't = {:0.3e}'.format(dd[var_name]['t1'])
        line_1 = '_{' + line_1 + '}'
    elif opt_var == 'potsc':
        t_point = one_signal['t-point']
        var_name = rd.potsc_t1(dd, t_point)
        tit_var = '\Phi'
        line_1 = 't = {:0.3e}'.format(dd[var_name]['t1'])
        line_1 = '_{' + line_1 + '}'
    elif opt_var == 'potsc-antenna':
        type_potsc = 'potsc_ant'
        t_point = one_signal['t-point']
        var_name = rd.potsc_t1(dd, t_point, type=type_potsc)
        tit_var = '\Phi'
        line_1 = 't = {:0.3e}'.format(dd[var_name]['t1'])
        line_1 = '_{ant, ' + line_1 + '}'
    else:
        mix.error_mes('Wrong name of a signal.')

    vvar = dd[var_name]['data']
    R    = dd[type_potsc + '_grids']['r']
    Z    = dd[type_potsc + '_grids']['z']
    s    = dd[type_potsc + '_grids']['s']
    chi  = dd[type_potsc + '_grids']['chi']

    tit_var += line_1

    res.update({
        'data': vvar.T,  # (s, chi)
        'r': R,
        'z': Z,
        's': s,
        'chi': chi,
        'tit': tit_var
    })

    return res


def choose_one_var_tchi(one_signal):
    dd = one_signal['dd']
    opt_var = one_signal['variable']
    s_point = one_signal['s-point']

    var_name, tit_var, line_chi = '', '', ''
    type_potsc = 'potsc'
    if opt_var == 'potsc':
        var_name = rd.potsc_s1(dd, s_point)
        tit_var  = '\Phi'
        line_chi = 's = {:0.3f}'.format(dd[var_name]['s1'])
        line_chi = '_{' + line_chi + '}'
    elif opt_var == 'potsc-antenna':
        type_potsc = 'potsc_ant'
        var_name = rd.potsc_s1(dd, s_point, type=type_potsc)
        tit_var = '\Phi'
        line_chi = 's = {:0.3f}'.format(dd[var_name]['s1'])
        line_chi = '_{ant, ' + line_chi + '}'
    elif opt_var == 'phinz':
        var_name = phinz_s1(dd, s_point)
        tit_var  = '\widetilde{\Phi}'
        line_chi = 's = {:0.3f}'.format(dd[var_name]['s1'])
        line_chi = '_{' + line_chi + '}'
    elif opt_var == 'ernz_r':
        var_name = ernz_r_s1(dd, s_point)
        tit_var  = '\widetilde{E}_r'
        line_chi = 's = {:0.3f}'.format(dd[var_name]['s1'])
        line_chi = '_{' + line_chi + '}'
    else:
        mix.error_mes('Wrong name of a signal.')

    vvar = dd[var_name]['data']
    chi  = dd[type_potsc + '_grids']['chi']
    t    = dd[type_potsc + '_grids']['t']

    tit_var += line_chi

    res = {
        'data': vvar,
        'chi': chi,
        't': t,
        'tit': tit_var
    }

    return res








