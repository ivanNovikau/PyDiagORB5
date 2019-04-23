import Mix as mix
import numpy as np
import scipy.signal
from scipy import stats
from scipy.optimize import curve_fit
from scipy import constants
from scipy import interpolate
from scipy.fftpack import fft as sci_fft

#-> avoid plotting here


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)


def estimate_w(x_init, y_init, oo={}):
    # additional parameters
    _, ids_x = mix.get_array_oo(oo, x_init, 'x')
    x = np.array(mix.get_slice(x_init, ids_x))
    y = np.array(mix.get_slice(y_init, ids_x))

    # find peaks of the signal
    ids_peaks, _ = scipy.signal.find_peaks(abs(y), height=0)
    if np.size(ids_peaks) is 0:
        return None

    x_peaks = x[ids_peaks]
    y_peaks = abs(y[ids_peaks])

    # frequency
    w = 2 * np.pi * 0.5 / np.mean(np.diff(x_peaks))

    x_fit = x[ids_peaks[0]:ids_peaks[-1] + 1] - x[ids_peaks[0]]
    y_fit = y_peaks[0] * np.cos(w * x_fit)

    # results
    out = {'x_peaks': x_peaks, 'y_peaks': y_peaks,
           'w': w,
           'y_fit': y_fit, 'x_fit': x_fit + x[ids_peaks[0]]}
    return out


def estimate_wg(x_init, y_init, oo={}):
    # additional parameters
    _, ids_x = mix.get_array_oo(oo, x_init, 'x')
    x = np.array(mix.get_slice(x_init, ids_x))
    y = np.array(mix.get_slice(y_init, ids_x))

    # find peaks of the signal
    ids_peaks, _ = scipy.signal.find_peaks(abs(y), height=0)
    if np.size(ids_peaks) is 0:
        return None

    x_peaks = x[ids_peaks]
    y_peaks = abs(y[ids_peaks])

    # frequency
    w = 2 * np.pi * 0.5 / np.mean(np.diff(x_peaks))
    g, b0,  r_value, p_value, std_err = stats.linregress(
                                        x_peaks - x[ids_peaks[0]],
                                        np.log(y_peaks)
                                    )

    x_fit = x[ids_peaks[0]:ids_peaks[-1]+1] - x[ids_peaks[0]]
    y_fit = np.cos(w * x_fit) * np.exp(b0 + g * x_fit)

    # results
    out = {'x_peaks': x_peaks, 'y_peaks': y_peaks,
          'w': w, 'g': g,
          'y_fit': y_fit, 'x_fit': x_fit + x[ids_peaks[0]]}
    return out


def advanced_wg(x_init, y_init, oo={}):
    # additional parameters
    _, ids_x = mix.get_array_oo(oo, x_init, 'x')
    x = np.array(mix.get_slice(x_init, ids_x))
    y = np.array(mix.get_slice(y_init, ids_x))

    # estimation of w,g
    wg_est = oo.get('est', None)
    if wg_est is None:
        wg_est = estimate_wg(x, y)

    # test function
    F = lambda x_loc, A0, w, g: A0 * np.cos(w*x_loc) * np.exp(g*x_loc)

    # initial guess
    p0 = [y[0], wg_est['w'], wg_est['g']]

    # nonlinear fitting
    popt, pcov = curve_fit(F, x - x[0], y, p0)

    # fitted function
    y_fit = F(x - x[0], *popt)

    # results
    out = {'est': wg_est,
           'w': popt[1], 'g': popt[2],
           'y_fit': y_fit, 'x_fit': x}
    return out


def estimate_g(x_init, y_init, oo={}):
    # additional parameters
    _, ids_x = mix.get_array_oo(oo, x_init, 'x')
    x = np.array(mix.get_slice(x_init, ids_x))
    y = np.array(mix.get_slice(y_init, ids_x))

    # initialize output:
    out = {}

    # find peaks of the signal
    ids_peaks, _ = scipy.signal.find_peaks(abs(y), height=0)
    if len(ids_peaks) is not 0:
        x_peaks = x[ids_peaks]
        y_peaks = abs(y[ids_peaks])
        g, b0, r_value, p_value, std_err = stats.linregress(
                                            x_peaks - x_peaks[0],
                                            np.log(y_peaks)
                                        )
        out['x_peaks'] = x_peaks
        out['y_peaks'] = y_peaks

        x_fit = x[ids_peaks[0]:ids_peaks[-1] + 1] - x[ids_peaks[0]]
        y_fit = np.exp(b0 + g * x_fit)
        x_fit += x[ids_peaks[0]]
    else:
        out['x_peaks'] = None
        out['y_peaks'] = None

        g, b0, r_value, p_value, std_err = stats.linregress(
                                                    x - x[0],
                                                    np.log(y)
                                                )
        x_fit = x - x[0]
        y_fit = np.exp(b0 + g * x_fit)
        x += x[0]

    # results
    out['g'] = g
    out['y_fit'] = y_fit
    out['x_fit'] = x_fit

    return out


def advanced_g(x_init, y_init, oo={}):
    # additional parameters
    _, ids_x = mix.get_array_oo(oo, x_init, 'x')
    x = np.array(mix.get_slice(x_init, ids_x))
    y = np.array(mix.get_slice(y_init, ids_x))

    # estimation of w,g
    wg_est = oo.get('est', None)
    if wg_est is None:
        wg_est = estimate_wg(x, y)

    # test function
    F = lambda x, A0, g: A0 * np.exp(g*x)

    # initial guess
    p0 = [y[0], wg_est['g']]

    # nonlinear fitting
    popt, pcov = curve_fit(F, x - x[0], y, p0)

    # fitted function
    y_fit = F(x - x[0], *popt)

    # results
    out = {'est': wg_est,
           'g': popt[1],
           'y_fit': y_fit, 'x_fit': x}
    return out


# def w_for_fft(x):
#     nx = np.size(x)
#     nx2_ceil  = np.ceil(np.log2(nx))
#     nx = 2**nx2_ceil
#
#     nx_half = np.int(nx / 2)
#
#     x_new = np.linspace(x[0], x[-1], nx)
#     dx = np.min(np.diff(x_new))
#     freq_max = 1. / dx
#     dfreq = freq_max / nx
#
#     w = np.array([dfreq * i for i in range(nx_half+1)])
#
#     int_shift = 0
#     if np.mod(nx, 2) == 0:
#         int_shift = 1
#     left_a  = -np.flipud(w)
#     right_a = w[2:np.size(w) + 1 - int_shift]
#     w2 = np.concatenate((left_a, right_a))
#
#     return w, w2, nx
#
#
# def prepare_y_for_fft(x, y, nw, axis=0):
#     ny = int(nw)
#     if ny != np.shape(y)[axis]:
#         x_res = np.linspace(x[0], x[-1], ny)
#         f_interp = interpolate.interp1d(x, y, axis=axis)
#         y_res = f_interp(x_res)
#     else:
#         x_res = x
#         y_res = y
#     return y_res, x_res
#
#
# def fft_y_old(y, nw, axis=-1, oo={}):
#     flag_f2_arranged = oo.get('flag_f2_arranged', False)
#
#     ny = int(nw)
#     ny_half = np.int(np.floor(ny/2.))
#
#     f2_raw = np.fft.fft(y, n=ny, axis=axis)  # two-sided FFT
#     f2 = np.abs(f2_raw / ny)
#
#     if axis is -1:
#         f = f2[range(ny_half + 1)]
#         f[1:np.size(f)-1] = 2 * f[1:np.size(f)-1]
#     if axis is 0:
#         f = f2[0:ny_half + 1, :]
#         f[1:np.size(f) - 1, :] = 2 * f[1:np.size(f) - 1, :]
#     if axis is 1:
#         f = f2[:, 0:ny_half + 1]
#         f[:, 1:np.size(f) - 1] = 2 * f[:, 1:np.size(f) - 1]
#
#     if flag_f2_arranged:
#         left_a  = np.flipud(f2[0:ny_half+1])
#         right_a = np.flipud(f2[ny_half+1:np.size(f2) + 1])
#         f2_arranged = np.concatenate((left_a, right_a))
#         return f, f2_raw, f2_arranged
#     else:
#         return f, f2_raw


def fft_y(x, y=None, oo={}):
    # for a 1d or 2d signal y

    # additional parameters:
    flag_f2_arranged = oo.get('flag_f2_arranged', False)
    axis    = oo.get('axis', 0)  # 0, 1

    flag_yw = True
    if y is None:
        flag_yw = False

    # calculate frequency:
    nx2_ceil = np.ceil( np.log2(np.size(x)) )
    nx = int(2 ** nx2_ceil)
    nx_half = np.int(nx / 2)

    x_new = np.linspace(x[0], x[-1], nx)
    dx = np.min(np.diff(x_new))
    freq_max = 1. / dx
    dfreq = freq_max / nx

    w = np.array([dfreq * i for i in range(nx_half + 1)])

    left_a  = - np.flipud(w)
    right_a = w[2:np.size(w)]
    w2 = np.concatenate((left_a, right_a))

    # first part of results
    res = {
        'w': w, 'w2': w2, 'x_new': x_new
    }

    # calculate FFT of the signal y
    if not flag_yw:
        return res

    if nx != np.shape(y)[axis]:
        f_interp = interpolate.interp1d(x, y, axis=axis)
        y_new = f_interp(x_new)
    else:
        x_new, y_new = x, y
    f2_raw = np.fft.fft(y_new, n=nx, axis=axis)  # two-sided FFT
    f2 = np.abs(f2_raw / nx)

    f = None
    if np.size(np.shape(y)) == 1:
        f = f2[range(nx_half + 1)]
        f[1:np.size(f) - 1] = 2 * f[1:np.size(f) - 1]
    elif axis is 0:
        f = f2[0:nx_half + 1, :]
        f[1:np.size(f) - 1, :] = 2 * f[1:np.size(f) - 1, :]
    elif axis is 1:
        f = f2[:, 0:nx_half + 1]
        f[:, 1:np.size(f) - 1] = 2 * f[:, 1:np.size(f) - 1]

    if flag_f2_arranged:
        left_a = np.flipud(f2[0:nx_half + 1])
        right_a = np.flipud(f2[nx_half + 1:np.size(f2) + 1])
        f2_arranged = np.concatenate((left_a, right_a))

        # update results
        res.update({'f2_arranged': f2_arranged})

    # update results:
    res.update({'f': f, 'f2_raw': f2_raw})

    return res


def filtering(x, y, oo):
    def res_None(x_loc, y_loc, norm_w_loc):
        ffres = fft_y(x_loc, y_loc, oo={'flag_f2_arranged': True})
        res = {'x': x,
               'filt': y_loc,
               'fft_init': ffres['f'],
               'fft_filt': None,
               'fft_init_2': ffres['f2_arranged'],
               'fft_filt_2': None,
               'w':  ffres['w']  * norm_w_loc,
               'w2': ffres['w2'] * norm_w_loc}
        return res

    # determine a filter type
    sel_filt = oo.get('sel_filt', 'rough')  # None, 'rough', 'smooth', 'fft_smooth'
    norm_w   = oo.get('norm_w', 1)

    # take the signal in a some predefined interval:
    x_filt, ids_x_filt = mix.get_array_oo(oo, x, 'x')
    y_filt = y[ids_x_filt[0]:ids_x_filt[-1] + 1]

    # no any filtering
    if sel_filt is None:
        out = res_None(x_filt, y_filt, norm_w)
    if sel_filt is 'rough':
        out = rough_filter(x_filt, y_filt, oo)
    if sel_filt is 'smooth':
        out = smooth_filter(x_filt, y_filt, oo)
    if sel_filt is 'fft_smooth':
        out = fft_smooth_filter(x_filt, y_filt, oo)

    out['x'] = x_filt

    return out


def rough_filter(x, y, oo):
    w_interval = oo.get('w_interval', [0, 0])
    norm_w = oo.get('norm_w', 1)

    ffres = fft_y(x, y, oo={'flag_f2_arranged': True})
    w,  w2  = ffres['w'], ffres['w2']

    w = w * norm_w
    w2 = w2 * norm_w

    nw, nw2 = np.size(w), np.size(w2)
    yw2 = ffres['f2_raw']

    # filtering
    if w_interval[-1] is 0:
        yw2[0] = 0.0
    else:
        id_w1, _ = mix.find(w, w_interval[0])
        id_w2, _ = mix.find(w, w_interval[-1])
        yw2[id_w1:id_w2 + 1] = 0.0
        yw2[nw2 - id_w2:nw2 - id_w1] = 0.0

    # new signal and its FFT
    y_new     = np.fft.irfft(yw2, np.size(ffres['x_new']))
    ffres_new = fft_y(ffres['x_new'], y_new, oo={'flag_f2_arranged': True})

    # build the new signal along the initial x-grid
    y_new = np.interp(x, ffres['x_new'], y_new)

    # results
    out = {'filt': y_new,
           'fft_init': ffres['f'],
           'fft_filt': ffres_new['f'],
           'fft_init_2': ffres['f2_arranged'],
           'fft_filt_2': ffres_new['f2_arranged'],
           'w': w, 'w2': w2}

    return out


def smooth_filter(x, y, oo):
    norm_w = oo.get('norm_w', 1)
    wind = int(oo.get('wind', 3))  # windows, which has to be an ODD number

    # smoothing
    y_filt = smooth(y, wind)

    # find FFT
    ffres      = fft_y(x, y,      oo={'flag_f2_arranged': True})
    ffres_filt = fft_y(x, y_filt, oo={'flag_f2_arranged': True})

    # results
    out = {'filt': y_filt,
           'fft_init': ffres['f'],
           'fft_filt': ffres_filt['f'],
           'fft_init_2': ffres['f2_arranged'],
           'fft_filt_2': ffres_filt['f2_arranged'],
           'w': ffres['w'] * norm_w,
           'w2': ffres['w2'] * norm_w
           }
    return out


def fft_smooth_filter(x, y, oo):
    norm_w = oo.get('norm_w', 1)
    w_interval = oo.get('w_interval', [0, 0])
    wind = oo.get('wind', [0, 0])

    ffres = fft_y(x, y, oo={'flag_f2_arranged': True})
    w, w2   = ffres['w'], ffres['w2']
    nw, nw2 = np.size(w), np.size(w2)
    yw2     = ffres['f2_raw']

    w = w * norm_w
    w2 = w2 * norm_w

    # smoothing of the filtering
    id_w1, _ = mix.find(w, w_interval[0])
    id_w2, _ = mix.find(w, w_interval[-1])
    yw2[id_w1:id_w2 + 1] = smooth(yw2[id_w1:id_w2 + 1], wind)
    yw2[nw2 - id_w2:nw2 - id_w1] = smooth(yw2[nw2 - id_w2:nw2 - id_w1], wind)

    # new signal and its FFT
    y_new     = np.fft.irfft(yw2, np.size(ffres['x_new']))
    ffres_new = fft_y(ffres['x_new'], y_new, oo={'flag_f2_arranged': True})

    # build the new signal along the initial x-grid
    if nw != np.size(y):
        y_new = np.interp(x, ffres['x_new'], y_new)

    # results
    out = {'filt': y_new,
           'fft_init': ffres['f'],
           'fft_filt': ffres_new['f'],
           'fft_init_2': ffres['f2_arranged'],
           'fft_filt_2': ffres_new['f2_arranged'],
           'w': w,
           'w2': w2
           }
    return out


def smooth(a, WSZ):
    # taken from
    # https://stackoverflow.com/questions/40443020/matlabs-smooth-implementation-n-point-moving-average-in-numpy-python

    # a: NumPy 1-D array containing the data to be smoothed
    # WSZ: smoothing window size needs, which must be ODD number,
    # as in the original MATLAB implementation
    out0 = np.convolve(a, np.ones(WSZ, dtype=int), 'valid') / WSZ
    r = np.arange(1, WSZ-1, 2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((start, out0, stop))


def find_wc(B0, m, Z):
    # B0 in Tesla, m in kg, T in J
    wc = Z * constants.elementary_charge * B0 / m  # in rad/s
    return wc


def find_Te(Lx, a, B0, mf, Zf):
    # find temperature Ti of a species i:
    # tau_i = Ti / Te
    # mf - mass of the first species
    # Zf - charge (without e) of the first species
    # Remark: in ORB5 tau_e = 1 always
    wc = find_wc(B0, mf, Zf)
    Te = (2*a*wc/Lx)**2 * mf  # in J
    return Te


def find_Ti(tau_i, Lx, a, B0, mf, Zf):
    # find electron temperature Te:
    # mf - mass of the first species
    # Zf - charge (without e) of the first species
    # Remark: in ORB5 tau_e = 1 always
    Te = find_Te(Lx, a, B0, mf, Zf)
    Ti = Te * tau_i
    return Ti


def find_Lx(a, Te, B0, mf, Zf):
    # mf - mass of the first species (kg)
    # Zf - charge (without e) of the first species
    # Te - electron temperature (J)
    # a - minor radius (m)
    # B0 - magnetic field (T)
    wc = find_wc(B0, mf, Zf)
    cs = find_cs(Te, mf)
    Lx = 2*a*wc / cs
    return Lx


def find_vt(T, m):
    # m in kg, T in J
    vt = np.sqrt(T/m)
    return vt


def find_cs(Te, m):
    # m in kg, Te in J
    cs = np.sqrt(Te/m)
    return cs


def find_k(krpert, Tf, B0, Lwork, mf, Zf):
    # Find normalised radial wavenumber
    # mf - mass of the first species (kg)
    # Zf - charge (without e) of the first species
    # Tf - temperature of the first species (J)
    # B0 - magnetic field (Tesla)
    # krpert - radial vector from an input file
    rhoL = find_rhoL(Tf, B0, mf, Zf)
    k = krpert * (np.pi / Lwork) * rhoL
    return k


def find_rhoL(T, B0, m, Z):
    # Find a value of the Larmor radius
    # m - mass of a species (kg)
    # Z - charge (without e) of a species
    # T - temperature of a species (J)
    # B0 - magnetic field (Tesla)
    wc = find_wc(B0, m, Z)
    vt = find_vt(T, m)
    rhoLarmor = np.sqrt(2) * vt / wc
    return rhoLarmor


















