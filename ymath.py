import Mix as mix
import numpy as np
import scipy.signal
from scipy import stats
from scipy.optimize import curve_fit

#-> avoid plotting here


def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)


def estimate_wg(x_init, y_init, oo={}):
    # additional parameters
    _, ids_x = mix.get_array_oo(oo, x_init, 'x')
    x = np.array(mix.get_slice(x_init, ids_x))
    y = np.array(mix.get_slice(y_init, ids_x))

    # find peaks of the signal
    ids_peaks, _ = scipy.signal.find_peaks(abs(y), height=0)
    x_peaks = x[ids_peaks]
    y_peaks = abs(y[ids_peaks])

    # frequency
    w = 2 * np.pi * 0.5 / np.mean(np.diff(x_peaks))
    g, b0,  r_value, p_value, std_err = stats.linregress(x_peaks, np.log(y_peaks))

    x_fit = x
    y_fit = np.cos(w * x_fit) * np.exp(b0 + g * x_fit)

    # results
    out = {'x_peaks': x_peaks, 'y_peaks': y_peaks,
          'w': w, 'g': g,
          'y_fit': y_fit, 'x_fit': x_fit}
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
    F = lambda x, A0, w, g: A0 * np.cos(w*x) * np.exp(g*x)

    # initial guess
    p0 = [y[0], wg_est['w'], wg_est['g']]

    # nonlinear fitting
    popt, pcov = curve_fit(F, x, y, p0)

    # fitted function
    y_fit = F(x, *popt)

    # results
    x_fit = x
    out = {'est': wg_est,
           'w': popt[1], 'g': popt[2],
           'y_fit': y_fit, 'x_fit': x_fit}
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
        g, b0, r_value, p_value, std_err = \
            stats.linregress(x_peaks, np.log(y_peaks))
        out['x_peaks'] = x_peaks
        out['y_peaks'] = y_peaks
    else:
        g, b0, r_value, p_value, std_err = \
            stats.linregress(x, np.log(y))

    x_fit = x
    y_fit = np.exp(b0 + g * x_fit)

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
    popt, pcov = curve_fit(F, x, y, p0)

    # fitted function
    y_fit = F(x, *popt)

    # results
    x_fit = x
    out = {'est': wg_est,
           'g': popt[1],
           'y_fit': y_fit, 'x_fit': x_fit}
    return out


def w_for_fft(x):
    nx = np.size(x)
    nx_half = np.int(np.floor(nx / 2.))
    dx = np.min(np.diff(x))
    freq_max = 1. / dx
    dfreq = freq_max / nx

    w = np.array([dfreq * i for i in range(nx_half+1)])
    w2 = []
    return w, w2


def fft_y(y, axis=-1):
    if axis == -1:
        ny = np.size(y)
    else:
        n_shape  = np.shape(y)
        ny = n_shape[axis]
    ny_half = np.int(np.floor(ny/2.))
    f2 = np.fft.fft(y, axis=axis)  # two-sided FFT
    f2 = np.abs(f2 / ny)

    if axis is -1:
        f = f2[range(ny_half + 1)]
        f[1:np.size(f)-1] = 2 * f[1:np.size(f)-1]
    if axis is 0:
        f = f2[0:ny_half + 1, :]
        f[1:np.size(f) - 1, :] = 2 * f[1:np.size(f) - 1, :]
    if axis is 1:
        f = f2[:, 0:ny_half + 1]
        f[:, 1:np.size(f) - 1] = 2 * f[:, 1:np.size(f) - 1]

    f2 = []

    return f, f2
















