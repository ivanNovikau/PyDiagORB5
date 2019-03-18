import Mix as mix
import numpy as np
from scipy.signal import find_peaks
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
    ids_peaks, _ = find_peaks(abs(y), height=0)
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
    # x_to_fit = x - x[0]
    x_to_fit = x
    popt, pcov = curve_fit(F, x_to_fit, y, p0)

    # fitted function
    y_fit = F(x_to_fit, *popt)

    # results
    x_fit = x
    out = {'est': wg_est,
           'w': popt[1], 'g': popt[2],
           'y_fit': y_fit, 'x_fit': x_fit}
    return out












