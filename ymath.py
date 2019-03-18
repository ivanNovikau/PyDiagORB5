import Mix as mix
import numpy as np
from scipy.signal import find_peaks
from scipy import stats

#-> avoid plotting here

def reload():
    # Important: put here all modules that you want to reload
    mix.reload_module(mix)

def estimate_wg(x, y):
    # find peaks of the signal
    ids_peaks, _ = find_peaks(abs(y), height=0)
    x_peaks = x[ids_peaks]
    y_peaks = abs(y[ids_peaks])

    # frequency
    w = 2 * np.pi * 0.5 / np.mean(x_peaks)
    g, b0,  r_value, p_value, std_err = stats.linregress(x_peaks, np.log(y_peaks))
    y_fit = np.exp(b0 + g * x)

    # results
    out={}
    out['x_peaks'] = x_peaks
    out['y_peaks'] = y_peaks
    out['w'] = w
    out['g'] = g
    out['y_fit'] = y_fit
    return out


def advanced_wg(x, y, oo={}):
    













