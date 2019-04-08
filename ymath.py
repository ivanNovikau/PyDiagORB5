import Mix as mix
import numpy as np
import scipy.signal
from scipy import stats
from scipy.optimize import curve_fit
from scipy import constants
from scipy.fftpack import fft as sci_fft

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
    nx2_ceil  = np.ceil(np.log2(nx))
    nx = 2**nx2_ceil

    nx_half = np.int(np.floor(nx / 2.))
    dx = np.min(np.diff(x))
    freq_max = 1. / dx
    dfreq = freq_max / nx

    w = np.array([dfreq * i for i in range(nx_half+1)])

    w2 = []
    return w, w2, nx


def fft_y(y, nw, axis=-1):
    ny = int(nw)
    ny_half = np.int(np.floor(ny/2.))

    f2 = np.fft.fft(y, n=ny, axis=axis)  # two-sided FFT
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


















