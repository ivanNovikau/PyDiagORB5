from IPython.core.getipython import get_ipython
import numpy as np


def reload():
    return


# ---------------------------------------------------------------------------
MIN_N_PEAKS = 3
COEF_ERR = 1.96          # alpha-quantile of the standard normal distribution
CONFIDENCE_PERC = 0.95   # (2*alpha - 1)-confidence interval,
                         # that corresponds to the alpha quantile

# ---------------------------------------------------------------------------
if 'Terminal' in get_ipython().__class__.__name__:
    FLAG_LATEX = True
else:
    FLAG_LATEX = False

DASHES_FORMAT = [0.5, 0.4]
MARKER_EDGE_WIDTH_COEF = 0.5
ERR_WIDTH_COEF = 0.33
DEF_MAXLOCATOR = 6
COLORMAP_LEVELS = 10
DEF_ONE_COLOR = 'blue'
DEF_ONE_STYLE = '-'
DEF_COLORMAP = 'jet'
if FLAG_LATEX:
    FLAG_LATEX = True
    FIG_SIZE_W = 15
    FIG_SIZE_H = 9.5
    LEG_SCALE = 1.0
    FONT_SIZE = 28
    SCALE_LABELS = 1.8
    SCALE_TICKS  = 1.5
    SCALE_ORDER  = 1.2
    SCALE_TITLE  = 1.3
    LINE_WIDTH = 6
    MARKER_SIZE = 14
else:
    FLAG_LATEX = False
    FIG_SIZE_W = 10
    FIG_SIZE_H = 6
    LEG_SCALE = 0.5
    FONT_SIZE = 22
    SCALE_LABELS = 0.6
    SCALE_TICKS = 0.6
    SCALE_ORDER = 0.6
    SCALE_TITLE = 0.6
    LINE_WIDTH = 6
    MARKER_SIZE = 14

# ---------------------------------------------------------------------------
NONE_FILTER = {'sel_filt': None}
DEF_COLORS = ['b', 'r', 'g', 'black', 'm', 'c',
              'y', 'k', 'cyan', 'Purple', 'gray', 'lightcoral']
DEF_STYLES = ['-', ':', '-.', ':']

DEF_PLOT_FORMAT = {  # describe format of a plot
    'xlabel': None,
    'ylabel': None,
    'zlabel': None,
    'wlabel': None,
    'title': None,
    'flag_semilogy': False,
    'flag_norm': False,
    'flag_colorbar': True,
    'fontS': FONT_SIZE,
    'xlimits': None,
    'ylimits': None,
    'zlimits': None,
    'xticks': np.nan,
    'yticks': np.nan,
    'xticks_labels': np.nan,
    'yticks_labels': np.nan,
    'flag_legend': True,
    'legend_position': 'best',  # 'upper right', 'center left'
    'legend_fcol': 'lightgray',
    'flag_diff_styles': False,
    'x_style': 'sci',  # 'sci', 'plain'
    'y_style': 'sci',  # 'sci', 'plain'
    'flag_maxlocator': False,
    'maxlocator': DEF_MAXLOCATOR,
}

DEF_CURVE_FORMAT = {  # describe format of a curve
    'legend': None,
    'style': '-',
    'width': LINE_WIDTH,
    'color': DEF_ONE_COLOR,
    'markersize': MARKER_SIZE,
    'markerfacecolor': "None",
    'colormap': DEF_COLORMAP,  # hot, jet, pink, hot_r, jet_r etc.
    'colormap_center': None,
    'levels': COLORMAP_LEVELS,  # for contour plot
    'pr_alpha': 1,
    'flag_errorbar': False,
    'flag_hist': False,
}


def new_color(count_color):
    if count_color < len(DEF_COLORS):
        one_color = DEF_COLORS[count_color]
    else:
        one_color = ','.join('{}'.format(*k) for k in enumerate(np.random.rand(3, 1)))
        one_color = 'rgb({})'.format(one_color)
    return one_color


def new_style(count_style):
    if count_style < len(DEF_STYLES):
        one_style = DEF_STYLES[count_style]
    else:
        one_style = ':'
    return one_style


# ---------------------------------------------------------------------------
# --- DEFAULT VARIABLE DEFINITIONS ---
# ---------------------------------------------------------------------------
def_erbar_ts = {
    'type':             'zonal',
    'variable':         'er',
    'plane':            'ts',
    'avr_operation':    'point-s',
    'avr_domain':       0.5,
}
def_safety_factor = {
    'type': 'equ-profile',
    'variable': 'q',
    'plane': 'ts',
    'avr_operation': 'point-t',
    'avr_domain': 0,
}
def_Teq_deuterium = {
    'type': 'equ-profile',
    'variable': 'T-equ',
    'species_name': 'deuterium',
    'plane': 'ts',
    'avr_operation': 'point-t',
    'avr_domain': 0,
}


def create_signal(default_signal, dd):
    res_signal = dict(default_signal)
    res_signal.update({
        'dd': dd
    })
    return [res_signal]